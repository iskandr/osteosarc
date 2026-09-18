"""Offline parsers. Retain source rows; normalize only documented identity fields."""

from __future__ import annotations

import csv
import gzip
import io
import json
import re
from collections import defaultdict
from pathlib import Path

from bs4 import BeautifulSoup

from .errors import SchemaError
from .models import Variant, Variants

PARSE_FORMATS = frozenset(("json", "csv", "tsv", "fasta", "fa", "fna", "faa"))


class Table:
    """Source-preserving table. Missing strings and numeric zero remain distinct."""

    def __init__(self, rows, *, columns=None, source=None):
        self.rows = tuple(rows)
        self.columns = tuple(columns if columns is not None else
                             dict.fromkeys(k for row in self.rows for k in row))
        self.source = source

    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def where(self, predicate):
        return Table((r for r in self if predicate(r)), columns=self.columns, source=self.source)

    def select(self, **fields):
        unknown = set(fields) - set(self.columns)
        if unknown:
            raise KeyError(f"Unknown columns: {sorted(unknown)}")
        return self.where(lambda row: all(row.get(k) == v for k, v in fields.items()))

    def to_dataframe(self):
        """Optional pandas conversion; performs no type or missing-value coercion."""
        import pandas as pd
        return pd.DataFrame(self.rows, columns=self.columns)


def read_text(path):
    """Read UTF-8 text, recognizing gzip by magic bytes, not a cache filename."""
    path = Path(path)
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return handle.read()


def parse_table(text, *, delimiter="\t", source=None, required=()):
    """Parse CSV/TSV with a header; preserve raw strings and reject ragged rows."""
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    columns = reader.fieldnames
    if not columns or len(set(columns)) != len(columns) or not set(required) <= set(columns):
        raise SchemaError(f"Missing or duplicate table columns; required={tuple(required)}")
    rows = []
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            raise SchemaError(f"Ragged table row {reader.line_num}")
        rows.append(row)
    return Table(rows, columns=columns, source=source)


def parse_file(path, *, format=None, source=None):
    """Parse JSON, CSV, TSV, or FASTA without changing source annotations."""
    suffix = Path(str(path).removesuffix(".gz")).suffix.lstrip(".").lower()
    rsem = str(path).removesuffix(".gz").endswith((".genes.results", ".isoforms.results"))
    format = format or ("tsv" if rsem else suffix)
    if format not in PARSE_FORMATS:
        raise ValueError(f"No built-in parser for {format!r}; download the original asset")
    text = read_text(path)
    if format == "json":
        return json.loads(text)
    if format in ("csv", "tsv"):
        return parse_table(text, delimiter="," if format == "csv" else "\t", source=source)
    if format in ("fasta", "fa", "fna", "faa"):
        rows, name, sequence = [], None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if name is not None:
                    rows.append(dict(name=name, sequence="".join(sequence)))
                name, sequence = line[1:], []
            elif line.strip():
                if name is None:
                    raise SchemaError("FASTA sequence precedes its header")
                sequence.append(line.strip())
        if name is not None:
            rows.append(dict(name=name, sequence="".join(sequence)))
        return Table(rows, columns=("name", "sequence"), source=source)
    raise ValueError(f"No built-in parser for {format!r}; download the original asset")


def parse_variant_index(html):
    """Read website identities and counts; fail on an incompatible table schema."""
    soup = BeautifulSoup(html, "html.parser")
    rows = []
    for row in soup.select("tr[data-vaccines]"):
        cells = row.find_all("td", recursive=False)
        link = row.select_one('a[href^="/variant/"]')
        if len(cells) != 9 or link is None:
            raise SchemaError("Variant index schema changed")
        rows.append(dict(id=link["href"].rstrip("/").split("/")[-1],
                         gene=cells[0].get_text(" ", strip=True),
                         location=cells[1].get_text(" ", strip=True),
                         protein_label=cells[2].get_text(" ", strip=True),
                         vaccine_count=int(row["data-vaccines"])))
    if not rows or len({r["id"] for r in rows}) != len(rows):
        raise SchemaError("Missing or duplicate variant IDs")
    return rows


def parse_variants(index_html, vafs, *, source_variants=(), vaccine_overlap=None, source=None):
    """Join exact IDs; preserve unresolvable entries and source disagreements.

    The GRCh38 label is the site's assertion, not independent REF validation.
    This does not normalize indels, infer alleles from protein names, or lift.
    """
    entries = {row["id"]: row for row in parse_variant_index(index_html)}
    required = {"variant_id", "gene", "chrom", "pos", "ref", "alt"}
    if not required <= set(vafs.columns):
        raise SchemaError(f"VAF table missing fields: {required - set(vafs.columns)}")
    alleles, annotations = defaultdict(set), {}
    for row in vafs:
        vid = row["variant_id"]
        try:
            allele = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
        except ValueError as error:
            raise SchemaError(f"Invalid variant position: {vid}") from error
        alleles[vid].add(allele)
        annotations.setdefault(vid, {k: row.get(k) for k in
                                     ("consequence", "protein_change", "variant_type")})
        if vid not in entries:
            entries[vid] = dict(id=vid, gene=row["gene"], on_site=False)
    source_by_id = {}
    for record in source_variants:
        if record["id"] in source_by_id:
            raise SchemaError(f"Duplicate source variant ID: {record['id']}")
        source_by_id[record["id"]] = record
        entries.setdefault(record["id"], dict(id=record["id"], gene=record["gene"], on_site=False))
    overlap = defaultdict(list)
    for record in (vaccine_overlap or {}).get("mutations", []):
        if record.get("chrom") and record.get("pos"):
            overlap[(record["gene"], record["chrom"], int(record["pos"]))].append(record)
    result = []
    locus_ids = defaultdict(set)
    for vid, candidates in alleles.items():
        if len(candidates) == 1:
            chrom, pos, _, _ = next(iter(candidates))
            locus_ids[(entries[vid]["gene"], chrom, pos)].add(vid)
    for vid, entry in entries.items():
        candidates = tuple(sorted(alleles[vid]))
        status = "missing_literal_allele" if not candidates else (
            "ambiguous_literal_allele" if len(candidates) > 1 else "ready")
        record = source_by_id.get(vid, {})
        extra = dict(annotations.get(vid, {}), index=entry, source_record=record)
        source_membership = tuple(sorted(k for k, value in record.get("vaccines", {}).items() if value))
        extra["source_vaccines"] = source_membership
        membership = ()
        if len(candidates) == 1:
            chrom, pos, ref, alt = candidates[0]
            if pos < 1 or not re.fullmatch("[ACGT]+", ref) or not re.fullmatch("[ACGT]+", alt):
                status = "non_literal_allele"
            elif entry.get("location", f"{chrom}:{pos}") != f"{chrom}:{pos}":
                status = "conflicting_coordinates"
            # A source JSON allele does not silently replace an exported allele.
            if record.get("ref") and record.get("alt") and record.get("pos"):
                original = (record.get("chr"), int(record["pos"]), record["ref"], record["alt"])
                if original != candidates[0]:
                    extra["source_allele_conflict"] = original
                    status = "conflicting_coordinates"
            matches = overlap.get((entry["gene"], chrom, pos), [])
            extra["vaccine_overlap_records"] = matches
            if len(matches) == 1 and len(locus_ids[(entry["gene"], chrom, pos)]) == 1:
                overlap_names = tuple(sorted(k for k, value in matches[0]["vaccines"].items() if value))
                membership = overlap_names
                if record and source_membership != overlap_names:
                    extra["vaccine_membership_conflict"] = dict(source=source_membership, overlap=overlap_names)
            elif matches:
                extra["ambiguous_vaccine_join"] = True
        count = entry.get("vaccine_count")
        if count is None and record:
            count = len(source_membership)
        result.append(Variant(vid, entry["gene"], "GRCh38", candidates, status,
                              entry.get("on_site", True), count, membership,
                              tuple(sorted(k for k, value in record.get("detection", {}).items() if value)),
                              extra))
    return Variants(result, source=source)
