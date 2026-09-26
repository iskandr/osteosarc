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

    def __init__(self, rows, *, columns=None, source=None, diagnostics=()):
        self.rows = tuple(rows)
        self.columns = tuple(columns if columns is not None else
                             dict.fromkeys(k for row in self.rows for k in row))
        self.source = source
        self.diagnostics = tuple(diagnostics)

    def __len__(self):
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def where(self, predicate):
        selected = [i for i, row in enumerate(self) if predicate(row)]
        positions = {old: new for new, old in enumerate(selected)}
        diagnostics = [dict(d, row=positions[d["row"]]) for d in self.diagnostics if d["row"] in positions]
        return Table((self.rows[i] for i in selected), columns=self.columns, source=self.source,
                     diagnostics=diagnostics)

    def select(self, **fields):
        unknown = set(fields) - set(self.columns)
        if unknown:
            raise KeyError(f"Unknown columns: {sorted(unknown)}")
        return self.where(lambda row: all(row.get(k) == v for k, v in fields.items()))

    def __repr__(self):
        from .display import preview
        shown = self.columns[:6]
        note = f" (first {len(shown)} shown)" if len(self.columns) > len(shown) else ""
        return preview(f"Table: {len(self):,} rows, {len(self.columns)} columns{note}",
                       self.rows, shown, lambda row: row)

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


def parse_table(text, *, delimiter="\t", source=None, required=(), strict=True):
    """Parse CSV/TSV; strict=False retains ragged rows and their diagnostics.

    Missing trailing fields are None; extra fields are retained in diagnostics.
    Every diagnostic has a zero-based row index and the original parsed fields.
    Missing/duplicate headers are errors in either mode.
    """
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    columns = reader.fieldnames
    if not columns or len(set(columns)) != len(columns) or not set(required) <= set(columns):
        raise SchemaError(f"Missing or duplicate table columns; required={tuple(required)}")
    rows, diagnostics = [], []
    for row in reader:
        if None in row or any(value is None for value in row.values()):
            if strict:
                raise SchemaError(f"Ragged table row {reader.line_num}")
            fields = tuple(row[c] for c in columns if row[c] is not None) + tuple(row.pop(None, ()))
            diagnostics.append(dict(row=len(rows), line=reader.line_num, code="ragged_row",
                                    message=f"Expected {len(columns)} fields, found {len(fields)}",
                                    fields=fields))
        rows.append(row)
    return Table(rows, columns=columns, source=source, diagnostics=diagnostics)


def parse_file(path, *, format=None, source=None):
    """Parse JSON, CSV, TSV, or FASTA without changing source annotations."""
    suffix = Path(str(path).removesuffix(".gz")).suffix.lstrip(".").lower()
    rsem = str(path).removesuffix(".gz").endswith((".genes.results", ".isoforms.results"))
    format = format or ("tsv" if rsem else suffix)
    if format not in PARSE_FORMATS:
        raise ValueError(f"No built-in parser for {format!r}; download the original file")
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
    raise ValueError(f"No built-in parser for {format!r}; download the original file")


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


def parse_variants(index, vafs, *, source_variants=(), vaccine_overlap=None, source=None):
    """Join exact IDs; preserve unresolvable entries and source disagreements.

    index is the variant index HTML or its already parsed rows. The GRCh38
    label is the site's assertion, not independent REF validation. This does
    not normalize indels, infer alleles from protein names, or lift. vafs may be
    TSV text or a Table. Malformed rows with an ID mark that entry non-ready;
    its annotations["parse_errors"] retain the diagnostics and source values.
    A reviewed GRCh38 allele_resolution can fill an allele absent from vafs.
    """
    rows = parse_variant_index(index) if isinstance(index, str) else index
    entries = {row["id"]: dict(row) for row in rows}
    required = {"variant_id", "gene", "chrom", "pos", "ref", "alt"}
    if isinstance(vafs, str):
        vafs = parse_table(vafs, required=required, strict=False)
    if len(set(vafs.columns)) != len(vafs.columns) or not required <= set(vafs.columns):
        raise SchemaError(f"Missing or duplicate VAF table columns; required={sorted(required)}")
    alleles, annotations, errors = defaultdict(set), {}, defaultdict(list)
    row_errors = defaultdict(list)
    for diagnostic in getattr(vafs, "diagnostics", ()):
        row_errors[diagnostic["row"]].append(diagnostic)
    for i, row in enumerate(vafs):
        vid = row.get("variant_id")
        if not isinstance(vid, str) or not vid.strip():
            raise SchemaError(f"VAF row {i + 1} has no recoverable variant ID")
        if vid not in entries:
            entries[vid] = dict(id=vid, gene=row.get("gene"), on_site=False)
        problems = list(row_errors[i])
        if not problems and any(row.get(field) is None for field in required):
            problems.append(dict(row=i, code="missing_field", message="Missing required VAF field"))
        if not problems and (not re.fullmatch(r"[0-9]+", str(row["pos"])) or int(row["pos"]) < 1):
            problems.append(dict(row=i, code="invalid_position", message=f"Invalid variant position: {row['pos']!r}"))
        if problems:
            errors[vid].extend(dict(p, source="vafs", values=dict(row)) for p in problems)
            continue
        allele = (row["chrom"], int(row["pos"]), row["ref"], row["alt"])
        alleles[vid].add(allele)
        annotations.setdefault(vid, {k: row.get(k) for k in
                                     ("consequence", "protein_change", "variant_type")})
    source_by_id = {}
    for record in source_variants:
        if record["id"] in source_by_id:
            raise SchemaError(f"Duplicate source variant ID: {record['id']}")
        source_by_id[record["id"]] = record
        entries.setdefault(record["id"], dict(id=record["id"], gene=record["gene"], on_site=False))
        # Reviewed source corrections can supply an allele even when the site
        # never counted it. Never replace VAF candidates or fabricate count rows.
        resolution = record.get("allele_resolution", {})
        if not alleles[record["id"]] and resolution.get("status") == "resolved":
            resolved = resolution["allele"]
            allele = (record.get("chr"), record.get("pos"), record.get("ref"), record.get("alt"))
            if (resolved["assembly"] == "GRCh38"
                    and allele == tuple(resolved[k] for k in ("chrom", "pos", "ref", "alt"))):
                alleles[record["id"]].add(allele)
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
        if "allele_resolution" in record:
            extra["allele_resolution"] = record["allele_resolution"]
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
        if errors[vid]:
            status = "malformed_source_row"
            extra["parse_errors"] = errors[vid]
        result.append(Variant(vid, entry["gene"], "GRCh38", candidates, status,
                              entry.get("on_site", True), count, membership,
                              tuple(sorted(k for k, value in record.get("detection", {}).items() if value)),
                              extra))
    return Variants(result, source=source)
