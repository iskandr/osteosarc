"""Reconcile all bucket objects with explicit, source-attributed metadata."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from pathlib import PurePosixPath
from urllib.parse import quote, unquote, urljoin, urlsplit

from .cache import stable_id
from .errors import SchemaError
from .models import Asset, Assets, SampleClaim

SITE = "https://osteosarc.com/"
BUCKET = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
SOURCE_REPO = "https://gitlab.com/slowkow/osteosarc.com/-/raw/main/"

SNAPSHOT_SOURCES = {
    "bams": SITE + "bams/bams.json",
    "bucket": SITE + "bucket_listing.json",
    "variant_index": SITE + "variants/",
    "vafs": SITE + "variants/variant_vafs_long.tsv",
    "vaf_columns": SITE + "variants/variant_vafs_long.columns.tsv",
    "vaccine_overlap": SITE + "data/vaccine_overlap.json",
    "source_variants": SOURCE_REPO + "src/data/variants.json",
    "bam_metadata": SOURCE_REPO + "scripts/data/bam-metadata-consolidated.tsv",
    "data_page": SITE + "data/",
}

TABLE_SOURCES = {
    "vafs": (SNAPSHOT_SOURCES["vafs"], "tsv"),
    "vaf_columns": (SNAPSHOT_SOURCES["vaf_columns"], "tsv"),
    "snv_top": (SITE + "oncoanalyser/tables/snv_top.tsv", "tsv"),
    "dna_fusions": (SITE + "fusions/tables/prioritized.tsv", "tsv"),
    "rna_fusions": (SITE + "ctat_lr_fusion/tables/lr_fusions_concordance.tsv", "tsv"),
}

ASSAYS = {"RNA": ("rna-seq", None), "WGS": ("wgs", None), "WES": ("wes", None),
          "scRNA ONT": ("scrna-seq", "ont"), "scRNA_ONT": ("scrna-seq", "ont"),
          "PacBio": ("scrna-seq", "pacbio"), "scRNA": ("scrna-seq", None),
          "Tumor scRNA": ("scrna-seq", None), "Blood scRNA": ("scrna-seq", None),
          "CITE": ("cite-seq", None)}


def normalize_provider(value):
    """Normalize spelling aliases without changing the underlying source rows."""
    return {"Boston Gene": "BostonGene"}.get(value, value) or None


def normalize_timepoint(value):
    """Separate the explicit T1-organoid label from its specimen annotation."""
    return {"T1-organoid": "T1"}.get(value, value) or None


def asset_type(key):
    """Classify file format separately from scientific interpretation."""
    lower = key.lower()
    if lower.endswith((".bai", ".csi", ".crai", ".tbi", ".fai")):
        return "index", lower.rsplit(".", 1)[-1]
    base = lower.removesuffix(".gz").removesuffix(".bgz")
    if base.endswith((".genes.results", ".isoforms.results")):
        return "expression", "tsv"
    suffix = base.rsplit(".", 1)[-1] if "." in base else ""
    if suffix in ("bam", "cram", "sam"):
        return "alignment", suffix
    if suffix in ("vcf", "bcf"):
        return "variants", suffix
    if suffix in ("fastq", "fq"):
        return "reads", "fastq"
    if suffix in ("fasta", "fa", "fna", "faa"):
        return "reference", "fasta"
    if suffix in ("gtf", "gff", "gff3", "bed"):
        return "annotation", suffix
    if suffix in ("tsv", "csv", "json"):
        return "table", suffix
    if suffix in ("h5", "h5ad", "h5mu", "rds", "mtx", "cloupe"):
        return "expression", suffix
    return "other", suffix


def bucket_url(key, base=BUCKET):
    """Encode an exact S3 object key once; '+' and spaces remain distinct."""
    return base.rstrip("/") + "/" + quote(key, safe="/")


def _key(value, base):
    if value.startswith(base):
        return unquote(value[len(base):])
    if urlsplit(value).scheme:
        raise SchemaError(f"Object is outside the catalog bucket: {value}")
    return value


def build_assets(listing, bams, metadata, vafs, path_claims=()):
    """Retain every listed object, enriching exact paths before basename matches.

    Basename joins are used only when unique among alignment objects. Path
    inferences are explicitly marked and excluded from default metadata filters.
    The global viewer genome is retained as a claim, not assigned as assembly.
    """
    base = listing.get("download_base", BUCKET)
    objects = {}
    for row in listing["files"]:
        if len(row) < 3 or row[0] in objects:
            raise SchemaError("Duplicate or malformed bucket object")
        objects[row[0]] = dict(size=int(row[1]), modified=row[2])
    catalog, claims, extra = {}, defaultdict(list), defaultdict(dict)
    for category in bams["categories"]:
        for row in category["bams"]:
            key = _key(row["url"], bams["baseUrl"])
            if key in catalog:
                raise SchemaError(f"Duplicate catalog alignment: {key}")
            catalog[key] = row
            objects.setdefault(key, dict(size=None, modified=None))
            assay, platform = ASSAYS.get(category["name"], (None, None))
            match = re.match(r"(T\d+)\b", row["name"])
            claims[key].append(SampleClaim("bams", row["name"], match[1] if match else None,
                                           assay=assay, platform=platform, tissue=row.get("tissue")))
            extra[key].update(catalog=row, category=category["name"],
                              catalog_genome_assertion=bams.get("genome"))
    for row in metadata:
        if not row.get("s3_path"):
            continue
        key = _key(row["s3_path"], base)
        objects.setdefault(key, dict(size=None, modified=None))
        assay, platform = ASSAYS.get(row.get("assay"), (None, None))
        claims[key].append(SampleClaim("bam_metadata", row["display_name"], normalize_timepoint(row.get("timepoint")),
                                       row.get("sample_date") or None, assay, platform,
                                       row.get("tissue", "").lower() or None, normalize_provider(row.get("provider"))))
        extra[key].setdefault("metadata_rows", []).append(row)
    counts = Counter(PurePosixPath(key).name for key in objects if asset_type(key)[0] == "alignment")
    by_basename = defaultdict(set)
    for row in vafs:
        by_basename[row["bam_file"]].add(tuple(row.get(k, "") for k in
            ("sample_label", "timepoint", "sample_date", "assay_type", "tissue", "data_source")))
    assets = []
    for key, object_metadata in objects.items():
        kind, format = asset_type(key)
        info = dict(extra.get(key, {}))
        records = list(claims.get(key, ()))
        for prefix, claim in path_claims:
            if key == prefix.rstrip("/") or key.startswith(prefix.rstrip("/") + "/"):
                records.append(claim)
        if kind == "alignment":
            basename = PurePosixPath(key).name
            if counts[basename] == 1:
                for label, timepoint, date, assay_name, tissue, provider in sorted(by_basename[basename]):
                    assay, platform = ASSAYS.get(assay_name, (None, None))
                    records.append(SampleClaim("vafs", label, normalize_timepoint(timepoint), date or None,
                                               assay, platform, tissue or None, normalize_provider(provider)))
            elif by_basename[basename]:
                info["ambiguous_vaf_basename"] = basename
        # This inference is useful for discovery, but never establishes identity.
        points = sorted(set(re.findall(r"(?:^|[/_ .-])(T[0-3])(?=[/_ .-]|$)", key)))
        libraries = sorted(set(re.findall(r"\b(?:BG\d{6}|SARC\d{4}|TL-\d{2}-[A-Z0-9]+)\b", key)))
        if kind != "other":
            for point in points:
                records.append(SampleClaim("bucket_path", key, timepoint=point, basis="inferred"))
            for library in libraries:
                records.append(SampleClaim("bucket_path", key, library=library, basis="inferred"))
        suffixes = [key + ".bai", key[:-4] + ".bai", key + ".csi"] if format == "bam" else (
            [key + ".crai", key[:-5] + ".crai", key + ".csi"] if format == "cram" else
            [key + ".tbi", key + ".csi"] if format in ("vcf", "bcf") else [])
        indexes = tuple(bucket_url(k, base) for k in suffixes if k in objects)
        url = bucket_url(key, base)
        assets.append(Asset(stable_id(url), key, url, kind, format, index_urls=indexes,
                            claims=tuple(records), metadata=info, **object_metadata))
    for name, (url, format) in TABLE_SOURCES.items():
        assets.append(Asset(stable_id(url), "site/" + name, url, "table", format,
                            metadata={"resource": name}))
    return Assets(assets)


def parse_data_paths(html):
    """Extract the data page's explicit path-to-sample claims, including FASTQs.

    Assays/platforms come from section headings, tissues/providers from cells.
    Composite tissues and timepoints remain unresolved rather than guessed.
    """
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(html, "html.parser")
    headings, result = {}, []
    for element in soup.find_all(["h2", "h3", "h4", "table"]):
        if element.name != "table":
            level = int(element.name[1])
            headings = {k: v for k, v in headings.items() if k < level}
            headings[level] = element.get_text(" ", strip=True)
            continue
        context = " / ".join(headings.values())
        title = context.lower()
        assay = ("wgs" if "whole genome" in title or re.search(r"\bwgs\b", title) else
                 "wes" if "whole exome" in title or re.search(r"\bwes\b", title) else
                 "scrna-seq" if any(x in title for x in ("single cell", "single-cell")) else
                 "rna-seq" if "bulk rna" in title else None)
        platform = ("ont" if "nanopore" in title else "pacbio" if "pacbio" in title else
                    "illumina" if "illumina" in title else None)
        headers = [c.get_text(" ", strip=True) for c in element.select("thead th")]
        if not headers:
            headers = [c.get_text(" ", strip=True) for c in element.select("tr:first-child th")]
        for row in element.select("tr"):
            cells = row.find_all("td", recursive=False)
            if len(cells) != len(headers):
                continue
            values = dict(zip(headers, (c.get_text(" ", strip=True) for c in cells)))
            if "Bucket Path" not in values:
                continue
            points = set(re.findall(r"\bT[0-3]\b", values.get("Timepoint", "")))
            tissue = values.get("Tissue", "").lower()
            tissue = {"tumor": "tumor", "normal": "normal", "normal (blood)": "blood",
                      "blood normal": "blood", "organoid": "organoid"}.get(tissue)
            label = " | ".join(v for k, v in values.items() if k not in ("Bucket Path", "Size", "Files"))
            result.append((values["Bucket Path"].strip("`"), SampleClaim(
                "data_page", context + " / " + label,
                timepoint=next(iter(points)) if len(points) == 1 else None,
                assay=assay, platform=platform, tissue=tissue, provider=normalize_provider(values.get("Provider")))))
    return tuple(result)


def table_links(html, page_url):
    """Discover downloadable tables from an actual source page."""
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(html, "html.parser")
    return tuple(sorted({urljoin(page_url, a["href"]) for a in soup.select("a[href]")
                         if asset_type(urlsplit(a["href"]).path)[0] in
                         ("table", "variants", "annotation", "reference")}))
