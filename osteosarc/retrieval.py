"""Native coordinates for historical retrieval fixtures; no implicit liftover."""

import copy
import json
from dataclasses import asdict
from pathlib import Path

import pysam

from .bundles import safe_path, verify_manifest_files
from .cache import stable_id
from .errors import CoordinateError, IntegrityError
from .models import Region
from .reads import assembly_from_header, normalize_assembly, resolve_regions


def annotate_retrieval_case(case, bam_path):
    """Copy a historical case with a verified native query and reference identity.

    ``variant_id`` is a catalogue label, never a coordinate to parse. The original
    variant and its reference-validation receipt survive unchanged. The BAM must
    have an index and contain at least one alignment overlapping the native locus.
    """
    result = copy.deepcopy(case)
    variant = case["variant"]
    with pysam.AlignmentFile(bam_path) as bam:
        header = bam.header.to_dict()
        assembly = assembly_from_header(header)
        if assembly is None or assembly != normalize_assembly(variant["assembly"]):
            raise CoordinateError(f"Case/header assembly differs: {case['case_id']}")
        validation = variant["genomic_validation"]
        genome = validation["genome"]
        mitochondrial = variant["chrom"].removeprefix("chr") in ("M", "MT")
        # hg19's 16,571-base chrM is not the 16,569-base rCRS used by hg38
        # and some GRCh37 alignments. The frozen reference validation determines
        # this choice; nuclear assembly and name aliases cannot determine it.
        length = {"hg19": 16571, "hg38": 16569}.get(genome) if mitochondrial else None
        if mitochondrial and length is None:
            raise CoordinateError("Unknown mitochondrial validation reference")
        if not mitochondrial and normalize_assembly(genome) != assembly:
            raise CoordinateError("Nuclear validation reference differs from the BAM assembly")
        pos, ref = variant["pos"], variant["ref"]
        offset = pos - 1 - validation["start"]
        if offset < 0 or validation["sequence"][offset:offset + len(ref)].upper() != ref.upper():
            raise IntegrityError(f"Reference validation does not match the native allele: {case['case_id']}")
        region, = resolve_regions([Region(variant["chrom"], pos - 1, pos - 1 + len(ref),
                                          assembly, length)], header)
        if not any(bam.fetch(region.contig, region.start, region.end)):
            raise IntegrityError(f"No indexed reads at the native locus: {case['case_id']}")
        sq = next(row for row in header["SQ"] if row["SN"] == region.contig)
        dictionary = [{k: row[k] for k in ("SN", "LN", "M5") if k in row} for row in header["SQ"]]
    result["catalogue_identity"] = copy.deepcopy(variant.get("original_identity", {
        k: variant[k] for k in ("assembly", "chrom", "pos", "ref", "alt")}))
    result["native_variant"] = dict(assembly=assembly, chrom=region.contig, pos=pos,
                                    ref=ref, alt=variant["alt"], coordinates="one-based")
    result["query_region"] = asdict(region)
    result["reference"] = dict(nuclear_assembly=assembly, validation_genome=genome, sequence=sq,
        mitochondrial_reference=("rCRS" if length == 16569 else "UCSC-hg19-chrM") if mitochondrial else None,
        dictionary_sha256=stable_id(dictionary))
    return result


def annotate_retrieval_manifest(manifest, directory):
    """Verify original assets, then return a new, separately versioned manifest."""
    directory = Path(directory)
    verify_manifest_files(directory, {a["filename"]: a for a in manifest["assets"]})
    result = copy.deepcopy(manifest)
    result["parent"] = dict(data_version=manifest["data_version"], manifest_content_sha256=stable_id(manifest))
    result["data_version"] = "native-retrieval-v1"
    result["cases"] = [annotate_retrieval_case(case, safe_path(directory, case["bam"]))
                       for case in manifest["cases"]]
    result["merge_policy"] = "Check reference dictionaries before merging; query native coordinates, never variant_id."
    return result


def assert_compatible_retrieval_cases(cases):
    """Refuse a raw BAM merge across differing builds, contig names or references.

    Matching dictionaries only establish coordinate compatibility. They do not
    establish biological independence or make balanced fixtures quantitative.
    """
    identities = {(c["reference"]["nuclear_assembly"], c["reference"]["dictionary_sha256"]) for c in cases}
    if len(identities) > 1:
        raise CoordinateError("Retrieval cases use different reference dictionaries; query each BAM with its "
                              "query_region, or explicitly harmonize references and coordinates before merging")


def load_retrieval_cases():
    """Reviewed native metadata for the 49 historical shared-v1 retrieval cases."""
    return json.loads((Path(__file__).with_name("data") / "retrieval-cases.json").read_text())
