import copy
import json
from pathlib import Path

import pytest

from osteosarc import CoordinateError, IntegrityError
from osteosarc.retrieval import (
    annotate_retrieval_case,
    annotate_retrieval_manifest,
    assert_compatible_retrieval_cases,
    load_retrieval_cases,
)

DATA = Path(__file__).parent / "data/retrieval"


def test_native_coordinates_and_contigs_reproduce_reviewed_cases():
    original = json.loads((DATA / "historical-manifest.json").read_text())
    before = copy.deepcopy(original)
    annotated = annotate_retrieval_manifest(original, DATA)
    shipped = {c["case_id"]: c for c in load_retrieval_cases()["cases"]}
    assert len(shipped) == 49
    assert annotated["data_version"] != original["data_version"]
    assert original == before
    for case in annotated["cases"]:
        assert case == shipped[case["case_id"]]
    normal, rcrs, hg19, nr2f2 = annotated["cases"]
    assert normal["native_variant"]["chrom"] == "chr7"
    assert rcrs["native_variant"]["chrom"] == "MT"
    assert rcrs["native_variant"]["assembly"] == "GRCh37"
    assert rcrs["native_variant"]["pos"] == 12994
    assert rcrs["reference"]["mitochondrial_reference"] == "rCRS"
    assert rcrs["reference"]["validation_genome"] == "hg38"
    assert hg19["query_region"]["reference_length"] == 16571
    assert hg19["native_variant"]["pos"] == 12995
    assert nr2f2["native_variant"]["pos"] == 96875528
    assert nr2f2["catalogue_identity"]["pos"] == 96332299
    assert_compatible_retrieval_cases([hg19, nr2f2])
    with pytest.raises(CoordinateError, match="different reference dictionaries"):
        assert_compatible_retrieval_cases(annotated["cases"])
    with pytest.raises(CoordinateError, match="different reference dictionaries"):
        assert_compatible_retrieval_cases([rcrs, hg19])


def test_bad_native_allele_or_reference_is_an_error():
    cases = json.loads((DATA / "historical-manifest.json").read_text())["cases"]
    changed = copy.deepcopy(cases[-1])
    changed["variant"]["pos"] = 96332299  # ID is GRCh38; this BAM is GRCh37.
    with pytest.raises(IntegrityError, match="native allele"):
        annotate_retrieval_case(changed, DATA / changed["bam"])
    changed = copy.deepcopy(cases[1])
    changed["variant"]["genomic_validation"]["genome"] = "hg19"
    with pytest.raises(CoordinateError, match="Reference length mismatch"):
        annotate_retrieval_case(changed, DATA / changed["bam"])
