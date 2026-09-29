"""Completeness and anti-false-negative checks on the frozen real evidence."""

import copy
from pathlib import Path

import pytest

from osteosarc import IntegrityError
from osteosarc.cache import digest
from osteosarc.corpus import (
    COUNT_FIELDS,
    count_matrix,
    counts,
    evidence_catalog,
    negative_controls,
    support_status,
    vaccine_membership,
)
from osteosarc.shared import read_json

FOLDER = Path(__file__).parent / "data/comprehensive"


@pytest.fixture(scope="module")
def inputs():
    return read_json(FOLDER / "inputs.json.gz")


@pytest.fixture(scope="module")
def catalog(inputs):
    return evidence_catalog(inputs)


def test_real_catalogue_rebuilds_exactly_and_inputs_are_pinned(catalog):
    for name, checksum in read_json(FOLDER / "checksums.json").items():
        assert digest(FOLDER / name) == checksum
    assert catalog == read_json(FOLDER / "catalog.json.gz")
    assert catalog["summary"] == read_json(FOLDER / "summary.json")


def test_every_vaccine_assertion_is_preserved(inputs, catalog):
    vaccinated = {vid for vid, v in catalog["variants"].items() if v["vaccine"]["included"]}
    assert len(vaccinated) == 51
    index = {v["id"] for v in inputs["variants"] if (v.get("vaccine_count") or 0) > 0}
    assert len(index) == 44 and index <= vaccinated
    assert vaccinated - index == {
        "ABI3BP-chr3-100850671", "ANKRD17-chr4-73097139", "ASPM-chr1-197102716",
        "EXD3-chr9-137373504", "PPP1R3F-chrX-49270761", "PRRC2C-chr1-171540477", "SPG11-chr15-44660435"}
    assert len(catalog["vaccine_rows"]) == len(inputs["vaccine_rows"]) == 38
    tecpr1 = catalog["variants"]["TECPR1-chr7-98241127"]["vaccine"]
    assert len(tecpr1["overlap_rows"]) == 2
    assert set(tecpr1["overlap_union"]) == {"mRNA", "JLF V1", "JLF V2", "JLF V3", "CeGaT"}
    assert catalog["variants"]["ABCF2-chr7-151218156"]["vaccine"]["vaccines_union"] == ["Cure 2024"]


def test_matrix_has_every_variant_product_pair_including_missing(catalog):
    matrix = catalog["allele_support"]
    expected = {(vid, bam) for vid in catalog["variants"] for bam in catalog["count_sources"]}
    assert len(matrix) == len(expected) == 7474
    assert {(r["variant_id"], r["bam_file"]) for r in matrix} == expected
    assert all(r["total_reads"] is None for r in matrix if r["status"] == "not_measured")
    assert all(r["total_reads"] == 0 for r in matrix if r["status"] == "no_coverage")
    assert all(r["ref_reads"] > 0 and r["alt_reads"] == 0 for r in matrix if r["status"] == "covered_no_alt")
    assert {r["timepoint"] for r in matrix if r["assay_type"] in {"RNA", "scRNA", "scRNA_ONT"}} == {
        "T0", "T1", "T2", "T3"}


def test_negative_controls_have_matched_dna_and_explicit_rna_coverage(catalog):
    assert len(catalog["negative_controls"]) == 4
    for control in catalog["negative_controls"]:
        assert control["tumor_rna_products_checked"] == 15
        assert control["comparisons"]
    classes = {c["variant_id"]: c["control_class"] for c in catalog["negative_controls"]}
    assert classes == {
        "ACE-chr17-63497361": "covered_rna_no_alt",
        "BTN3A3-chr6-26448446": "covered_rna_no_alt",
        "CTSE-chr1-206022947": "insufficient_rna_coverage",
        "DBH-chr9-133652279": "insufficient_rna_coverage",
    }
    for c in catalog["negative_controls"]:
        assert all(bool(m["covered_rna"]) == (c["control_class"] == "covered_rna_no_alt")
                   for m in c["comparisons"])
    btn = next(c for c in catalog["negative_controls"] if c["variant_id"].startswith("BTN3A3-"))
    assert [(r["tissue"], r["assay"], r["depth"], r["vaf"]) for r in btn["nonzero_rna_summary_exceptions"]] == [
        ("blood", "CITE", 435, 0.0023)]


def test_historical_rna_prevents_false_negative_claims(inputs, catalog):
    disagreements = catalog["known_rna_disagreements"]
    assert {r["variant_id"] for r in disagreements} == {"NR2F2-chr15-96332299"}
    assert any(r["observed_before_sampling"]["alt"] == 7 for r in disagreements)
    changed = copy.deepcopy(inputs)
    row = copy.deepcopy(disagreements[0])
    row["variant_id"] = "ACE-chr17-63497361"
    changed["historical_rna"].append(row)
    with pytest.raises(IntegrityError, match="historical RNA support"):
        evidence_catalog(changed)


def test_other_tumor_rna_summaries_cannot_be_ignored(inputs):
    changed = copy.deepcopy(inputs)
    variant = next(v for v in changed["variants"] if v["id"] == "ACE-chr17-63497361")
    row = next(r for r in variant["annotations"]["source_record"]["vafs"]
               if r["assay"] == "RNA" and r["tissue"] == "tumor")
    row["vaf"] = 0.01
    with pytest.raises(IntegrityError, match="nonzero tumor RNA summary"):
        evidence_catalog(changed)


@pytest.mark.parametrize("change", [{"alt_reads": 1}, {"alt_reads": None}, {"status": "not_measured"}])
def test_a_negative_control_cannot_hide_positive_or_missing_rna(inputs, catalog, change):
    matrix = copy.deepcopy(catalog["allele_support"])
    row = next(r for r in matrix if r["variant_id"] == "ACE-chr17-63497361" and r["assay_type"] == "RNA")
    row.update(change)
    with pytest.raises(IntegrityError, match="alternate or unmeasured"):
        negative_controls(inputs["variants"], matrix)


def test_negative_control_requires_covered_rna_at_dna_timepoint(inputs, catalog):
    matrix = copy.deepcopy(catalog["allele_support"])
    for row in matrix:
        if row["variant_id"] == "ACE-chr17-63497361" and row["timepoint"] == "T1" and row["assay_type"] in {"RNA", "scRNA", "scRNA_ONT"}:
            row.update(ref_reads=0, alt_reads=0, other_reads=0, total_reads=0, status="no_coverage")
    with pytest.raises(IntegrityError, match="No matched"):
        negative_controls(inputs["variants"], matrix)


def test_missing_zero_other_only_and_invalid_counts_are_distinct():
    assert support_status(counts({})) == "not_measured"
    assert support_status(counts(dict.fromkeys(COUNT_FIELDS, "0"))) == "no_coverage"
    assert support_status(counts(dict(ref_reads=0, alt_reads=0, other_reads=3, total_reads=3))) == "other_only"
    for bad in ("-1", "1.5", "nan", True):
        with pytest.raises(IntegrityError):
            counts(dict(ref_reads=bad))
    with pytest.raises(IntegrityError, match="differs from depth"):
        counts(dict(ref_reads=10, alt_reads=2, other_reads=0, total_reads=10))


def test_ambiguous_vaccine_join_duplicate_evidence_and_wrong_allele_fail(inputs):
    duplicate = copy.deepcopy(inputs["variants"])
    duplicate.append(next(v for v in duplicate if v["id"] == "SMC5-chr9-70298024"))
    with pytest.raises(IntegrityError, match="joins 2 alleles"):
        vaccine_membership(duplicate, inputs["vaccine_rows"])
    with pytest.raises(IntegrityError, match="Duplicate evidence"):
        count_matrix(inputs["variants"], inputs["counts"] + inputs["counts"][:1])
    rows = copy.deepcopy(inputs["counts"])
    rows[0]["pos"] = str(int(rows[0]["pos"]) + 1)
    with pytest.raises(IntegrityError, match="another allele"):
        count_matrix(inputs["variants"], rows)


def test_structural_catalogue_keeps_unsupported_and_exact_allele_evidence(catalog):
    candidates = catalog["structural_candidates"]["targets"]
    assert len(candidates) == 637
    evidence = [e for v in candidates.values() for e in v.get("rna_evidence", [])]
    assert any(e["status"] == "no_observed_adjacency" for e in evidence)
    assert any(e.get("exact_DNA_allele_20bp_templates") for e in evidence)
    assert any(e.get("exact_DNA_allele_20bp_templates") is None for e in evidence)


def test_new_spec_preserves_old_release_and_sets_explicit_fusion_sides():
    from osteosarc.shared import structural_targets
    spec = read_json(FOLDER / "spec.json")
    assert spec["selection"]["split_depth"] == 1_000_000
    targets, _, _ = structural_targets(spec)
    for name, sides in {"GABBR1--SLC29A1": ["right", "right"], "OTUD7A--FMN1": ["right", "right"],
                        "PARD3B--CDKN2B-AS1-CDKN2B": ["left", "right"]}.items():
        assert [end["retained_side"] for end in targets[name]["breakends"]] == sides


def test_independent_snv_recount_collapses_mates_and_preserves_quality_unknowns():
    import pysam

    from osteosarc.corpus import snv_read_support
    header = pysam.AlignmentHeader.from_dict(dict(SQ=[dict(SN="chr1", LN=100)]))

    def read(name, base="T", quality=40, mapq=60, flag=0, nh=None, rg=None):
        r = pysam.AlignedSegment(header)
        r.query_name, r.reference_id, r.reference_start = name, 0, 0
        r.flag, r.mapping_quality, r.cigarstring = flag, mapq, "10M"
        r.query_sequence = "AAA" + base + "AAAAAA"
        r.query_qualities = None if quality is None else [quality] * 10
        if nh is not None:
            r.set_tag("NH", nh)
        if rg is not None:
            r.set_tag("RG", rg)
        return r

    records = [read("alt"), read("alt"), read("ref", "G"),
               read("conflict", "G"), read("conflict", "T"),
               read("same-name", rg="a"), read("same-name", rg="b"),
               read("star", mapq=255, nh=1), read("unknown-mapq", mapq=255),
               read("no-qual", quality=None), read("unknown-qual", quality=255),
               read("low-bq", quality=19), read("low-mq", mapq=19),
               read("duplicate", flag=1024), read("secondary", flag=256), read("supplementary", flag=2048)]
    result = snv_read_support(records, 4, "G", "T")
    assert result["templates"] == dict(alt=4, ref=1, other=1)
    assert result["filtered_records"] == dict(flag_filtered=3, low_base_quality=1, low_mapping_quality=1,
                                               missing_base_quality=2, unknown_mapping_quality=1)
    assert result["measured_before_sampling"] is True
    with pytest.raises(ValueError, match="single-base"):
        snv_read_support(records, 4, "GG", "T")
