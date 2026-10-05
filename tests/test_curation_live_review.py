"""Verified October source edits retain old layouts and reject later drift."""

import copy
import json
import warnings
from importlib.resources import files
from pathlib import Path

import pytest

from osteosarc import Table, parse_variants
from osteosarc.cache import digest
from osteosarc.curation import ASSAYS, CORRECTIONS, Curation, CurationWarning, normalize_tissue

DATA = Path(__file__).parent / "data"
REVIEW = json.loads((DATA / "curation-review-2026-10-05.json").read_text())
REVISIONS = json.loads(files("osteosarc").joinpath("data/curation-revisions.json").read_text())
RULES = {c.id: c for c in CORRECTIONS if c.id in REVISIONS["corrections"]}
ATTACHED = {"ush2a-transposed-duplicate", "allele-FAM157A-p_W70_Q71ins_14",
            "allele-COL3A1-Splice", "otud4-source-unavailable"}


def test_reviewed_layouts_have_no_stale_rules_or_source_mutations():
    raw = copy.deepcopy(REVIEW["sources"])
    before = copy.deepcopy(raw)
    curated = Curation(RULES.values(), raw.get)
    with warnings.catch_warnings():
        warnings.simplefilter("error", CurationWarning)
        report = {r["id"]: r["status"] for r in curated.report()}
        for source in raw:
            curated.records(source)
    assert report == {cid: "applied" if cid in ATTACHED else "fixed_upstream" for cid in RULES}
    assert raw == before
    assert curated.records("vafs")[0] == before["vafs"]
    source = {r["id"]: r for r in curated.records("source_variants")[0]}
    assert source["MAP2-chr2-209694768"]["variant_type"] == "delins"
    assert source["DCHS2-chr4-154323273"]["variant_type"] == "complex"
    assert source["GAPVD1-chr9-125301980"]["variant_type"] == "complex"
    assert source["USH2A-chr1-215650752"]["upstream_merge"]["retired_id"] == "USH2A-chr1-215560752"


@pytest.mark.parametrize("cid,source,key,value,field,changed", [
    ("gene-symbol-TRMO", "source_variants", "id", "TRMO-chr9-97910412", "ref", "G"),
    ("allele-DCHS2-chr4-154322488", "vafs", "variant_id", "DCHS2-chr4-154323273", "alt_reads", "1"),
    ("allele-GAPVD1-chr9-125299105", "source_variants", "id", "GAPVD1-chr9-125301980", "alt", "C"),
    ("allele-MAP2-chr2-209694768", "source_variants", "id", "MAP2-chr2-209694768", "pos", 209694769),
    ("fam157a-withdrawn-protein", "source_variants", "id", "FAM157A-p_W70_Q71ins_14", "note", "Protein verified"),
    ("allele-FAM157A-p_W70_Q71ins_14", "source_variants", "id", "FAM157A-p_W70_Q71ins_14", "refseq_id", "NM_000001"),
    ("allele-COL3A1-Splice", "vafs", "variant_id", "COL3A1-Splice", "pos", "189010890"),
    ("otud4-source-unavailable", "source_variants", "id", "OTUD4-p_A153del", "ref", "G"),
    ("ush2a-transposed-duplicate", "source_variants", "id", "USH2A-chr1-215650752", "genomic_ref_context", "AAAA"),
    ("specimen-T1-site", "specimens", "sample_id", "T1_tumor", "tissue_source", "Tumor resection"),
    ("specimen-T3-site", "specimens", "sample_id", "T3_tumor", "collection_site", "UCSF"),
    ("tempus-file-labels", "bam_metadata", "s3_path", "vendor/tempus/TL-24-ALMY2X4KMV/DNA/TL-24-ALMY2X4KMV_N.sorted.bam", "timepoint", "T0"),
    ("reyagel-end-date", "events", "title", "ReyaGel", "end_date", "2027-07-14"),
])
def test_unreviewed_future_changes_are_stale_and_atomic(cid, source, key, value, field, changed):
    raw = copy.deepcopy(REVIEW["sources"])
    next(r for r in raw[source] if r.get(key) == value)[field] = changed
    before = copy.deepcopy(raw)
    curated = Curation([RULES[cid]], raw.get)
    with pytest.warns(CurationWarning):
        assert curated.report()[0]["status"] == "stale"
    for name in raw:
        assert curated.records(name)[0] == before[name]


@pytest.mark.parametrize("cid,source,key,restored", [
    ("map2-split-representations", "variant_index", "id", "MAP2-chr2-209694772"),
    ("gene-symbol-TRMO", "source_variants", "id", "TMRO-chr9-97910412"),
    ("tempus-grch37-counts", "vafs", "bam_file", "TL-24-5GQLV9WSXQ_T.sorted.bam"),
])
def test_reintroduced_retired_records_do_not_match_a_reviewed_removal(cid, source, key, restored):
    raw = copy.deepcopy(REVIEW["sources"])
    raw[source].append({key: restored})
    curated = Curation([RULES[cid]], raw.get)
    with pytest.warns(CurationWarning):
        assert curated.report()[0]["status"] == "stale"


def test_public_otud4_allele_has_primary_call_and_reference_evidence_without_somatic_claim():
    folder = DATA / "otud4-review"
    provenance = json.loads((folder / "provenance.json").read_text())
    for filename, entry in provenance.items():
        assert digest(folder / filename) == entry["fixture_sha256"]
        if entry["selection"].startswith("whole"):
            assert entry["sha256"] == entry["fixture_sha256"]
    raw = copy.deepcopy(REVIEW["sources"])
    curated = Curation(RULES.values(), raw.get)
    variants = parse_variants(curated.records("variant_index")[0], Table(curated.records("vafs")[0]),
                              source_variants=curated.records("source_variants")[0])
    variant = variants["OTUD4-p_A153del"]
    assert variant.status == "ready" and variant.allele == ("chr4", 145155970, "TCAG", "T")
    evidence = variant.annotations["allele_resolution"]
    assert evidence["source"]["sha256"] == provenance["cegat-normal-indels.vcf"]["sha256"]
    assert evidence["annotation_source"]["sha256"] == provenance["cegat-normal-indels.annotated.tsv"]["sha256"]
    row, = (folder / "cegat-normal-indels.vcf").read_text().splitlines()
    fields = row.split("\t")
    assert fields[:5] == ["chr4", "146077122", ".", "TCAG", "T"]
    assert "DP=323" in fields[7] and "AD=7" in fields[7]
    annotation, = [line.split("\t") for line in (folder / "cegat-normal-indels.annotated.tsv").read_text().splitlines()
                   if "NM_001102653.1" in line]
    assert annotation[4] == evidence["source_filter"] == "FP/HET"
    assert annotation[12:14] == [evidence["source_hgvs_c"], evidence["source_hgvs_p"]]
    sequences = []
    for reference, allele in zip(evidence["references"], [evidence["source_allele"], evidence["allele"]]):
        sequence = "".join((folder / (reference["accession"] + ".fasta")).read_text().splitlines()[1:])
        offset = allele["pos"] - reference["start"]
        assert sequence[offset:offset + len(allele["ref"])] == allele["ref"]
        sequences.append(sequence)
    assert sequences[0] == sequences[1]
    mapping, = json.loads((folder / "mapping.json").read_text())["mappings"]
    assert mapping == evidence["mapping"]["mappings"][0]
    assert mapping["mapped"]["start"] == evidence["allele"]["pos"]
    assert evidence["original_natera_report"] == "unavailable"
    assert evidence["somatic_status"] == "not_established"
    label = json.loads((folder / "source-label.json").read_text())
    assert label["row"]["tissue"] == "Blood" and "normal DNA" in label["row"]["notes"]
    assert evidence["source_tissue"] == "blood"
    assert evidence["source_metadata"]["sha256"] == label["receipt"]["sha256"]


def test_reviewed_labels_preserve_assay_and_composite_tissue_distinctions():
    assert ASSAYS["Panel"] == ("panel", None)
    assert normalize_tissue("blood normal (drawn 2024-06-11)") == "blood"
    assert normalize_tissue("organoid (read 1)") == "organoid"
    assert normalize_tissue("tumor + blood (drawn 2024-09-21)") is None
    assert normalize_tissue("normal (blood, drawn 2027-01-01)") == "normal (blood, drawn 2027-01-01)"
