"""Catalogue resolutions checked against public VCFs and versioned RefSeq bases."""

import hashlib
import json
from copy import deepcopy
from importlib.resources import files
from pathlib import Path

import pytest

from osteosarc import CORRECTIONS, Dataset, Table, parse_variants
from osteosarc.curation import Curation, CurationWarning
from osteosarc.parsing import parse_table

DATA = Path(__file__).parent / "data" / "allele_resolution"
REVIEW = json.loads(files("osteosarc").joinpath("data/allele_resolutions.json").read_text())
FAM = "FAM157A-p_W70_Q71ins_14"
COL = "COL3A1-Splice"
MUC = "MUC3A-chr7-100953130"
USH = "USH2A-chr1-215560752"
TARGET = "USH2A-chr1-215650752"
OTUD = "OTUD4-p_A153del"


def sources():
    return dict(source_variants=json.loads((DATA / "source-variants.json").read_text()),
                variant_index=json.loads((DATA / "variant-index.json").read_text()),
                vafs=list(parse_table((DATA / "vafs.tsv").read_text())))


def curated(raw=None, enabled=True):
    raw = sources() if raw is None else raw
    # Other registry entries need sources outside this six-entry excerpt.
    corrections = [c for c in CORRECTIONS if any("allele_resolution" in ch.set for ch in c.changes)]
    curation = Curation(corrections, raw.get, enabled=enabled)
    variants = parse_variants(curation.records("variant_index")[0],
                              Table(curation.records("vafs")[0]),
                              source_variants=curation.records("source_variants")[0])
    return variants, curation


def sequence(reference):
    manifest = json.loads((DATA / "provenance.json").read_text())["files"]
    source = manifest[reference["accession"] + ".fasta"]
    assert {k: reference[k] for k in ("url", "sha256")} == {k: source[k] for k in ("url", "sha256")}
    return "".join((DATA / (reference["accession"] + ".fasta")).read_text().splitlines()[1:])


def vcf_record(evidence):
    name = "tempus.freebayes.vcf" if "freebayes" in evidence["source"]["url"] else "tempus.pindel.vcf"
    source = json.loads((DATA / "provenance.json").read_text())["files"][name]
    assert evidence["source"] == {k: source[k] for k in ("url", "sha256")}
    text = (DATA / name).read_text()
    assert "##reference=GRCh37" in text or "##reference=grch37" in text
    allele = evidence["source_allele"]
    matches = [line.split("\t") for line in text.splitlines() if not line.startswith("#")
               and line.split("\t")[:2] == [allele["chrom"].removeprefix("chr"), str(allele["pos"])]]
    row, = matches
    assert row[3:5] == [allele["ref"], allele["alt"]]
    return row


def test_fixture_checksums():
    for entry in json.loads((DATA / "provenance.json").read_text())["files"].values():
        checksum = hashlib.sha256((DATA / entry["fixture"]).read_bytes()).hexdigest()
        assert checksum == entry["fixture_sha256"]
        if entry["selection"] == "whole file":
            assert checksum == entry["sha256"]


@pytest.mark.parametrize("vid,delta", [(FAM, 42), (COL, -737)])
def test_resolved_alleles_match_vcf_refseq_and_equivalent_haplotypes(vid, delta):
    evidence = REVIEW[vid]["resolution"]
    row = vcf_record(evidence)
    assert evidence["transcript"] in row[7]
    assert evidence["source_hgvs_c"] in row[7]
    a, b = evidence["references"]
    old, new = evidence["source_allele"], evidence["allele"]
    assert sequence(a) == sequence(b)  # full flanks and deletion, not just the anchor
    for allele, reference in ((old, a), (new, b)):
        seq = sequence(reference)
        offset = allele["pos"] - reference["start"]
        assert seq[offset:offset + len(allele["ref"])] == allele["ref"]
    block, = evidence["mapping"]["blocks"]
    assert block["strand"] == "+"
    assert new["pos"] - old["pos"] == block["target"][1] - block["source"][1]
    assert (old["ref"], old["alt"]) == (new["ref"], new["alt"])
    assert len(new["alt"]) - len(new["ref"]) == delta
    seq = sequence(b)
    offset = new["pos"] - b["start"]
    ref, alt = new["ref"], new["alt"]
    haplotype = seq[:offset] + alt + seq[offset + len(ref):]
    equivalent = []
    # Exhaustively compare alternate haplotypes across the checked window.
    # This verifies left alignment, including rotation within FAM157A's repeat.
    for i in range(len(seq) - len(ref) + 1):
        alternative = haplotype[i:i + len(alt)]
        if (alternative[0] == seq[i]
                and seq[:i] + alternative + seq[i + len(ref):] == haplotype):
            equivalent.append(b["start"] + i)
    assert equivalent == evidence["normalization"]["equivalent_positions_grch38"]
    assert equivalent[0] == new["pos"]


def test_five_dispositions_preserve_ids_counts_and_protein_caveat():
    variants, curation = curated()
    assert len(variants) == 6
    assert {r["status"] for r in curation.report()} == {"applied"}
    for vid in (FAM, COL):
        variant = variants[vid]
        assert variant.status == "ready"
        allele = variant.annotations["allele_resolution"]["allele"]
        assert variant.allele == tuple(allele[k] for k in ("chrom", "pos", "ref", "alt"))
        assert variant.region().assembly == "GRCh38"
    assert curation.records("vafs")[0] == sources()["vafs"]  # no invented count rows
    protein = variants[FAM].annotations["allele_resolution"]["protein_interpretation"]
    assert protein["status"] == "withdrawn_protein_model"
    assert protein["current_transcript"] == "NR_146164.1"
    ncbi = (DATA / "FAM157A-ncbi.txt").read_text()
    assert "NM_001145248.1" in ncbi and "not for the protein" in ncbi
    for vid, outcome in ((MUC, "unresolved_grch38_placement"),
                         (OTUD, "source_unavailable"), (USH, "duplicate")):
        assert variants[vid].status != "ready"
        assert variants[vid].annotations["allele_resolution"]["status"] == outcome
        with pytest.raises(ValueError):
            variants[vid].region()
    assert variants[TARGET].allele == ("chr1", 215650752, "C", "A")
    assert variants[USH].alleles == (("chr1", 215560752, "A", "not_reported"),)
    raw, _ = curated(enabled=False)
    assert raw[FAM].status == raw[COL].status == "missing_literal_allele"
    assert all("allele_resolution" not in v.annotations for v in raw)


def test_muc3a_source_is_in_chain_gap_and_does_not_match_catalogue_reference():
    resolution = REVIEW[MUC]["resolution"]
    vcf_record(resolution)
    old = resolution["source_allele"]
    a, b = resolution["references"]
    assert sequence(a)[old["pos"] - a["start"]] == old["ref"]
    gap = resolution["mapping"]["source_gap"]
    assert gap[1] <= old["pos"] - 1 < gap[2]
    gap = resolution["mapping"]["target_gap"]
    catalog_pos = REVIEW[MUC]["expected_source"]["pos"]
    assert gap[1] <= catalog_pos - 1 < gap[2]
    assert sequence(b)[catalog_pos - b["start"]] != old["ref"]
    assert old["alt"][1:] not in sequence(b)
    assert "allele" not in resolution


def test_ush2a_evidence_supports_only_the_candidate_target():
    resolution = REVIEW[USH]["resolution"]
    relationship = resolution["relationship"]
    assert relationship["kind"] == "duplicate"
    assert relationship["variant_id"] == TARGET
    evidence = relationship["evidence"]
    row = vcf_record(evidence)
    assert evidence["source_hgvs_c"] in row[7]
    a, b = evidence["references"]
    assert sequence(a) == sequence(b)
    assert sequence(b)[215650752 - b["start"]] == "C"
    assert "allele" not in resolution  # A duplicate gets no second copy of the allele


def test_natera_source_has_no_otud4_allele_and_supports_col3a1_hgvs():
    rows = parse_table((DATA / "natera-catalogue.tsv").read_text())
    otud = rows.select(gene="OTUD4").rows[0]
    assert otud["Natera_2022"] == "yes" and otud["protein_change"] == "p.A153del"
    assert all(otud[k] == "" for k in ("refseq", "chr", "pos", "ref", "alt", "genomic_change_on_cdna"))
    col = rows.select(gene="COL3A1").rows[0]
    assert col["genomic_change_on_cdna"] == REVIEW[COL]["resolution"]["source_hgvs_c"]
    assert col["refseq"] == REVIEW[COL]["resolution"]["transcript"].split(".")[0]


@pytest.mark.parametrize("vid", list(REVIEW))
@pytest.mark.parametrize("source,field", [("source_variants", "genomic_change_on_cdna"),
                                         ("variant_index", "location")])
def test_upstream_drift_skips_the_whole_resolution(vid, source, field):
    raw = sources()
    row = next(r for r in raw[source] if r["id"] == vid)
    row[field] = "changed upstream"
    with pytest.warns(CurationWarning):
        variants, curation = curated(raw)
    assert "allele_resolution" not in variants[vid].annotations
    assert any(r["status"] == "stale" for r in curation.report())
    actual = next(r for r in curation.records("source_variants")[0] if r["id"] == vid)
    expected = next(r for r in raw["source_variants"] if r["id"] == vid)
    assert actual == expected


def test_ush2a_candidate_change_invalidates_relationship():
    raw = sources()
    next(r for r in raw["source_variants"] if r["id"] == TARGET)["alt"] = "G"
    with pytest.warns(CurationWarning):
        variants, _ = curated(raw)
    assert "allele_resolution" not in variants[USH].annotations


@pytest.mark.parametrize("pos,ref,alt,status", [
    ("NA", "G", "A", "malformed_source_row"),
    ("198153259", "G", "A", "conflicting_coordinates"),
])
def test_resolution_does_not_hide_malformed_or_conflicting_vafs(pos, ref, alt, status):
    raw = sources()
    raw["vafs"].append(dict(variant_id=FAM, gene="FAM157A", chrom="chr3", pos=pos, ref=ref, alt=alt))
    variants, _ = curated(raw)
    assert variants[FAM].status == status
    assert variants[COL].status == "ready"
    with pytest.raises(ValueError):
        variants[FAM].region()


def test_source_alleles_require_explicit_grch38_resolution():
    _, curation = curated()
    records = deepcopy(curation.records("source_variants")[0])
    index = curation.records("variant_index")[0]
    next(r for r in records if r["id"] == FAM).pop("allele_resolution")
    next(r for r in records if r["id"] == COL)["allele_resolution"]["allele"]["assembly"] = "GRCh37"
    variants = parse_variants(index, Table(sources()["vafs"]), source_variants=records)
    assert variants[FAM].status == variants[COL].status == "missing_literal_allele"


def test_dataset_resolution_provenance_and_missing_counts(dataset):
    for vid in (FAM, COL):
        v = dataset.variants()[vid]
        assert f"allele-{vid}" in v.annotations["corrections"]
        assert not dataset.vafs.select(variant_id=vid)
        raw = Dataset.open("fixture", cache=dataset.cache, corrections=False).variants()[vid]
        assert raw.annotations["source_record"]["pos"] is None
        assert "allele_resolution" not in raw.annotations
