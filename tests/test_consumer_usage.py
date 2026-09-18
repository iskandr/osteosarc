"""Small offline examples matching the four consumers' acquisition sites.

Native consumer integrations run when their optional packages are installed;
the consumers CI job installs the reviewed versions explicitly.
"""

from dataclasses import replace

import pytest

from osteosarc import (
    Asset,
    Assets,
    Cache,
    Region,
    Variant,
    Variants,
    extract_reads,
    inspect_alignment,
)
from osteosarc.cache import stable_id
from osteosarc.catalog import asset_type


def cached_asset(dataset, tmp_path, name, text):
    url = "https://example.test/" + name
    original = tmp_path / name
    original.write_text(text)
    dataset.cache.import_file(original, url)
    kind, format = asset_type(name)
    asset = Asset(stable_id(url), name, url, kind, format)
    dataset.assets = Assets([*dataset.assets, asset])
    return asset


def test_varcode_custom_reference_preserves_identity_and_metadata(dataset, tmp_path):
    pyensembl = pytest.importorskip("pyensembl")
    pytest.importorskip("varcode")
    genome = pyensembl.Genome(reference_name="GRCh38-osteosarc-six-transcript-subset",
                             annotation_name="fixture", annotation_version=87,
                             gtf_path_or_url=str(tmp_path / "not-downloaded.gtf"))
    selected = dataset.variants(status="ready", gene="SMC5")
    original = selected[0]
    # Two website entries may normalize to one native allele; neither ID is lost.
    selected = Variants([original, replace(original, id="second-entry")], source=selected.source)[:]
    native = selected.to_varcode(genome=genome, assembly="GRCh38")
    assert len(native) == 1
    variant = native[0]
    assert variant.ensembl is genome
    assert variant.ensembl.reference_name == "GRCh38-osteosarc-six-transcript-subset"
    metadata = native.metadata[variant]
    assert [entry["id"] for entry in metadata["entries"]] == [original.id, "second-entry"]
    assert metadata["source"]["snapshot_id"] == dataset.id
    with pytest.raises(ValueError, match="custom-named"):
        selected.to_varcode(genome=genome)
    with pytest.raises(ValueError, match="match variant"):
        selected.to_varcode(genome=genome, assembly="GRCh37")
    with pytest.raises(ValueError, match="conflicts"):
        selected.to_varcode(genome=pyensembl.EnsemblRelease(95), assembly="GRCh37")
    with pytest.raises(ValueError, match="missing_literal"):
        dataset.variants(gene="COL3A1").to_varcode(genome=genome, assembly="GRCh38")
    assert not (tmp_path / "not-downloaded.gtf").exists()


def test_isovar_accepts_native_variants_and_cached_alignments(bam, tmp_path):
    isovar = pytest.importorskip("isovar")
    pyensembl = pytest.importorskip("pyensembl")
    selected = Variants([Variant("toy", "GENE", "GRCh38", (("chr1", 106, "A", "C"),), "ready")])
    native = selected.to_varcode(genome=pyensembl.EnsemblRelease(95))
    subset = extract_reads(bam, selected.regions(padding=10), cache=tmp_path / "cache")
    with subset.open() as alignment:
        evidence = isovar.ReadCollector().read_evidence_for_variant(native[0], alignment)
    assert evidence.alt_read_names == {"repeated"}
    assert not evidence.ref_reads
    assert not evidence.other_reads


def test_topiary_rsem_download_and_raw_table(dataset, tmp_path):
    asset = cached_asset(dataset, tmp_path, "sample.genes.results",
                         "gene_id\ttranscript_id(s)\tTPM\nENSG000001\tENST000001\t0\n")
    assert asset.kind == "expression"
    assert list(dataset.table(asset.key)) == [{"gene_id": "ENSG000001", "transcript_id(s)": "ENST000001", "TPM": "0"}]
    path = dataset.download(asset)
    assert path.name == "sample.genes.results"
    pytest.importorskip("topiary")
    from topiary.rna.expression_loader import load_expression
    expression = load_expression(path)
    assert expression.to_dict("records") == [{"gene_id": "ENSG000001", "TPM": 0}]


@pytest.mark.parametrize("body", ["", "AAAAAAAAA\tHLA-A*02:01\t100\t0.5\tGENE\n"])
def test_topiary_pvac_original_columns_and_header_only_reports(dataset, tmp_path, body):
    asset = cached_asset(dataset, tmp_path, "sample.all_epitopes.aggregated.tsv",
                         "Best Peptide\tAllele\tIC50 MT\t%ile MT\tGene\n" + body)
    raw = dataset.table(asset)
    assert raw.columns == ("Best Peptide", "Allele", "IC50 MT", "%ile MT", "Gene")
    assert len(raw) == bool(body)
    topiary = pytest.importorskip("topiary")
    result = topiary.read_pvacseq(dataset.download(asset))
    assert result.metadata.extra["pvacseq_format"] == "aggregated"
    assert result.df.empty == (not body)
    if body:
        assert set(result.df["peptide"]) == {"AAAAAAAAA"}


def test_vaxrank_corpus_header_survey_and_paired_mate_recovery(bam, tmp_path):
    import pysam
    # Vaxrank surveys headers before choosing native/lifted regions, then asks
    # for mates outside the panel for short-read libraries. No ranking wrapper.
    paired = tmp_path / "paired.bam"
    with pysam.AlignmentFile(bam) as template, pysam.AlignmentFile(paired, "wb", template=template) as output:
        for flag, start, mate, length in [(99, 100, 1000, 940), (147, 1000, 100, -940)]:
            read = pysam.AlignedSegment(output.header)
            read.query_name = "pair"
            read.query_sequence = "ACGT" * 10
            read.query_qualities = pysam.qualitystring_to_array("I" * 40)
            read.flag, read.reference_id, read.reference_start = flag, 0, start
            read.mapping_quality, read.cigarstring = 60, "40M"
            read.next_reference_id, read.next_reference_start, read.template_length = 0, mate, length
            output.write(read)
    cache = Cache(tmp_path / "cache", offline=True)
    info = inspect_alignment(paired, cache=cache)  # survey works before indexing
    assert info.assembly == "GRCh38"
    pysam.index(str(paired))
    regions = [Region("1", 100, 101, info.assembly)]
    regional = extract_reads(paired, regions, cache=cache)
    with_mates = extract_reads(paired, regions, cache=cache, fetch_pairs=True)
    assert regional.receipt["records"] == 1
    assert with_mates.receipt["records"] == 2
    with with_mates.open() as reads:
        assert [(r.query_name, r.reference_start) for r in reads] == [("pair", 100), ("pair", 1000)]
    assert with_mates.receipt["scope"] == "regional_records_and_paired_mates"


def test_vaccine_design_inputs_keep_sequences_and_untested_state(dataset):
    # Existing Vaxrank/Isovar objects remain downstream; use this package for
    # public comparator peptides, vaccine membership, and experimental states.
    targets = dataset.variants("vaccine", status="ready")
    assert set(targets.select(gene="SMC5")[0].vaccines) >= {"mRNA"}
    peptides = dataset.vaccine_peptides("mRNA")
    assert len(peptides) > 0
    assert all(p["variant_id"] and p["gene"] for p in peptides)
    smc5 = dataset.vaccines.select(gene="SMC5").rows[0]
    assert smc5["elispot_status"] == "not_tested"
    assert smc5["elispot_response"] is None
