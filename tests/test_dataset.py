import json
from collections import Counter

import pytest

from osteosarc import Dataset, IntegrityError
from osteosarc.cache import stable_id
from osteosarc.cli import main


def test_real_source_identity_and_unresolved_variants(dataset):
    variants = dataset.variants()
    assert len(variants) == 5
    assert variants["DYNC1H1-chr14-101980529"].allele == ("chr14", 101980529, "G", "A")
    assert variants["COL3A1-Splice"].status == "missing_literal_allele"
    assert len(dataset.variants(status="ready")) == 3
    deletion = variants["GTF3C5-chr9-133057893"].region()
    # The original anchored REF has 13 bases: anchor plus twelve deleted bases.
    assert (deletion.start, deletion.end) == (133057892, 133057905)
    with pytest.raises(ValueError):
        variants.regions()


def test_vaccine_peptides_annotations_and_missing_states(dataset):
    smc5 = dataset.variants(gene="SMC5")[0]
    assert {"mRNA", "JLF V3"} <= set(smc5.vaccines)
    assert dataset.variants(vaccine="mRNA", gene="SMC5")[0] == smc5
    assert not dataset.variants(vaccine="mRNA", vaccine_source="source_variants", gene="SMC5")
    assert "vaccine_membership_conflict" in smc5.annotations
    row = next(r for r in dataset.vaccines if r["gene"] == "SMC5")
    assert row["elispot_status"] == "not_tested"
    assert row["elispot_response"] is None
    assert len(dataset.vaccine_peptides("mRNA")) > 0
    assert smc5.annotations["source_record"]["consequence"] == "missense_variant"


def test_sample_conflicts_do_not_silently_become_replicates(dataset):
    raw = Dataset.open("fixture", cache=dataset.cache, corrections=False)
    assets = raw.assets.select(kind="alignment")
    assert len(assets) == 3
    source = assets.select(contains="BG009368")[0]
    assert set(source.values("timepoint")) == {"T0", "T1"}
    assert source.resolved("timepoint") is None
    assert not assets.select(contains="BG009368", timepoint="T0")
    assert len(assets.select(contains="BG009368", timepoint="T1", include_conflicts=True)) == 1
    assert source.metadata["catalog_genome_assertion"] == "hg38"
    assert "assembly" not in source.metadata
    assert len(dataset.samples) > len(assets)
    # The viewer label's stale date is surfaced alongside its timepoint.
    assert source.conflicts["date"] == ("2022-12", "2024-06")
    # The central correction for that stale label resolves it, traceably.
    corrected = dataset.asset(source.key)
    assert corrected.resolved("timepoint") == "T1" and not corrected.conflicts
    assert corrected.metadata["corrections"] == ("viewer-label-BG009368",)


def test_catalog_normalization_does_not_invent_conflicts():
    from osteosarc import SampleClaim
    from osteosarc.catalog import BUCKET, build_assets
    from osteosarc.parsing import Table
    cite, rna = "kamil/blood/Nov2025_CITE/possorted_genome_bam.bam", "vendor/cegat/P2/P3.bam"
    listing = dict(files=[[cite, 1, 0], [rna, 1, 0], ["vendor/cegat/P2/P2.1.fastq.gz", 1, 0]])
    bams = dict(baseUrl=BUCKET, genome="hg38", categories=[dict(name="Blood scRNA", bams=[
        dict(name="Hudson Lab Blood CITE 2025-11 (Pool 1)", url=BUCKET + cite, tissue="blood")])])
    metadata = Table([dict(s3_path=BUCKET + cite, display_name="CITE pool", assay="CITE",
                           timepoint="", sample_date="2025-11-06", tissue="Blood", provider="Hudson Lab")])
    vafs = Table([], columns=("bam_file",))
    # A directory row for WES FASTQs must not override the RNA BAM's own row.
    paths = (("vendor/cegat/P2", SampleClaim("data_page", "WES FASTQ", assay="wes")),
             (rna, SampleClaim("data_page", "RNA BAM", assay="rna-seq")))
    assets = {a.key: a for a in build_assets(listing, bams, metadata, vafs, paths)}
    assert assets[cite].values("assay") == ("cite-seq",)
    assert assets[cite].values("provider") == ("Hudson Lab",)
    assert assets[cite].values("date") == ("2025-11-06",)  # month precision agrees
    assert not assets[cite].conflicts
    assert assets[rna].values("assay") == ("rna-seq",)
    assert assets["vendor/cegat/P2/P2.1.fastq.gz"].values("assay") == ("wes",)


def test_reopen_is_offline_and_uses_pinned_receipts(dataset):
    reopened = Dataset.open("fixture", cache=dataset.cache)
    assert reopened.cache.offline
    assert reopened.id == dataset.id
    assert reopened.variants().to_records() == dataset.variants().to_records()
    path = reopened.download("vafs")
    assert len(reopened.table("vafs")) > 0
    assert all(isinstance(r["alt_reads"], str) for r in reopened.table("vafs"))
    assert Counter(v.status for v in reopened.variants()) == Counter(v.status for v in dataset.variants())
    path.write_text("corrupt")
    with pytest.raises(IntegrityError):
        Dataset.open("fixture", cache=dataset.cache)


def test_snapshot_and_cli_use_same_selection(dataset, capsys):
    assert main(["--cache", str(dataset.cache.root), "variants", "fixture", "--gene", "SMC5"]) == 0
    assert json.loads(capsys.readouterr().out)[0]["id"] == dataset.variants(gene="SMC5")[0].id
    assert main(["--cache", str(dataset.cache.root), "assets", "fixture", "--kind", "alignment", "--limit", "1"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["total"] == 3
    assert len(result["assets"]) == 1
    with pytest.raises(FileExistsError):
        Dataset.sync("fixture", cache=dataset.cache, refresh=True)


def test_acquired_data_remains_pinned_after_url_refresh(dataset, tmp_path):
    from osteosarc import Asset
    url = "https://example.test/results.tsv"
    path = tmp_path / "results.tsv"
    path.write_text("id\tcount\na\t0\n")
    asset = Asset(stable_id(url), "results.tsv", url, "table", "tsv")
    old = dataset.cache.import_file(path, url)
    pinned = dataset.download(asset)
    path.write_text("id\tcount\na\t9\n")
    new = dataset.cache.import_file(path, url)
    assert old.sha256 != new.sha256
    assert dataset.download(asset) == pinned
    assert dataset.parse(asset).rows[0]["count"] == "0"
    with pytest.raises(ValueError, match="new snapshot"):
        dataset.download(asset, refresh=True)


def test_unsupported_parser_refuses_before_downloading_alignment(dataset, monkeypatch):
    asset = dataset.assets.select(format="bam")[0]
    monkeypatch.setattr(dataset, "download", lambda *a: pytest.fail("downloaded unsupported format"))
    with pytest.raises(ValueError, match="No built-in parser"):
        dataset.parse(asset)


def test_vcf_annotations_and_indexed_subsetting(dataset, tmp_path):
    import pysam

    from osteosarc import Asset, Assets
    source = tmp_path / "calls.vcf"
    source.write_text('##fileformat=VCFv4.2\n'
                      '##contig=<ID=chr1,length=1000>\n'
                      '##INFO=<ID=GENE,Number=1,Type=String,Description="gene">\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n'
                      'chr1\t11\tfirst\tA\tC,G\t.\tPASS\tGENE=ONE\tGT\t1/2\n'
                      'chr1\t50\tsecond\tT\tG\t.\tPASS\tGENE=TWO\tGT\t0/1\n')
    compressed = tmp_path / "calls.vcf.gz"
    pysam.tabix_compress(str(source), str(compressed))
    pysam.tabix_index(str(compressed), preset="vcf")
    url = "https://example.test/calls.vcf.gz"
    index_url = url + ".tbi"
    asset = Asset(stable_id(url), "calls.vcf.gz", url, "variants", "vcf", index_urls=(index_url,))
    index = Asset(stable_id(index_url), "calls.vcf.gz.tbi", index_url, "index", "tbi")
    dataset.assets = Assets([asset, index])
    dataset.cache.import_file(compressed, url)
    dataset.cache.import_file(str(compressed) + ".tbi", index_url)
    with dataset.open_variants(asset) as calls:
        rows = list(calls.fetch("chr1", 10, 11))
    assert len(rows) == 1
    assert rows[0].alts == ("C", "G")
    assert rows[0].info["GENE"] == "ONE"
    assert rows[0].samples["tumor"]["GT"] == (1, 2)
