import json
from copy import deepcopy

import pytest

from osteosarc import SchemaError, parse_file, parse_table, parse_variants
from osteosarc.catalog import build_assets, parse_data_paths


def test_missing_values_zero_and_ragged_tables(tmp_path):
    table = parse_table('id\tdepth\tpeptide\na\t0\t"line1\nline2"\nb\tNA\t\n')
    assert table.rows[0] == dict(id="a", depth="0", peptide="line1\nline2")
    assert table.rows[1] == dict(id="b", depth="NA", peptide="")
    assert len(table.select(depth="0")) == 1
    with pytest.raises(KeyError):
        table.select(unknown="x")
    for text in ("id\tid\na\tb", "id\tdepth\na", "id\na\tb"):
        with pytest.raises(SchemaError):
            parse_table(text)
    path = tmp_path / "empty.tsv"
    path.write_text("id\tdepth\n")
    empty = parse_file(path)
    assert len(empty) == 0 and empty.columns == ("id", "depth")


def test_conflicting_alleles_and_unknown_fields_are_retained(dataset):
    rows = list(dataset.vafs)
    original = next(r for r in rows if r["variant_id"] == "SMC5-chr9-70298024")
    rows.append(dict(original, alt="T"))
    from osteosarc import Table
    variants = parse_variants(dataset.source_path("variant_index").read_text(), Table(rows))
    item = variants["SMC5-chr9-70298024"]
    assert item.status == "ambiguous_literal_allele"
    assert len(item.alleles) == 2
    with pytest.raises(ValueError):
        item.region()


def test_vaccine_join_does_not_assign_ambiguous_loci(dataset):
    rows = list(dataset.vafs)
    original = next(r for r in rows if r["variant_id"] == "SMC5-chr9-70298024")
    rows.append(dict(original, variant_id="another-SMC5-allele", alt="T"))
    from osteosarc import Table
    variants = parse_variants(dataset.source_path("variant_index").read_text(), Table(rows),
                              vaccine_overlap=json.loads(dataset.source_path("vaccine_overlap").read_text()))
    assert not variants["SMC5-chr9-70298024"].vaccines
    assert variants["SMC5-chr9-70298024"].annotations["ambiguous_vaccine_join"]


def test_all_objects_and_ambiguous_basenames_remain_visible(dataset):
    listing = dict(download_base="https://example.test/", files=[
        ["a/same.bam", 1, 1], ["b/same.bam", 2, 1], ["a/same.bam.bai", 3, 1],
        ["new+folder/raw.fastq.gz", 4, 1], ["unknown.extension", 5, 1]])
    from osteosarc import Table
    vafs = Table([dict(bam_file="same.bam", sample_label="T0", timepoint="T0")])
    assets = build_assets(listing, dict(baseUrl="https://example.test/", categories=[]), [], vafs)
    bams = assets.select(kind="alignment")
    assert len(bams) == 2
    assert not bams.select(timepoint="T0")
    assert all(a.metadata["ambiguous_vaf_basename"] == "same.bam" for a in bams)
    assert bams[0].index_urls == ("https://example.test/a/same.bam.bai",)
    assert assets.select(kind="reads")[0].url == "https://example.test/new%2Bfolder/raw.fastq.gz"
    assert assets.select(kind="other")[0].key == "unknown.extension"
    duplicate = deepcopy(listing)
    duplicate["files"].append(duplicate["files"][0])
    with pytest.raises(SchemaError):
        build_assets(duplicate, dict(categories=[]), [], vafs)


def test_data_page_links_enrich_raw_files_without_provider_guessing(dataset):
    claims = parse_data_paths(dataset.source_path("data_page").read_text())
    claim = next(c for prefix, c in claims if prefix == "rna-seq/fastq/bostongene_2022")
    assert claim.timepoint == "T0" and claim.assay == "rna-seq"
    assert claim.platform is None
