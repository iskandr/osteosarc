import json
from copy import deepcopy

import pytest

from osteosarc import SchemaError, parse_variants
from osteosarc.catalog import build_files, parse_data_paths
from osteosarc.parsing import parse_file, parse_table


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
    assets = build_files(listing, dict(baseUrl="https://example.test/", categories=[]), [], vafs)
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
        build_files(duplicate, dict(categories=[]), [], vafs)


def test_data_page_links_enrich_raw_files_without_provider_guessing(dataset):
    claims = parse_data_paths(dataset.source_path("data_page").read_text())
    claim = next(c for prefix, c in claims if prefix == "rna-seq/fastq/bostongene_2022")
    assert claim.timepoint == "T0" and claim.assay == "rna-seq"
    assert claim.platform is None
def test_unsupported_format_does_not_read_large_file(tmp_path):
    from osteosarc.parsing import parse_file
    # The nonexistent path must not be opened before deciding parser support.
    with pytest.raises(ValueError, match="No built-in parser"):
        parse_file(tmp_path / "huge.bam")


@pytest.mark.parametrize("bad", [
    "bad\tBAD\tchr1\t10\n",
    "bad\tBAD\tchr1\t10\tA\tC\textra\n",
    "bad\tBAD\tchr1\tNA\tA\tC\n",
    "bad\tBAD\tchr1\t\tA\tC\n",
    "bad\tBAD\tchr1\t10.0\tA\tC\n",
])
@pytest.mark.parametrize("bad_first", [True, False])
def test_malformed_catalogue_rows_preserve_each_entry(dataset, tmp_path, bad, bad_first):
    from osteosarc import Dataset

    index = [dict(id="bad", gene="BAD", location="chr1:10", vaccine_count=0),
             dict(id="good", gene="GOOD", location="chr2:20", vaccine_count=0)]
    good = "good\tGOOD\tchr2\t20\tG\tT\n"
    text = "variant_id\tgene\tchrom\tpos\tref\talt\n" + (bad + good if bad_first else good + bad)

    def check(variants):
        assert variants["good"].allele == ("chr2", 20, "G", "T")
        unavailable = variants["bad"]
        assert unavailable.status == "malformed_source_row" and not unavailable.alleles
        error, = unavailable.annotations["parse_errors"]
        assert error["source"] == "vafs" and error["row"] == (0 if bad_first else 1)
        assert error["values"]["pos"] == bad.split("\t")[3].rstrip("\n")
        if error["code"] == "ragged_row":
            assert error["fields"] == tuple(bad.rstrip("\n").split("\t"))
        with pytest.raises(ValueError, match="malformed_source_row"):
            unavailable.region()

    # Text and explicit tolerant Table are both supported public entry points.
    check(parse_variants(index, text))
    table = parse_table(text, strict=False)
    check(parse_variants(index, table))
    bad_only = table.select(variant_id="bad")
    if table.diagnostics:
        assert bad_only.diagnostics[0]["row"] == 0
        assert not table.select(variant_id="good").diagnostics

    original = tmp_path / "malformed.tsv"
    original.write_text(text)
    dataset.cache.import_file(original, dataset.manifest["sources"]["vafs"]["url"])
    # Use the actual sync/reopen/correction paths. These two export-only IDs
    # must survive alongside the fixture's ordinary site entries.
    data = Dataset.sync("malformed", cache=dataset.cache)
    check(data.variants("all"))
    bad_row = data.vafs.select(variant_id="bad").rows[0]
    assert bad_row["ref"] == (None if len(bad.rstrip("\n").split("\t")) == 4 else "A")
    assert data.vafs.diagnostics == table.diagnostics
    check(Dataset.open("malformed", cache=data.cache, corrections=False).variants("all"))


@pytest.mark.parametrize("header", [
    "variant_id\tgene\tchrom\tpos\tref",  # missing ALT
    "variant_id\tgene\tchrom\tpos\tref\talt\talt",
])
def test_malformed_catalogue_headers_still_fail(header):
    with pytest.raises(SchemaError):
        parse_variants([], header + "\n")


def test_malformed_row_for_otherwise_valid_entry_stays_nonready():
    text = ("variant_id\tgene\tchrom\tpos\tref\talt\n"
            "entry\tGENE\tchr1\t10\tA\tC\n"
            "entry\tGENE\tchr1\tNA\tA\tC\n")
    variants = parse_variants([], text)
    assert variants["entry"].alleles == (("chr1", 10, "A", "C"),)
    assert variants["entry"].status == "malformed_source_row"
    assert not variants.select(status="ready")
