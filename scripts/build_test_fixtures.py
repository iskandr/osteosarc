"""Make small public-source excerpts from a locally acquired metadata directory.

No network. Usage: python scripts/build_test_fixtures.py .cache/audit
The input filenames match docs/validation.md. Each source digest is recorded.
"""

import csv
import hashlib
import io
import json
import sys
from pathlib import Path

from bs4 import BeautifulSoup


def main(source):
    destination = Path(__file__).resolve().parents[1] / "tests" / "data"
    destination.mkdir(parents=True, exist_ok=True)
    receipts = {}

    def read(name):
        data = (source / name).read_bytes()
        receipts[name] = dict(source_sha256=hashlib.sha256(data).hexdigest())
        return data.decode()

    def save(name, value):
        if not isinstance(value, str):
            value = json.dumps(value, indent=2) + "\n"
        path = destination / name
        path.write_text(value)
        receipts[name]["fixture_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()

    ids = {"DYNC1H1-chr14-101980529", "SMC5-chr9-70298024", "COL3A1-Splice",
           "GTF3C5-chr9-133057893", "FAM157A-p_W70_Q71ins_14"}
    soup = BeautifulSoup(read("variants.html"), "html.parser")
    rows = [str(row) for row in soup.select("tr[data-vaccines]")
            if row.select_one('a[href^="/variant/"]')["href"].rstrip("/").split("/")[-1] in ids]
    save("variants.html", "<table><tbody>" + "\n".join(rows) + "</tbody></table>")
    variants = json.loads(read("source-variants.json"))
    save("source-variants.json", [{k: v for k, v in row.items() if k != "vafs"}
                                  for row in variants if row["id"] in ids])
    vaccines = json.loads(read("vaccine_overlap.json"))
    vaccines["mutations"] = [r for r in vaccines["mutations"] if r["gene"] in ("SMC5", "DYNC1H1", "GTF3C5")]
    save("vaccine_overlap.json", vaccines)
    bams = json.loads(read("bams.json"))
    categories = []
    for category in bams["categories"]:
        selected = [row for row in category["bams"] if any(token in row["url"] for token in
                    ("BG003082", "BG009368", "SARC0277"))]
        if selected:
            categories.append(dict(category, bams=selected))
    bams["categories"] = categories
    save("bams.json", bams)
    keys = {row["url"] for c in categories for row in c["bams"]}
    listing = json.loads(read("bucket_listing.json"))
    listing["files"] = [row for row in listing["files"] if row[0] in keys
                         or row[0] in {key + ".bai" for key in keys}]
    save("bucket_listing.json", listing)
    for name in ("vafs.tsv", "bam-metadata.tsv", "vafs-columns.tsv"):
        reader = csv.DictReader(io.StringIO(read(name)), delimiter="\t")
        rows = list(reader)
        if name == "vafs.tsv":
            rows = [row for row in rows if row["variant_id"] in ids]
        elif name == "bam-metadata.tsv":
            rows = [row for row in rows if row["s3_path"] in keys]
        out = io.StringIO()
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        save(name, out.getvalue())
    soup = BeautifulSoup(read("data.html"), "html.parser")
    save("data.html", "<h2>Bulk RNA Sequencing</h2>" + next(str(table) for table in soup.find_all("table")
         if "rna-seq/fastq/bostongene_2022" in table.get_text()))
    (destination / "provenance.json").write_text(json.dumps(dict(
        accessed="2026-09-18", selection="Five named variants, three reprocessed bulk RNA products; parsed excerpts, not byte-identical source files", files=receipts), indent=2) + "\n")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
