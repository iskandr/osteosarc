"""Make small public-source excerpts from a locally acquired metadata directory.

No network. Usage: python scripts/build_test_fixtures.py .cache/audit
The input filenames match docs/testing.md. Each source digest is recorded.
"""

import csv
import hashlib
import io
import json
import sys
from pathlib import Path

from bs4 import BeautifulSoup

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # run from a checkout
from osteosarc.timeline import normalize_date  # noqa: E402


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
    save("vafs.tsv", tsv(read("vafs.tsv"), lambda row: row["variant_id"] in ids))
    save("bam-metadata.tsv", tsv(read("bam-metadata.tsv"), lambda row: row["s3_path"] in keys))
    save("vafs-columns.tsv", tsv(read("vafs-columns.tsv"), lambda row: True))
    soup = BeautifulSoup(read("data.html"), "html.parser")
    save("data.html", "<h2>Bulk RNA Sequencing</h2>" + next(str(table) for table in soup.find_all("table")
         if "rna-seq/fastq/bostongene_2022" in table.get_text()))
    timeline_excerpts(read, save)
    (destination / "provenance.json").write_text(json.dumps(dict(
        accessed="2026-09-18", selection="Five named variants, three reprocessed bulk RNA products, and "
        "timeline rows around the four timepoints; parsed excerpts, not byte-identical source files",
        files=receipts), indent=2) + "\n")


def tsv(text, keep, delimiter="\t"):
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter=delimiter, lineterminator="\n")
    writer.writeheader()
    writer.writerows(row for row in reader if keep(row))
    return out.getvalue()


def timeline_excerpts(read, save):
    """Rows around the timepoints, including every record a timeline correction targets."""
    events = json.loads(read("events.json"))
    titles = {"T0", "T1", "T2", "T3", "Tempus xT", "Tempus xE", "Tempus xR", "Apheresis", "Trabectedin",
              "SQ3370 (Doxorubicin Biogel click)", "Proton therapy", "mRNA neoantigen vaccine", "ReyaGel"}
    end = events["date_range"]["end"]
    events["events"] = [e for e in events["events"] if e["title"] in titles or e.get("end_date") == end
                        or e["date"] in ("2024-06-06", "2024-06-11", "2025-01-28", "2025-04-17")]
    save("events.json", events)
    kept = {(e["title"], e["date"]) for e in events["events"]}

    save("timeline.csv", tsv(read("timeline.csv"), delimiter=",",
                             keep=lambda row: (row["Title"], normalize_date(row["Start date"])) in kept))
    mrd = json.loads(read("mrd.json"))
    mrd["measurements"] = [m for m in mrd["measurements"] if m["date"] in ("2023-01-30", "2024-06-11", "2024-11-19")
                           or isinstance(m["value"], dict)][:12]
    save("mrd.json", mrd)
    specimens = {"T1_tumor", "T2_tumor", "T3_tumor", "T3_tumor_CD45neg", "blood_2025-06-26",
                 "blood_2025-07-24", "blood_2025-08-21", "blood_2025-09-22"}
    save("samples-consolidated.tsv", tsv(read("samples-consolidated.tsv"), lambda r: r["sample_id"] in specimens))
    save("samples.json", json.loads(read("samples.json")))
    save("fastqs-consolidated.tsv", tsv(read("fastqs-consolidated.tsv"),
                                        lambda r: r["sample_id"] in ("T1_tumor", "T2_tumor")))
    flow = json.loads(read("flow-manifest.json"))
    save("flow-manifest.json", dict(date_source=flow["date_source"], samples=[
        s for s in flow["samples"] if s["draw_date"] in ("2025-06-24", "2025-07-22", "2026-06-22")]))
    imaging = json.loads(read("dicom-studies.json"))
    save("dicom-studies.json", dict(generated=imaging["generated"], studies=[
        {k: v for k, v in s.items() if k != "series"} for s in imaging["studies"]
        if s["studyDate"] in ("20240611", "20250128")]))
    pathology = json.loads(read("pathology-slides.json"))
    save("pathology-slides.json", dict(groups=pathology["groups"], slides=[]))
    save("lab_results.tsv", tsv(read("lab_results.tsv"), lambda r: r["date"] in ("2024-06-05", "2025-01-20")))
    save("cytometry.tsv", tsv(read("cytometry.tsv"), lambda r: r["date"] == "2024-07-15"))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
