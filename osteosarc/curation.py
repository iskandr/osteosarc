"""Every hand-written interpretation of the osteosarc sources, in one place.

Two layers live here:

* Vocabulary (always on): maps source spellings to the names used by
  selection filters, such as "Boston Gene" -> "BostonGene" or a "CITE" label
  -> cite-seq. Mapping is lossless: claims keep their original labels and
  metadata rows keep their original values. Values outside this vocabulary
  are reported by ``unrecognized_values`` instead of being silently dropped.

* Corrections (optional): specific source records that are known to be wrong
  or misleading. Each correction names the records it touches, the published
  values it was written against, what it replaces (or that it only flags the
  records), and its evidence. ``Dataset.open(name, corrections=False)`` uses
  the published sources unchanged; ``corrections=[...]`` supplies your own.

Corrections anticipate upstream edits. Every load re-checks each correction
against the snapshot's sources:

``applied``         expectations hold; the change was made (or the flag attached)
``fixed_upstream``  corrected values or a reviewed removal; nothing to do
``stale``           the source changed some other way or the record is gone;
                    the correction is NOT applied and a CurationWarning is issued
``unavailable``     the snapshot predates a source the correction needs
``disabled``        corrections were turned off for this Dataset

A correction applies atomically: if any of its changes is stale, none is made.
"""

from __future__ import annotations

import copy
import fnmatch
import functools
import json
import re
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from importlib.resources import files

from .urls import BUCKET, SOURCE_REPO

# ---------------------------------------------------------------------------
# Vocabulary (always on)
# ---------------------------------------------------------------------------

#: Source assay labels (viewer categories, consolidated/VAF assay columns).
ASSAYS = {"RNA": ("rna-seq", None), "WGS": ("wgs", None), "WES": ("wes", None),
          "scRNA ONT": ("scrna-seq", "ont"), "scRNA_ONT": ("scrna-seq", "ont"),
          "PacBio": ("scrna-seq", "pacbio"), "scRNA": ("scrna-seq", None),
          "Tumor scRNA": ("scrna-seq", None), "Blood scRNA": ("scrna-seq", None),
          "CITE": ("cite-seq", None)}

#: Viewer labels that name a more specific assay than their category.
#: CITE-seq libraries are filed under the viewer's "Blood scRNA" category.
LABEL_ASSAYS = {"CITE": "cite-seq"}

#: Provider names recognized in free-text labels and provider columns.
PROVIDERS = ("BostonGene", "CeGaT", "Hudson Lab", "Natera", "Personalis", "Tempus", "UCLA", "UCSF")
PROVIDER_ALIASES = {"Boston Gene": "BostonGene"}

#: "T1-organoid" is timepoint T1; the organoid specimen is recorded as tissue.
TIMEPOINT_ALIASES = {"T1-organoid": "T1"}

#: Every normal specimen in this dataset is a blood normal, and the sources use
#: "normal" and "blood" interchangeably, so both map to "blood".
#: A composite "Tumor + Normal" row stays unresolved rather than guessed.
TISSUES = {"tumor": "tumor", "blood": "blood", "normal": "blood", "normal (blood)": "blood",
           "blood normal": "blood", "blood/normal": "blood", "organoid": "organoid",
           "tumor + normal": None}

#: Detection and vaccine keys present when the curation below was written.
#: New keys are not errors, but they are reported as possible upstream drift.
PIPELINES = ("BG 2024", "DRAGEN", "LENS 2022", "LENS 2024", "Mutect2 2024", "Mutect2 2025",
             "Natera 2022", "Tempus 2022", "oncoanalyser", "pVACtools 2025")
VACCINES = ("CeGaT", "Cure 2024", "JLF V1", "JLF V2", "JLF V3", "mRNA")


def normalize_provider(value):
    value = PROVIDER_ALIASES.get(value, value)
    return value or None


def normalize_timepoint(value):
    return TIMEPOINT_ALIASES.get(value, value) or None


def normalize_tissue(value):
    """Map a tissue label; unknown labels are kept (lowercased), not guessed."""
    value = (value or "").strip().lower()
    return TISSUES[value] if value in TISSUES else value or None


def label_assay(label, default):
    """Refine a viewer category's assay when its label names a narrower assay."""
    for token, assay in LABEL_ASSAYS.items():
        if re.search(rf"\b{re.escape(token)}\b", label):
            return assay
    return default


def label_claims(label):
    """Date and provider stated in a viewer label such as 'T0 BostonGene Tumor RNA 2022-12'.

    A provider is taken only when exactly one known provider name appears.
    """
    dates = re.findall(r"\b(\d{4}-\d{2}(?:-\d{2})?)\b", label)
    providers = [p for p in PROVIDERS if re.search(rf"\b{p}\b", label)]
    return dict(date=dates[0] if len(dates) == 1 else None,
                provider=providers[0] if len(providers) == 1 else None)


def unrecognized_values(*, bams=None, metadata=(), vafs=(), data_page=(), source_variants=()):
    """Source values outside the vocabulary above, with record counts.

    An empty result means every assay, tissue, provider, timepoint, pipeline and
    vaccine label is one this module interprets. New values usually mean the
    website changed; review them before trusting filters on those fields.
    """
    seen = Counter()
    known_assays, known_tissues = set(ASSAYS), set(TISSUES)
    known_providers = set(PROVIDERS) | set(PROVIDER_ALIASES)

    def check(source, name, value, known):
        if value and value not in known:
            seen[(source, name, value)] += 1

    def timepoint(source, value):
        if value and not re.fullmatch(r"T\d+", value) and value not in TIMEPOINT_ALIASES:
            seen[(source, "timepoint", value)] += 1

    for category in (bams or {}).get("categories", ()):
        check("bams", "category", category["name"], known_assays)
        for row in category["bams"]:
            check("bams", "tissue", (row.get("tissue") or "").lower(), known_tissues)
    for row in metadata:
        check("bam_metadata", "assay", row.get("assay"), known_assays)
        check("bam_metadata", "tissue", (row.get("tissue") or "").lower(), known_tissues)
        check("bam_metadata", "provider", row.get("provider"), known_providers)
        timepoint("bam_metadata", row.get("timepoint"))
    for row in vafs:
        check("vafs", "assay_type", row.get("assay_type"), known_assays)
        check("vafs", "tissue", (row.get("tissue") or "").lower(), known_tissues)
        check("vafs", "data_source", row.get("data_source"), known_providers)
        timepoint("vafs", row.get("timepoint"))
    for values in data_page:
        check("data_page", "Tissue", (values.get("Tissue") or "").lower(), known_tissues)
    for record in source_variants:
        for key in record.get("detection", {}):
            check("source_variants", "detection", key, set(PIPELINES))
        for key in record.get("vaccines", {}):
            check("source_variants", "vaccines", key, set(VACCINES))
    return [dict(source=s, field=f, value=v, records=n) for (s, f, v), n in sorted(seen.items())]


# ---------------------------------------------------------------------------
# Corrections (optional)
# ---------------------------------------------------------------------------

class CurationWarning(UserWarning):
    """A correction no longer matches its source and was not applied."""


@dataclass(frozen=True)
class Glob:
    """Match a string field with shell-style wildcards (``*`` spans ``/``)."""

    pattern: str

    def matches(self, value):
        return (isinstance(value, str) and value.startswith(_literal_prefix(self.pattern))
                and fnmatch.fnmatchcase(value, self.pattern))


@functools.lru_cache(maxsize=None)
def _literal_prefix(pattern):
    return re.split(r"[*?\[]", pattern, maxsplit=1)[0]


def glob(pattern):
    return Glob(pattern)


@dataclass(frozen=True)
class Change:
    """One edit to the records of one snapshot source.

    ``match`` selects records by field equality (or ``glob``). ``expect`` lists
    published values the correction was written against. ``set`` gives the
    replacement values; leave it empty to flag the records without editing.
    Dotted field names reach into nested objects ("detection.pVACtools 2025").
    ``absent=True`` requires that no records match, for a reviewed removal.
    """

    source: str
    match: dict
    expect: dict = field(default_factory=dict)
    set: dict = field(default_factory=dict)
    absent: bool = False

    def __post_init__(self):
        if self.absent and (self.expect or self.set):
            raise ValueError("An absence check cannot expect or set record values")


@dataclass(frozen=True)
class Correction:
    """An atomic group of changes, with optional reviewed source layouts.

    Exactly one of ``changes`` and the groups in ``alternatives`` must match.
    """

    id: str
    summary: str
    changes: tuple[Change, ...]
    evidence: tuple[str, ...] = ()
    verified: str = ""
    alternatives: tuple[tuple[Change, ...], ...] = ()

    @property
    def versions(self):
        return (self.changes, *self.alternatives)


def _get(record, name):
    value = record
    for part in name.split("."):
        if not isinstance(value, dict) or part not in value:
            return _MISSING
        value = value[part]
    return value


def _set(record, name, value):
    *parents, last = name.split(".")
    for part in parents:
        record = record.setdefault(part, {})
    record[last] = value


_MISSING = object()


def _matches(record, match):
    for name, wanted in match.items():
        value = _get(record, name)
        if isinstance(wanted, Glob):
            if not wanted.matches(value):
                return False
        elif value != wanted:
            return False
    return True


class Curation:
    """Evaluate corrections against a snapshot's parsed source records.

    ``load(source)`` returns the published records (a list of dicts) for a
    source name. Records are never modified in place.
    """

    def __init__(self, corrections, load, *, enabled=True):
        self.corrections = tuple(corrections)
        ids = [c.id for c in self.corrections]
        if len(set(ids)) != len(ids):
            raise ValueError("Correction IDs must be unique")
        self.enabled = enabled
        self._load = load
        self._raw = {}
        self._evaluated = {}
        self._matched = {}
        self._selected = {}
        self._records = {}

    def raw(self, source):
        if source not in self._raw:
            self._raw[source] = self._load(source)
        return self._raw[source]

    def evaluate(self, correction):
        """Return (status, per-change details) for one correction."""
        if correction.id in self._evaluated:
            return self._evaluated[correction.id]
        candidates = [self._evaluate_changes(changes) for changes in correction.versions]
        valid = [i for i, (status, _, _) in enumerate(candidates)
                 if status in ("applied", "fixed_upstream")]
        # Prefer the matching published layout for diagnostics as well as edits.
        selected = valid[0] if valid else max(range(len(candidates)), key=lambda i:
            sum(d["records"] for d in candidates[i][1]))
        status, details, matched = candidates[selected]
        if len(valid) > 1:
            status = "stale"
            details = [*details, dict(source="correction", match=correction.id,
                                     state="ambiguous_versions", records=0, differing=[])]
        self._selected[correction.id] = correction.versions[selected]
        for n, indices in enumerate(matched):
            self._matched[correction.id, n] = indices
        if status == "stale" and self.enabled:
            problems = "; ".join(f"{d['source']} {d['match']}: {d['state']}"
                                 + (f" ({', '.join(d['differing'])})" if d["differing"] else "")
                                 for d in details if d["state"] in
                                 ("missing", "changed", "unexpected", "ambiguous_versions"))
            warnings.warn(f"Correction {correction.id!r} was not applied because its source changed: "
                          f"{problems}. Review data.corrections.", CurationWarning, stacklevel=3)
        self._evaluated[correction.id] = status, details
        return status, details

    def _evaluate_changes(self, changes):
        details, matches = [], []
        for change in changes:
            records = self.raw(change.source)
            if records is None:  # the snapshot predates this source
                details.append(dict(source=change.source, match=_describe(change.match),
                                    state="unavailable", records=0, differing=[]))
                matches.append([])
                continue
            matched = [i for i, record in enumerate(records) if _matches(record, change.match)]
            matches.append(matched)
            differing = sorted({name for i in matched for name, value in change.expect.items()
                                if _get(records[i], name) != value})
            if change.absent:
                state = "unexpected" if matched else "absent"
            elif not matched:
                state = "missing"
            elif (change.set and not any(k not in change.set for k in differing)
                  and all(_get(records[i], k) == v for i in matched for k, v in change.set.items())):
                state = "already_correct"
            elif differing:
                state = "changed"
            else:
                state = "pending"
            details.append(dict(source=change.source, match=_describe(change.match), state=state,
                                records=len(matched), differing=differing))
        states = [d["state"] for d in details]
        if "unavailable" in states:
            return "unavailable", details, matches
        edits = [d["state"] for d, c in zip(details, changes) if c.set]
        if not states or any(s not in ("pending", "already_correct", "absent") for s in states):
            status = "stale"
        elif not edits and "absent" in states:
            status = "fixed_upstream"
        elif edits and all(s == "already_correct" for s in edits):
            status = "fixed_upstream"
        else:
            status = "applied"
        return status, details, matches

    def records(self, source):
        """Corrected records plus, per record index, the IDs of corrections touching it."""
        if source in self._records:
            return self._records[source]
        records, touched = list(self.raw(source) or ()), defaultdict(list)
        if self.enabled:
            for correction in self.corrections:
                if not any(c.source == source for group in correction.versions for c in group):
                    continue
                if self.evaluate(correction)[0] != "applied":
                    continue
                for n, change in enumerate(self._selected[correction.id]):
                    if change.source != source:
                        continue
                    # Matches are made against the published records.
                    for i in self._matched[correction.id, n]:
                        if change.set:
                            record = copy.deepcopy(records[i])
                            for name, value in change.set.items():
                                _set(record, name, value)
                            records[i] = record
                        if correction.id not in touched[i]:
                            touched[i].append(correction.id)
        self._records[source] = records, {i: tuple(ids) for i, ids in touched.items()}
        return self._records[source]

    def release(self, source):
        """Drop cached records of a large source; they are reloaded if needed again."""
        self._raw.pop(source, None)
        self._records.pop(source, None)

    def applied(self):
        return tuple(c.id for c in self.corrections if self.enabled and self.evaluate(c)[0] == "applied")

    def report(self):
        rows = []
        for correction in self.corrections:
            status, details = self.evaluate(correction)
            rows.append(dict(id=correction.id, status=status if self.enabled else "disabled",
                             evaluation=status, summary=correction.summary,
                             action="edit" if any(c.set for c in self._selected[correction.id]) else "flag",
                             changes=details, evidence=list(correction.evidence),
                             verified=correction.verified))
        return rows


def _describe(match):
    return ", ".join(f"{k}={v.pattern if isinstance(v, Glob) else v}" for k, v in match.items())


# ---------------------------------------------------------------------------
# The correction registry
# ---------------------------------------------------------------------------

_SITE_REPO, _BUCKET = SOURCE_REPO, BUCKET
_CONSOLIDATED = _SITE_REPO + "scripts/data/bam-metadata-consolidated.tsv"
_PVAC = "neoantigen_prediction/pvactools/"


_ALMY = _BUCKET + "vendor/tempus/TL-24-ALMY2X4KMV/DNA/TL-24-ALMY2X4KMV.soma."
_COUNTS = dict.fromkeys(("ref_reads", "alt_reads", "other_reads", "total_reads", "vaf"), "")
_RESOLUTIONS = json.loads(files("osteosarc").joinpath("data/allele_resolutions.json").read_text())
_COMMITS = SOURCE_REPO.split("/-/")[0] + "/-/commit/"
_COUNT_FIX = _COMMITS + "7e80e3a6432bb4c8d5c5260fdd4e18f25ba5b9fb"
_COORDINATE_FIX = _COMMITS + "cf2f59169aa65e851269544f78bbfdf4520f896a"
_PILEUP_FIX = _COMMITS + "f6c59f0fc01201c2158b4503fe1b821fbbb5b415"
_USH2A_FIX = _COMMITS + "bc138897b77323f595acffa46cbfce5c45720a1d"


def _absent_variant(variant_id):
    return tuple(Change(source, {key: variant_id}, absent=True)
                 for source, key in (("source_variants", "id"), ("variant_index", "id"),
                                     ("vafs", "variant_id")))


def _kind(ref, alt):
    return "insertion" if len(alt) > len(ref) else "deletion" if len(alt) < len(ref) else (
        "snv" if len(ref) == 1 else "mnv")


def _catalogue_resolution(variant_id, correction_id):
    """Attach reviewed evidence; fill missing alleles without inventing count rows."""
    entry = _RESOLUTIONS[variant_id]
    resolution = entry["resolution"]
    values, index_values = {"allele_resolution": resolution}, {}
    if resolution["status"] == "resolved":
        allele = resolution["allele"]
        chrom, pos, ref, alt = (allele[k] for k in ("chrom", "pos", "ref", "alt"))
        values.update(chr=chrom, pos=pos, ref=ref, alt=alt,
                      variant_type=_kind(ref, alt), genomic_location=f"{chrom}:{pos}")
        index_values["location"] = f"{chrom}:{pos}"
    changes = [
        Change("source_variants", {"id": variant_id}, expect=entry["expected_source"], set=values),
        Change("variant_index", {"id": variant_id}, expect=entry["expected_index"], set=index_values),
    ]
    if entry["expected_source"]["pos"] is not None:
        old = entry["expected_source"]
        changes.append(Change("vafs", {"variant_id": variant_id},
                              expect=dict(chrom=old["chr"], pos=str(old["pos"]),
                                          ref=old["ref"], alt=old["alt"])))
    if "expected_related_source" in entry:
        changes.append(Change("source_variants", {"id": resolution["relationship"]["variant_id"]},
                              expect=entry["expected_related_source"]))
    evidence = [resolution["source"]["url"]]
    evidence.extend(r["url"] for r in resolution.get("references", ()))
    if "mapping" in resolution:
        evidence.append(resolution["mapping"]["url"])
    return Correction(correction_id, resolution["summary"], tuple(changes),
                      evidence=tuple(evidence), verified=entry["verified"])


def _allele(variant_id, old, new, summary, evidence, *, extra_source=None):
    """Replace a catalogue allele everywhere it is published.

    old/new are (chrom, pos, ref, alt). The count rows were measured for the old
    allele, so their counts become missing ("") rather than being reattributed.
    """
    (chrom, pos, ref, alt), (new_chrom, new_pos, new_ref, new_alt) = old, new
    kind = _kind(new_ref, new_alt)
    changes = [
        Change("source_variants", {"id": variant_id},
               expect=dict(chr=chrom, pos=pos, ref=ref, alt=alt),
               set=dict(chr=new_chrom, pos=new_pos, ref=new_ref, alt=new_alt, variant_type=kind,
                        genomic_location=f"{new_chrom}:{new_pos}", **(extra_source or {}))),
        Change("vafs", {"variant_id": variant_id},
               expect=dict(chrom=chrom, pos=str(pos), ref=ref, alt=alt),
               set=dict(chrom=new_chrom, pos=str(new_pos), chr_pos=f"{new_chrom}:{new_pos}", ref=new_ref,
                        alt=new_alt, change=f"{new_ref}>{new_alt}", variant_type=kind, **_COUNTS)),
    ]
    if (chrom, pos) != (new_chrom, new_pos):
        changes.append(Change("variant_index", {"id": variant_id}, expect={"location": f"{chrom}:{pos}"},
                              set={"location": f"{new_chrom}:{new_pos}"}))
    return Correction(f"allele-{variant_id}", summary, tuple(changes), evidence=evidence,
                      verified="2026-09-18")


def _tempus_relocation(variant_id, old, new, grch37, caller="pindel"):
    original = _allele(
        variant_id, old, new,
        f"The catalogue places this Tempus call at {old[0]}:{old[1]} without a literal allele. "
        f"The original Tempus TL-24-ALMY2X4KMV record (GRCh37 {grch37}) lifts to "
        f"{new[0]}:{new[1]}, where its REF matches GRCh38; indel normalization cannot bridge "
        f"the gap. The website position came from a wrong transcript offset. Counts measured "
        f"at the old position are cleared.",
        (_ALMY + caller + ".vcf", f"https://rest.ensembl.org/map/human/GRCh37/{grch37}..{grch37.split(':')[1]}:1/GRCh38"))
    # September 20: coordinates/IDs were fixed, but the alleles used for the
    # regenerated counts are still placeholders. Do not reuse those counts.
    new_id = f"{variant_id.split('-')[0]}-{new[0]}-{new[1]}"
    relocated = _allele(new_id, (new[0], new[1], old[2], old[3]), new,
                        original.summary, original.evidence)
    current = (*relocated.changes,
               Change("variant_index", {"id": new_id}, expect={"location": f"{new[0]}:{new[1]}"}),
               *_absent_variant(variant_id))
    return replace(original, changes=(*original.changes, *_absent_variant(new_id)),
                   alternatives=(current,), verified="2026-09-21",
                   summary=original.summary + " The renamed upstream entry still needs its literal "
                           "allele; counts regenerated for the placeholder are also cleared.",
                   evidence=(*original.evidence, _COORDINATE_FIX, _PILEUP_FIX))


def _ush2a_resolution():
    old_id, new_id = "USH2A-chr1-215560752", "USH2A-chr1-215650752"
    original = _catalogue_resolution(old_id, "ush2a-transposed-duplicate")
    current = (
        *_absent_variant(old_id),
        Change("source_variants", {"id": new_id},
               expect=dict(gene="USH2A", chr="chr1", pos=215650752, ref="C", alt="A",
                           genomic_change_on_cdna="c.14183G>T", refseq_id="NM_206933",
                           protein_change="p.Cys4728Phe", genomic_location="chr1:215560752"),
               set=dict(genomic_location="chr1:215650752", upstream_merge=dict(
                   retired_id=old_id, retained_id=new_id, evidence=_USH2A_FIX,
                   original_source_identity="unconfirmed"))),
        Change("variant_index", {"id": new_id},
               expect=dict(gene="USH2A", location="chr1:215650752", protein_label="p.Cys4728Phe")),
        Change("vafs", {"variant_id": new_id},
               expect=dict(gene="USH2A", chrom="chr1", pos="215650752", ref="C", alt="A")),
    )
    return replace(original, alternatives=(current,), verified="2026-09-21",
                   evidence=(*original.evidence, _USH2A_FIX),
                   summary=original.summary + " After the site's merge, preserve its retained "
                           "allele/counts, record the merge provenance, and fix the old location label.")


def _stale_viewer_label(library, old, new, note):
    key = f"rna-seq/reprocessed/{library}/{library}.Aligned.sortedByCoord.out.md.bam"
    return Correction(
        f"viewer-label-{library}",
        f"The BAM viewer still labels {library} '{old}'; the consolidated metadata "
        f"deliberately re-assigned it ({note}).",
        (Change("bams", {"url": key}, expect={"name": old}, set={"name": new}),),
        evidence=(_CONSOLIDATED + " (notes column)", "https://osteosarc.com/bams/bams.json"),
        verified="2026-09-18")


CORRECTIONS = (
    Correction(
        "tempus-grch37-counts",
        "The website counted reads in the GRCh37 (b37) Tempus WES BAM TL-24-5GQLV9WSXQ at GRCh38 "
        "coordinates, so every count for this BAM describes an unrelated locus. PDZRN4's 7/7 "
        "(100% VAF) are reference reads elsewhere; at the lifted position depth is 0. Recounted "
        "at lifted positions, six empty sites have coverage (e.g. KMT2D 0/1851, AKT2 2/1217). "
        "Counts are cleared; recount from the BAM on GRCh37 coordinates if needed.",
        (Change("vafs", {"bam_file": "TL-24-5GQLV9WSXQ_T.sorted.bam"},
                expect={"sample_label": "T1 Tempus Tumor WES 2024-06 (TL-24-5GQLV9WSXQ)"}, set=_COUNTS),
         # Witness: the known-wrong value. If the site recomputes counts, this goes stale.
         Change("vafs", {"bam_file": "TL-24-5GQLV9WSXQ_T.sorted.bam", "variant_id": "PDZRN4-chr12-41572809"},
                expect={"alt_reads": "7", "total_reads": "7"}, set=_COUNTS)),
        evidence=(_BUCKET + "vendor/tempus/TL-24-5GQLV9WSXQ/DNA/TL-24-5GQLV9WSXQ_T.sorted.bam (header: "
                  "AS:human_g1k_v37)", _SITE_REPO + "scripts/variants/run_pileup-json.sh",
                  _SITE_REPO + "crates/pileup-json/src/main.rs (resolve_chrom only toggles 'chr')",
                  _COUNT_FIX),
        verified="2026-09-21",
        alternatives=((
            Change("vafs", {"bam_file": glob("*TL-24-5GQLV9WSXQ*")}, absent=True),
            Change("bams", {"url": glob("*TL-24-5GQLV9WSXQ*")}, absent=True),
            Change("bam_metadata", {"s3_path": "vendor/tempus/TL-24-5GQLV9WSXQ/DNA/"
                                                "TL-24-5GQLV9WSXQ_T.sorted.bam"},
                   expect=dict(provider="Tempus", assay="WES", timepoint="T1", tissue="Tumor")),
            Change("vafs", {"variant_id": "PDZRN4-chr12-41572809"}),
        ),)),
    _stale_viewer_label(
        "BG009368", "T0 BostonGene Tumor RNA 2022-12 (BG009368 reprocessed)",
        "T1 BostonGene Tumor RNA 2024-06 (BG009368 reprocessed)",
        "oncoanalyser T1 RNA was built from its 2024-06-11 FASTQ; 123/133 variants VAF-match"),
    _stale_viewer_label(
        "SARC0277", "T0 BostonGene Tumor RNA 2022-12 (SARC0277 reprocessed)",
        "T2 UCLA Tumor RNA 2025-01 (SARC0277 reprocessed)",
        "oncoanalyser T2 RNA was built from its 2025-01-28 FASTQ; 114/120 variants VAF-match"),
    Correction(
        "provider-IPISRC044-T1-rna",
        "The oncoanalyser T1 RNA BAM was built from BostonGene's BG009368 FASTQ (per the "
        "consolidated metadata's own note), but its consolidated row and viewer label say UCLA. "
        "The VAF export already says BostonGene.",
        (Change("bam_metadata",
                {"s3_path": "kamil/oncoanalyser/IPISRC044_T1_ucla/alignments/rna/"
                            "IPISRC044_tumor_T1_ucla_rna.md.bam"},
                expect={"provider": "UCLA"}, set={"provider": "BostonGene"}),
         Change("bams", {"url": "kamil/oncoanalyser/IPISRC044_T1_ucla/alignments/rna/"
                                "IPISRC044_tumor_T1_ucla_rna.md.bam"},
                expect={"name": "T1 UCLA Tumor RNA oncoanalyser"},
                set={"name": "T1 BostonGene Tumor RNA oncoanalyser"})),
        evidence=(_CONSOLIDATED, "https://osteosarc.com/variants/variant_vafs_long.tsv",
                  "https://osteosarc.com/bams/bams.json"),
        verified="2026-09-18"),
    Correction(
        "gene-symbol-TRMO",
        "Catalogue gene symbol TMRO is a typo for TRMO (HGNC, and the symbol pVACseq reports "
        "at the identical chr9:97910412 C>T allele). Gene-symbol joins otherwise miss it.",
        tuple(Change(source, {key: "TMRO-chr9-97910412"}, expect={"gene": "TMRO"}, set={"gene": "TRMO"})
              for source, key in (("variant_index", "id"), ("source_variants", "id"),
                                  ("vafs", "variant_id"))),
        evidence=("https://www.genenames.org/data/gene-symbol-report/#!/symbol/TRMO",
                  _BUCKET + _PVAC + "2025.04.27.sg.curated.neoantigen.predictions/MHC_Class_I/"
                                    "SG.WGS_SG.WGS.UCLA.2025.01.tumor.tsv"),
        verified="2026-09-18"),
    Correction(
        "pvac-2025-detection",
        "The pVACtools 2025 detection flag reflects only the 2025-04-25 run. The 2025-04-27 "
        "runs add exactly CDC40, PIP5K1A, SMC5 and TECPR1 at identical alleles, with "
        "hundreds of epitope rows each.",
        tuple(Change("source_variants", {"id": vid}, expect={"detection.pVACtools 2025": False},
                     set={"detection.pVACtools 2025": True})
              for vid in ("CDC40-chr6-110228838", "PIP5K1A-chr1-151242178",
                          "SMC5-chr9-70298024", "TECPR1-chr7-98241127")),
        evidence=(_BUCKET + _PVAC + "2025.04.27.sg.curated.neoantigen.predictions/MHC_Class_I/"
                                    "SG.WGS_SG.WGS.UCLA.2025.01.tumor.tsv",),
        verified="2026-09-18"),
    Correction(
        "pvac-header-only-filtered-reports",
        "Seven pVACseq filtered reports contain only a header (byte-identical across runs). "
        "An empty filtered report is not evidence that no epitope passed review.",
        tuple(Change("bucket", {"key": glob(_PVAC + f"*/{group}/*.filtered.tsv")}, expect={"size": size})
              for group, size in (("MHC_Class_I", 2021), ("MHC_Class_II", 1351), ("combined", 2515))),
        evidence=(_BUCKET + _PVAC,), verified="2026-09-18"),
    Correction(
        "pvac-rna-fields-na",
        "RNA depth, RNA VAF, gene and transcript expression are NA in every row of all 21 "
        "published pVACseq reports (131,209 rows) and in their input TSVs; RNA support was "
        "never supplied to pVACseq.",
        (Change("bucket", {"key": glob(_PVAC + "*.tsv")}),),
        evidence=(_BUCKET + _PVAC,), verified="2026-09-18"),
    Correction(
        "pvac-extended-run-no-class-i",
        "Despite its name, the MHCI.extended pVACtools run has no Class-I or combined final "
        "reports; only its Class-II reports exist.",
        (Change("bucket", {"key": glob(_PVAC + "2025.04.27.sg.curated.neoantigen.predictions."
                                              "MHCI.extended/*")}),),
        evidence=(_BUCKET + _PVAC + "2025.04.27.sg.curated.neoantigen.predictions.MHCI.extended/",),
        verified="2026-09-18"),
    _catalogue_resolution("MUC3A-chr7-100953130", "muc3a-grch38-placement"),
    _tempus_relocation("CABLES1-chr18-23135500", ("chr18", 23135500, "G", "dup"),
                       ("chr18", 23135764, "T", "TGGCGGC"), "18:20715728"),
    _tempus_relocation("CCDC40-chr17-80058951", ("chr17", 80058951, "A", "not_reported"),
                       ("chr17", 80090148, "A", "AGAACAACACGGGACGCGCGCAGGCACGTGCAC"), "17:78063947"),
    _tempus_relocation("DCHS2-chr4-154322488", ("chr4", 154322488, "T", "not_reported"),
                       ("chr4", 154323273, "GTTTTTTGCACGACTGCTTCCCAAATGCTGTTTTTCCCTTCAGAGGCATAGGTCTAGCTGCC",
                        "TTTGCACGACTGCTTCCCAAATGCTGTTTTTCCCTACCGAGGCATATGTCTAGCTGCCAA"), "4:155244425"),
    _tempus_relocation("GAPVD1-chr9-125299105", ("chr9", 125299105, "A", "not_reported"),
                       ("chr9", 125301980, "TAGTGC", "ATTGG"), "9:128064259"),
    _tempus_relocation("GOLGA6L2-chr15-23441121", ("chr15", 23441121, "T", "dup"),
                       ("chr15", 23440197, "C", "CTCCCGCATCTTCTCCACCTGCTGCCACATCTTCTGCTCCCGCATTCTCTCCTCCTTC"
                        "TCCCGCAGCCTCTCGTCCTGCTCCCACATCCTCTCCTTCTGGTCCCACATCTTCTGCTCCTGA"), "15:23685344"),
    _allele(
        "MAP2-chr2-209694768", ("chr2", 209694768, "CCTGGGCTACTGTGTGTTCAATA", "C"),
        ("chr2", 209694768, "CCTGGGCTACTGTGTGTTCAATAAGTACACAGT", "CAGGG"),
        "The curated 22-bp deletion (a JLF/mRNA vaccine target) is not the observed allele. Tempus "
        "(freebayes and pindel) and CeGaT describe one net -28 bp complex replacement, "
        "c.2599_2630delinsAGGG; the catalogue's own protein sequence (…DSQLEDRAHCHHLF…) "
        "translates from it, not from the 22-bp deletion. It is written here anchored at the "
        "same position. The vaccine peptide lies downstream in the shared frame and is unaffected.",
        (_ALMY + "pindel.vcf", _ALMY + "freebayes.vcf",
         _BUCKET + "vendor/cegat/P116686_2_S000048/P116686_2_somatic.tsv"),
        extra_source={"genomic_change_on_cdna": "c.2599_2630delinsAGGG"}),
    Correction(
        "map2-split-representations",
        "MAP2-chr2-209694769 (CT>AG, off-site) plus MAP2-chr2-209694772 (28-bp deletion) are "
        "pieces of the same complex event as MAP2-chr2-209694768; do not count them as "
        "independent variants.",
        (Change("vafs", {"variant_id": "MAP2-chr2-209694769"}, expect={"ref": "CT", "alt": "AG"}),
         Change("variant_index", {"id": "MAP2-chr2-209694772"}),
         Change("vafs", {"variant_id": "MAP2-chr2-209694772"},
                expect={"ref": "GGCTACTGTGTGTTCAATAAGTACACAGT", "alt": "G"})),
        evidence=(_ALMY + "pindel.vcf", _BUCKET + "vendor/cegat/P116686_2_S000048/P116686_2_somatic.tsv"),
        verified="2026-09-18"),
    Correction(
        "transcript-DCHS2",
        "DCHS2's RefSeq accession is missing two leading zeros and its version.",
        (Change("source_variants", {"id": "DCHS2-chr4-154322488"}, expect={"refseq_id": "NM_1142552"},
                set={"refseq_id": "NM_001142552.1"}),
         Change("source_variants", {"id": "DCHS2-chr4-154323273"}, absent=True)),
        alternatives=((
            Change("source_variants", {"id": "DCHS2-chr4-154323273"},
                   expect={"refseq_id": "NM_001142552"}, set={"refseq_id": "NM_001142552.1"}),
            Change("source_variants", {"id": "DCHS2-chr4-154322488"}, absent=True),
        ),),
        evidence=(_ALMY + "pindel.vcf", _COORDINATE_FIX), verified="2026-09-21"),
    Correction(
        "fam157a-withdrawn-protein",
        "The 14-residue insertion is annotated on NM_001145248.1, which NCBI has suppressed "
        "(transcript supported, protein not); FAM157A is now only lncRNA NR_146164.1. Tempus "
        "gives a literal allele (GRCh37 3:197880130 G>G+42, GRCh38 chr3:198153259) inside a "
        "low-complexity repeat. The genomic allele is supplied separately from this protein caveat.",
        (Change("variant_index", {"id": "FAM157A-p_W70_Q71ins_14"}),),
        evidence=(_ALMY + "pindel.vcf", "https://www.ncbi.nlm.nih.gov/nuccore/NM_001145248.1",
                  "https://www.ncbi.nlm.nih.gov/gene/728262"),
        verified="2026-09-18"),
    _ush2a_resolution(),
    _catalogue_resolution("FAM157A-p_W70_Q71ins_14", "allele-FAM157A-p_W70_Q71ins_14"),
    _catalogue_resolution("COL3A1-Splice", "allele-COL3A1-Splice"),
    _catalogue_resolution("OTUD4-p_A153del", "otud4-source-unavailable"),
    # --- Timeline and specimen registry -------------------------------------
    Correction(
        "specimen-T1-site",
        "The registry records T1 as a UCSF resection. The homepage summary ('UCLA (Th4 biopsy)'), "
        "the timeline (UCLA biopsies and pathology) and the pathology archive all describe a UCLA "
        "biopsy. T1 actually spans two UCLA biopsies (2024-06-06 T4 and 2024-06-11 T4/5); the "
        "BostonGene libraries (e.g. BG009368) come from the 2024-06-11 procedure, so the date is "
        "left as published.",
        (Change("specimens", {"sample_id": "T1_tumor"},
                expect={"collection_site": "UCSF", "tissue_source": "Recurrent tumor resection"},
                set={"collection_site": "UCLA", "tissue_source": "Recurrent tumor biopsy"}),),
        evidence=(_SITE_REPO + "src/data/samples.json", "https://osteosarc.com/data/events.json",
                  _SITE_REPO + "src/data/pathology-slides.json"),
        verified="2026-09-18"),
    Correction(
        "specimen-T2-date-site",
        "The registry dates T2 2025-01-06 (the bucket delivery folder genomics-bulk/2025.01.06). "
        "The homepage summary, the timeline's time point, biopsy and pathology events, and the "
        "pathology archive (2025-01-28_Biopsy) all give a UCLA biopsy on 2025-01-28.",
        (Change("specimens", {"sample_id": "T2_tumor"},
                expect={"collection_date": "2025-01-06", "collection_site": "UCSF"},
                set={"collection_date": "2025-01-28", "collection_site": "UCLA"}),),
        evidence=(_SITE_REPO + "src/data/samples.json", "https://osteosarc.com/data/events.json",
                  _SITE_REPO + "src/data/pathology-slides.json"),
        verified="2026-09-18"),
    Correction(
        "specimen-T3-site",
        "The registry records T3 as a UCSF biopsy. The homepage summary ('MSKCC (Th4-Th5 excision)'), "
        "the timeline (T4/5 resection surgery, MSKCC pathology) and the pathology archive "
        "(2025-04-17_Resection) describe an MSKCC resection; UCSF performed the single-cell assays.",
        (Change("specimens", {"sample_id": "T3_tumor"},
                expect={"collection_site": "UCSF", "tissue_source": "Tumor biopsy"},
                set={"collection_site": "MSKCC", "tissue_source": "Tumor resection"}),
         Change("specimens", {"sample_id": "T3_tumor_CD45neg"}, expect={"collection_site": "UCSF"},
                set={"collection_site": "MSKCC"})),
        evidence=(_SITE_REPO + "src/data/samples.json", "https://osteosarc.com/data/events.json",
                  _SITE_REPO + "src/data/pathology-slides.json"),
        verified="2026-09-18"),
    Correction(
        "pbmc-capture-dates",
        "These PBMC single-cell specimens are dated by their capture files. The flow-cytometry "
        "workbook (documented as authoritative, and stating that file tokens are not draw dates) "
        "gives draws two to four days earlier: 2025-06-24, 07-22, 08-20 and 09-18. MRD draws fall "
        "on the same flow dates.",
        tuple(Change("specimens", {"sample_id": f"blood_{day}"}, expect={"collection_date": day})
              for day in ("2025-06-26", "2025-07-24", "2025-08-21", "2025-09-22")),
        evidence=("https://osteosarc.com/data/flow/manifest.json",
                  _SITE_REPO + "scripts/flow_data/sample_dates.tsv", "https://osteosarc.com/data/mrd.json"),
        verified="2026-09-18"),
    Correction(
        "events-duplicate-rows",
        "SQ3370 (2023-07-17..2024-01-26) and Trabectedin (2024-05-22) each appear twice as "
        "identical timeline rows; count each course once.",
        (Change("events", {"title": "SQ3370 (Doxorubicin Biogel click)", "date": "2023-07-17"}),
         Change("events", {"title": "Trabectedin", "date": "2024-05-22"})),
        evidence=("https://osteosarc.com/data/events.json",), verified="2026-09-18"),
    Correction(
        "tempus-timepoint",
        "The timeline dates Tempus xT/xE/xR at T0 (2022-12-16), but every public Tempus BAM and VCF "
        "is a TL-24 accession labelled T1 2024-06; the only 2022 Tempus objects are RNA FASTQs. "
        "The catalogue's 'Tempus 2022' detection flag comes from TL-24-ALMY2X4KMV.",
        tuple(Change("events", {"title": title, "date": "2022-12-16"})
              for title in ("Tempus xT", "Tempus xE", "Tempus xR")),
        evidence=(_BUCKET + "vendor/tempus/", _ALMY + "pindel.vcf"), verified="2026-09-18"),
    Correction(
        "apheresis-date",
        "The timeline dates the apheresis 2024-05-14; the ELISPOT records label the same PBMCs "
        "'Apheresis (2024/05/15)'.",
        (Change("events", {"title": "Apheresis", "date": "2024-05-14"}),),
        evidence=("https://osteosarc.com/data/events.json", _SITE_REPO + "src/data/variants.json"),
        verified="2026-09-18"),
    Correction(
        "reyagel-end-date",
        "The timeline sheet gives ReyaGel (started 2026-06-15) the end date '7/14' with no year; "
        "the site's build drops it, so the course appears as a single day. It probably ended "
        "2026-07-14, but the year is not stated, so the end is not supplied here.",
        (Change("events", {"title": "ReyaGel", "date": "2026-06-15"}),
         # Witness: the unreadable published value. A fixed sheet makes this stale.
         Change("events_sheet", {"Title": "ReyaGel", "Start date": "06/15/2026"}, expect={"End date": "7/14"})),
        evidence=(_SITE_REPO + "scripts/timeline/timeline.csv", "https://osteosarc.com/data/events.json"),
        verified="2026-09-18"),
    Correction(
        "natera-alleles-unavailable",
        "COL3A1 (splice) and OTUD4 (p.Ala153del) come from the Natera 2022 report, which is not "
        "public. An independent public Tempus call now resolves COL3A1; OTUD4 remains unavailable.",
        (Change("variant_index", {"id": "COL3A1-Splice"}),
         Change("variant_index", {"id": "OTUD4-p_A153del"})),
        evidence=(_SITE_REPO + "scripts/variants/source_data/SS%20neoantigen%20_%20mutations%20-%20Mutations.tsv",
                  _BUCKET + "vendor/natera/manifest.tsv"),
        verified="2026-09-18"),
)
