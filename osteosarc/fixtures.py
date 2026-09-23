"""Versioned fixture panels and deterministic source-scoped record selection.

Evidence assignments are pinned inputs from scientific producers. This module
executes their policy; it does not infer allele, junction or protein support.
"""

from __future__ import annotations

import copy
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from .cache import digest, stable_id
from .errors import IntegrityError, SchemaError
from .models import Region
from .reads import ReadSubset, assembly_from_header, normalize_assembly, resolve_regions
from .records import RECORD_ENCODING, read_records, sam_digest


def _count(value, label):
    if type(value) is not int or value < 0:
        raise SchemaError(f"{label} must be a nonnegative integer")


def validate_recipe(recipe):
    """Validate and copy a v1 recipe; coordinates and evidence must be explicit."""
    if not isinstance(recipe, dict):
        raise SchemaError("A fixture recipe must be an object")
    recipe = copy.deepcopy(recipe)
    if (type(recipe.get("schema_version")) is not int or recipe["schema_version"] != 1
            or not isinstance(recipe.get("id"), str) or not recipe["id"]):
        raise SchemaError("Expected a named fixture recipe with schema_version=1")
    for key in ("targets", "sources", "members"):
        if not isinstance(recipe.get(key), dict):
            raise SchemaError(f"Recipe requires a {key} mapping")
        if any(not isinstance(name, str) or not name or not isinstance(value, dict)
               for name, value in recipe[key].items()):
            raise SchemaError(f"Recipe {key} must map nonempty names to objects")
    for name, target in recipe["targets"].items():
        kind = target.get("kind")
        if kind == "unresolved":
            if not target.get("reason"):
                raise SchemaError(f"Unresolved target {name} requires a reason")
            continue
        if not target.get("assembly") or not target.get("reference"):
            raise SchemaError(f"Target {name} requires pinned assembly/reference identity")
        if kind == "small_variant":
            if target.get("coordinates") != "one-based" or not target.get("contig"):
                raise SchemaError("Small variants require one-based coordinates and a contig")
            _count(target.get("position"), "position")
            if target["position"] == 0 or any(not target.get(k) for k in ("ref", "alt")):
                raise SchemaError("Small variants require position >= 1 and explicit alleles")
        elif kind == "sv":
            if (target.get("coordinates") != "zero-based-interbase"
                    or not isinstance(target.get("breakends"), list)
                    or len(target["breakends"]) < 2
                    or any(not isinstance(end, dict) for end in target["breakends"])):
                raise SchemaError("SVs require at least two explicit interbase breakends")
            for end in target["breakends"]:
                _count(end.get("position"), "breakend position")
                if not end.get("contig") or end.get("orientation") not in ("+", "-", None):
                    raise SchemaError("Invalid oriented breakend")
        else:
            raise SchemaError(f"Unknown target kind: {kind}")
    for name, source in recipe["sources"].items():
        if (not isinstance(source.get("identity"), dict) or not source["identity"]
                or not source.get("assembly")):
            raise SchemaError(f"Source {name} requires identity and assembly")
        for field in ("sample", "library", "product"):
            if field not in source:
                raise SchemaError(f"Source {name} must declare {field} (null if unresolved)")
    for name, member in recipe["members"].items():
        if (not isinstance(member.get("target"), str) or not isinstance(member.get("source"), str)
                or member["target"] not in recipe["targets"] or member["source"] not in recipe["sources"]):
            raise SchemaError(f"Member {name} refers to an unknown target/source")
        policy = member.get("policy", {})
        if (not isinstance(policy, dict) or type(policy.get("version")) is not int
                or policy["version"] != 1 or policy.get("kind") not in (
                    "regional", "exact", "witnesses", "stratified", "empty", "omitted")):
            raise SchemaError(f"Unknown selection policy for {name}")
        if policy["kind"] == "omitted" and not policy.get("reason"):
            raise SchemaError("An omission must state its reason")
        if "cap" in policy:
            _count(policy["cap"], "cap")
        if not isinstance(policy.get("strata", {}), dict):
            raise SchemaError("Stratum caps must be a mapping")
        for n in policy.get("strata", {}).values():
            _count(n, "stratum cap")
        if policy.get("duplicate_policy", "preserve") not in ("preserve", "identical-record-once"):
            raise SchemaError("Unknown duplicate policy")
        if policy["kind"] == "exact":
            if policy.get("encoding", RECORD_ENCODING) not in (RECORD_ENCODING, "sam-text-v1"):
                raise SchemaError("Unsupported record identity encoding")
            if not isinstance(policy.get("records"), dict):
                raise SchemaError("Exact selection requires an explicit records mapping (empty if intentional)")
            for checksum, n in policy["records"].items():
                if (not isinstance(checksum, str) or len(checksum) != 64
                        or any(c not in "0123456789abcdef" for c in checksum)):
                    raise SchemaError("Invalid record checksum")
                _count(n, "record multiplicity")
                if not n:
                    raise SchemaError("Record multiplicity must be positive")
        assignments = policy.get("assignments", [])
        if not isinstance(assignments, list) or any(not isinstance(a, dict) for a in assignments):
            raise SchemaError("Evidence assignments must be a list of objects")
        for assignment in assignments:
            if (not assignment.get("reason") or not isinstance(assignment.get("producer"), dict)
                    or not assignment["producer"].get("name") or not assignment["producer"].get("version")):
                raise SchemaError("Witness/stratum assignments require reason and producer name/version")
            selector = assignment.get("selector", {})
            if (not isinstance(selector, dict) or not isinstance(selector.get("qname"), str)
                    or not selector["qname"] or "rg" not in selector):
                raise SchemaError("Witness selectors require qname and explicit rg (null if absent)")
            if selector["rg"] is not None and not isinstance(selector["rg"], str):
                raise SchemaError("Witness read groups must be strings or null")
            if "required" in assignment and type(assignment["required"]) is not bool:
                raise SchemaError("Witness required must be a boolean")
            if selector.get("segment") not in (None, 0, 64, 128, 192):
                raise SchemaError("Segment is the original FLAG & 0xc0")
        for field in ("regions", "context_regions"):
            if not isinstance(member.get(field, []), list):
                raise SchemaError(f"{field} must be a list of region objects")
            for region in member.get(field, []):
                try:
                    Region(**region)
                except (TypeError, ValueError) as error:
                    raise SchemaError(f"Invalid {field} in member {name}: {error}") from error
        if policy["kind"] == "regional" and not member.get("regions"):
            raise SchemaError("Regional selection requires explicit bounded regions")
    return recipe


def load_panel(name):
    """Load a shipped named target panel without network or downstream imports."""
    path = Path(__file__).with_name("data") / "panels.json"
    panels = json.loads(path.read_text())
    if name not in panels:
        raise KeyError(f"Unknown panel {name!r}; available: {', '.join(sorted(panels))}")
    return copy.deepcopy(panels[name])


@dataclass
class FixtureSelection:
    recipe: dict
    members: dict
    headers: dict
    records: dict
    receipts: dict

    @property
    def manifest(self):
        return dict(schema_version=1, recipe_sha256=stable_id(self.recipe),
                    record_encoding=RECORD_ENCODING, members=self.members,
                    sources=self.recipe["sources"], acquisition=self.receipts,
                    suitable_for_abundance=False)


def _matches(record, selector):
    return (record.template == (selector["rg"], selector["qname"])
            and (selector.get("segment") is None or record.segment == selector["segment"]))


def _overlaps(record, regions):
    r = record.read
    return any(r.reference_name == region.contig and r.reference_start < region.end
               and (r.reference_end or r.reference_start + 1) > region.start for region in regions)


def select_fixture_records(records, policy, *, regions=(), context_regions=()):
    """Execute one policy over FixtureRecords; return counts, reasons and status.

    Scope is one declared source. Duplicate occurrences survive; shared members
    refer to the same source records. Required witnesses bypass sampling caps.
    A cap limits optional templates, never records belonging to a kept template.
    """
    # Public callers may pass read_records() directly. Selection revisits the
    # input for witnesses, strata and context, so consume an iterator only once.
    records = tuple(records)
    regions, context_regions = tuple(regions), tuple(context_regions)
    available = Counter(r.digest for r in records)
    reasons = defaultdict(set)
    kind = policy["kind"]
    counts = Counter()
    truncated = False
    if kind in ("omitted", "empty"):
        return counts, {}, "omitted" if kind == "omitted" else "empty"
    if kind == "exact":
        required = Counter(policy.get("records", {}))
        if policy.get("encoding", RECORD_ENCODING) == "sam-text-v1":
            by_sam = defaultdict(Counter)
            for record in records:
                by_sam[sam_digest(record.read.to_string())][record.digest] += 1
            if required - Counter({k: sum(v.values()) for k, v in by_sam.items()}):
                raise IntegrityError("Missing pinned SAM-text records or duplicate occurrences")
            for key, n in required.items():
                if len(by_sam[key]) != 1:
                    raise IntegrityError("Ambiguous SAM identity maps to distinct binary records")
                counts[next(iter(by_sam[key]))] = n
        else:
            if missing := required - available:
                raise IntegrityError(f"Missing pinned records or duplicate occurrences: {dict(missing)}")
            counts = required
        for key in counts:
            reasons[key].add(policy.get("reason", "pinned historical record"))
    else:
        groups = defaultdict(list)
        for r in records:
            groups[r.template].append(r)
        selected = defaultdict(set)
        strata = defaultdict(set)
        for a in policy.get("assignments", []):
            matches = [r for r in records if _matches(r, a["selector"])]
            if a.get("required", True) and not matches:
                raise IntegrityError(f"Missing pinned witness: {a['selector']}")
            for r in matches:
                if a.get("required", True) or kind == "witnesses":
                    selected[r.digest].add(a["reason"])
                elif kind == "stratified":
                    strata[a.get("stratum", a["reason"])].add(r.template)
        if kind == "regional":
            candidates = {r.template for r in records if _overlaps(r, regions)}
            strata["regional background"] = candidates
        for stratum, templates in sorted(strata.items()):
            ordered = sorted(templates, key=lambda t: (stable_id([str(policy.get("seed", "0")), t]), str(t)))
            cap = policy.get("strata", {}).get(stratum, policy.get("cap", len(ordered)))
            truncated |= len(ordered) > cap
            for template in ordered[:cap]:
                for r in groups[template]:
                    if kind != "regional" or _overlaps(r, regions):
                        selected[r.digest].add(stratum)
        for key, why in selected.items():
            counts[key] = available[key]
            reasons[key].update(why)
    for r in records:
        if _overlaps(r, context_regions):
            counts[r.digest] = available[r.digest]
            reasons[r.digest].add("assembly context")
    if policy.get("duplicate_policy", "preserve") == "identical-record-once":
        counts = Counter(dict.fromkeys(counts, 1))
    elif policy.get("duplicate_policy", "preserve") != "preserve":
        raise SchemaError("Unknown duplicate policy")
    return counts, {k: sorted(reasons[k]) for k in sorted(counts)}, "truncated" if truncated else "selected" if counts else "empty"


def select_fixtures(recipe, sources):
    """Select from explicit local archived BAMs/ReadSubsets, without acquisition.

    ``sources`` maps recipe source IDs to paths or verified ReadSubsets. Remote
    inputs must first be acquired through indexed extraction. Local archives may
    pin a SHA256; source metadata and producer assignments enter the recipe hash.
    """
    import pysam
    recipe = validate_recipe(recipe)
    records, headers, receipts = {}, {}, {}
    members = {}
    for name, member in sorted(recipe["members"].items()):
        sid = member["source"]
        policy = member["policy"]
        target = recipe["targets"][member["target"]]
        if target["kind"] == "unresolved" or policy["kind"] == "omitted":
            members[name] = dict(source=sid, target=member["target"], records={}, reasons={}, record_count=0,
                                 status="unresolved" if target["kind"] == "unresolved" else "omitted",
                                 explanation=target.get("reason", policy.get("reason")))
            continue
        if sid not in records:
            value = sources[sid]
            path = Path(value.path if isinstance(value, ReadSubset) else value)
            expected = recipe["sources"][sid].get("archive_sha256")
            before = digest(path)
            if isinstance(value, ReadSubset):
                request = value.receipt.get("request", {})
                request = request.get("seed", request)
                acquired_source = request.get("source", "")
                declared = recipe["sources"][sid]["identity"]
                if acquired_source.startswith(("https://", "http://")):
                    if declared.get("url") and declared["url"] != acquired_source:
                        raise IntegrityError(f"Acquired source identity differs: {sid}")
                elif expected and request.get("source_sha256") != expected:
                    raise IntegrityError(f"Acquired archive source checksum differs: {sid}")
            elif expected and before != expected:
                raise IntegrityError(f"Archive checksum mismatch: {sid}")
            if isinstance(value, ReadSubset) and value.receipt.get("files", {}).get("reads.bam") != before:
                raise IntegrityError(f"Acquired BAM checksum mismatch: {sid}")
            with pysam.AlignmentFile(path) as bam:
                header = bam.header.to_dict()
            assembly = assembly_from_header(header)
            if assembly != normalize_assembly(recipe["sources"][sid]["assembly"]):
                raise IntegrityError(f"Source assembly cannot be established or differs: {sid}")
            headers[sid] = header
            records[sid] = list(read_records(path))
            if digest(path) != before:
                raise IntegrityError(f"Source changed during selection: {sid}")
            receipts[sid] = value.receipt if isinstance(value, ReadSubset) else dict(archive_sha256=before)
        if normalize_assembly(target["assembly"]) != normalize_assembly(recipe["sources"][sid]["assembly"]):
            raise IntegrityError(f"Target/source assembly mismatch for {name}")
        regions = resolve_regions([Region(**r) for r in member["regions"]], headers[sid]) if member.get("regions") else ()
        context = resolve_regions([Region(**r) for r in member["context_regions"]], headers[sid]) if member.get("context_regions") else ()
        counts, reasons, status = select_fixture_records(records[sid], policy, regions=regions, context_regions=context)
        acquisition = receipts[sid]
        if member.get("retain_partners") and acquisition.get("scope") == "bounded_mate_SA_context":
            templates = {r.template for r in records[sid] if r.digest in counts}
            available = Counter(r.digest for r in records[sid])
            for r in records[sid]:
                why = set(acquisition["reasons"].get(r.digest, [])) & {"paired mate", "SA-linked partner"}
                if r.template in templates and why:
                    counts[r.digest] = available[r.digest]
                    reasons[r.digest] = sorted(set(reasons.get(r.digest, [])) | why)
        members[name] = dict(source=sid, target=member["target"], records=dict(sorted(counts.items())),
                             acquisition_status=acquisition.get("status", "available-input"),
                             reasons=reasons, status=status, record_count=sum(counts.values()))
    return FixtureSelection(recipe, members, headers, records, receipts)
