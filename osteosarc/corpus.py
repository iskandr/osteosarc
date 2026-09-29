"""Auditable evidence tables for the comprehensive regression corpus.

Published read counts, published rounded VAFs, and balanced fixture template
counts are different measurements. Never pool them, or different processing
products of the same biological sample. A zero is only a measured zero.
"""

from collections import Counter, defaultdict

from .errors import IntegrityError

RNA_ASSAYS = {"RNA", "scRNA", "scRNA_ONT", "PacBio"}
LOW_COVERAGE_CONTROLS = ("CTSE-chr1-206022947", "DBH-chr9-133652279")
NEGATIVE_CONTROLS = ("ACE-chr17-63497361", "BTN3A3-chr6-26448446", *LOW_COVERAGE_CONTROLS)
DISCORDANCE_CONTROLS = ("NR2F2-chr15-96332299",)
AUDIT_SNVS = (*NEGATIVE_CONTROLS, *DISCORDANCE_CONTROLS)
CONTROL_THRESHOLDS = dict(tumor_depth=50, tumor_alt=10, tumor_vaf=0.10,
                          normal_depth=30, normal_max_vaf=0.01, rna_ref=20)
COUNT_FIELDS = ("ref_reads", "alt_reads", "other_reads", "total_reads")


def counts(row):
    """Validate published counts; unknown counts stay None, never zero."""
    result = {}
    for key in COUNT_FIELDS:
        value = row.get(key)
        if value in (None, ""):
            result[key] = None
        else:
            try:
                number = int(value)
            except (ValueError, TypeError) as error:
                raise IntegrityError(f"Invalid {key}: {value!r}") from error
            if str(number) != str(value) or number < 0:
                raise IntegrityError(f"Invalid {key}: {value!r}")
            result[key] = number
    if all(value is not None for value in result.values()):
        if sum(result[k] for k in COUNT_FIELDS[:3]) != result["total_reads"]:
            raise IntegrityError("Published ref + alt + other differs from depth")
    return result


def support_status(values):
    if any(values.get(k) is None for k in COUNT_FIELDS):
        return "not_measured"
    if values["total_reads"] == 0:
        return "no_coverage"
    if values["alt_reads"]:
        return "alt_observed"
    if values["ref_reads"]:
        return "covered_no_alt"
    return "other_only"


def _locus(chrom, pos):
    return (str(chrom).removeprefix("chr").replace("MT", "M"), int(pos))


def vaccine_membership(variants, vaccine_rows):
    """Union independent vaccine assertions, preserving their disagreements.

    Each overlap row must join exactly one literal allele by locus. Duplicate
    rows at the same locus (TECPR1) retain their own assertions. Neither gene
    names nor protein names alone establish allele identity.
    """
    loci = defaultdict(list)
    result = {}
    for v in variants:
        if v["status"] == "ready":
            chrom, pos, _, _ = v["alleles"][0]
            loci[_locus(chrom, pos)].append(v["id"])
        result[v["id"]] = dict(index_count=v.get("vaccine_count"),
                               source_variants=list(v["annotations"].get("source_vaccines", ())),
                               parsed_overlap=list(v.get("vaccines", ())), overlap_rows=[])
    joined = []
    for number, row in enumerate(vaccine_rows):
        matches = loci[_locus(row["chrom"], row["pos"])]
        if len(matches) != 1:
            raise IntegrityError(f"Vaccine row {number} ({row['gene']}) joins {len(matches)} alleles")
        vid = matches[0]
        result[vid]["overlap_rows"].append(number)
        joined.append(dict(row, variant_id=vid, row=number))
    for vid, entry in result.items():
        overlap = {name for i in entry["overlap_rows"] for name, used in vaccine_rows[i]["vaccines"].items()
                   if used}
        entry["overlap_union"] = sorted(overlap)
        entry["vaccines_union"] = sorted(overlap | set(entry["source_variants"]) | set(entry["parsed_overlap"]))
        entry["included"] = bool(entry["vaccines_union"] or (entry["index_count"] or 0) > 0)
        entry["membership_disagreement"] = set(entry["source_variants"]) != overlap
    return result, joined


def count_matrix(variants, rows):
    """Complete variant × reported-BAM matrix, with explicit missing cells."""
    by_key, sources = {}, {}
    source_fields = ("bam_file", "sample_label", "data_source", "pipeline", "assay_type",
                     "timepoint", "tissue", "sample_date")
    for row in rows:
        key = row["variant_id"], row["bam_file"]
        if key in by_key:
            raise IntegrityError(f"Duplicate evidence row: {key}")
        by_key[key] = row
        source = {k: row.get(k, "") for k in source_fields}
        if row["bam_file"] in sources and sources[row["bam_file"]] != source:
            raise IntegrityError(f"Conflicting source metadata: {row['bam_file']}")
        sources[row["bam_file"]] = source
    known = {v["id"] for v in variants}
    if set(vid for vid, _ in by_key) - known:
        raise IntegrityError("Count table contains variants missing from the catalogue")
    matrix = []
    for v in sorted(variants, key=lambda v: v["id"]):
        for bam, source in sorted(sources.items()):
            row = by_key.get((v["id"], bam))
            numbers = counts(row or {})
            status = support_status(numbers)
            if v["status"] != "ready":
                status = "unresolved_allele"
            elif row and status != "not_measured":
                actual = [row["chrom"], int(row["pos"]), row["ref"], row["alt"]]
                if actual != list(v["alleles"][0]):
                    raise IntegrityError(f"Counts are for another allele: {v['id']} / {bam}")
            matrix.append(dict(variant_id=v["id"], **source, **numbers, status=status,
                               reported_vaf=(row or {}).get("vaf") or None,
                               corrections=(row or {}).get("corrections", ""),
                               row_present=row is not None, evidence="published_allele_counts"))
    return matrix, sources


def negative_controls(variants, matrix, ids=NEGATIVE_CONTROLS):
    """Require DNA support and matched-normal evidence, with explicit RNA coverage.

    ALL reported tumor RNA products must have measured zero alternate reads.
    Products are never summed or counted as independent replicates.
    """
    by_variant = defaultdict(list)
    for row in matrix:
        by_variant[row["variant_id"]].append(row)
    by_id = {v["id"]: v for v in variants}
    result = []
    t = CONTROL_THRESHOLDS
    for vid in ids:
        v = by_id[vid]
        if v["status"] != "ready" or any(len(a) != 1 for a in v["alleles"][0][2:]):
            raise IntegrityError(f"Negative control is not a literal SNV: {vid}")
        rows = by_variant[vid]
        rna = [r for r in rows if r["assay_type"] in RNA_ASSAYS and r["tissue"] == "tumor"]
        if not rna or any(r["alt_reads"] != 0 or r["status"] == "not_measured" for r in rna):
            raise IntegrityError(f"Negative control has alternate or unmeasured RNA evidence: {vid}")
        low_coverage = vid in LOW_COVERAGE_CONTROLS
        if low_coverage and max(r["ref_reads"] for r in rna) >= t["rna_ref"]:
            raise IntegrityError(f"Low-coverage control now has covered RNA: {vid}")
        good_dna = [r for r in rows if r["assay_type"] in {"WGS", "WES"} and r["tissue"] == "tumor"
                    and r["status"] == "alt_observed" and r["total_reads"] >= t["tumor_depth"]
                    and r["alt_reads"] >= t["tumor_alt"]
                    and r["alt_reads"] / r["total_reads"] >= t["tumor_vaf"]]
        matches = []
        for dna in good_dna:
            normal = [r for r in rows if r["tissue"] == "blood" and r["assay_type"] == dna["assay_type"]
                      and r["timepoint"] == dna["timepoint"] and r["data_source"] == dna["data_source"]
                      and r["status"] in {"covered_no_alt", "alt_observed"}
                      and r["total_reads"] >= t["normal_depth"]
                      and r["alt_reads"] / r["total_reads"] <= t["normal_max_vaf"]]
            covered = [r for r in rna if r["timepoint"] == dna["timepoint"] and r["ref_reads"] >= t["rna_ref"]]
            if normal and (covered or low_coverage):
                matches.append(dict(tumor_dna=dna["bam_file"], timepoint=dna["timepoint"],
                                    normal_dna=sorted(r["bam_file"] for r in normal),
                                    covered_rna=sorted(r["bam_file"] for r in covered)))
        if not matches:
            raise IntegrityError(f"No matched DNA/normal/covered-RNA comparison for {vid}")
        result.append(dict(variant_id=vid, evidence="published_allele_counts", thresholds=t,
                           control_class="insufficient_rna_coverage" if low_coverage else "covered_rna_no_alt",
                           comparisons=matches, tumor_rna_products_checked=len(rna),
                           interpretation="No alternate RNA reads reported; not proof of absent expression."))
    return result


def snv_read_support(reads, position, ref, alt, *, min_mapping_quality=20, min_base_quality=20):
    """Independent SNV audit of full regional extracts, before balanced sampling.

    Counts templates (RG + QNAME), collapses paired ends, and calls conflicting
    observations 'other'. Missing QUAL is unknown. MAPQ 255 is accepted only with
    NH:i:1 (the STAR convention), never as a numerical quality of 255.
    """
    from .records import read_template
    if len(ref) != 1 or len(alt) != 1:
        raise ValueError("SNV audit requires single-base REF and ALT")
    shown = defaultdict(set)
    skipped = Counter()
    for read in reads:
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_duplicate:
            skipped["flag_filtered"] += 1
            continue
        if read.mapping_quality == 255:
            if not (read.has_tag("NH") and read.get_tag("NH") == 1):
                skipped["unknown_mapping_quality"] += 1
                continue
        elif read.mapping_quality < min_mapping_quality:
            skipped["low_mapping_quality"] += 1
            continue
        offsets = [q for q, r in read.get_aligned_pairs(matches_only=True) if r == position - 1]
        if len(offsets) != 1 or read.query_sequence is None:
            skipped["no_aligned_base"] += 1
            continue
        q = offsets[0]
        if read.query_qualities is None or read.query_qualities[q] == 255:
            skipped["missing_base_quality"] += 1
            continue
        if read.query_qualities[q] < min_base_quality:
            skipped["low_base_quality"] += 1
            continue
        base = read.query_sequence[q].upper()
        shown[read_template(read)].add("alt" if base == alt else "ref" if base == ref else "other")
    result = Counter(next(iter(values)) if len(values) == 1 else "other" for values in shown.values())
    return dict(templates={key: result[key] for key in ("alt", "ref", "other")},
                filtered_records=dict(sorted(skipped.items())), min_mapping_quality=min_mapping_quality,
                min_base_quality=min_base_quality, unit="RG+QNAME", measured_before_sampling=True)


def evidence_catalog(inputs):
    variants = inputs["variants"]
    memberships, vaccine_rows = vaccine_membership(variants, inputs["vaccine_rows"])
    matrix, sources = count_matrix(variants, inputs["counts"])
    controls = negative_controls(variants, matrix)
    historical = inputs.get("historical_rna", [])
    historical_positive = {e["variant_id"] for e in historical if e["observed_before_sampling"].get("alt", 0) > 0}
    if positive := historical_positive & set(NEGATIVE_CONTROLS):
        raise IntegrityError(f"Negative controls have historical RNA support: {sorted(positive)}")
    disagreements = [e for e in historical if e["variant_id"] in DISCORDANCE_CONTROLS
                     and e["observed_before_sampling"].get("alt", 0) > 0]
    if {e["variant_id"] for e in disagreements} != set(DISCORDANCE_CONTROLS):
        raise IntegrityError("The known cross-build RNA disagreement must be retained")
    targets = {}
    source_summaries = []
    for v in variants:
        a = v["annotations"]
        source = a.get("source_record", {})
        allele = list(v["alleles"][0]) if v["status"] == "ready" else None
        targets[v["id"]] = dict(gene=v["gene"], assembly=v["assembly"], allele=allele,
                                 status=v["status"], on_site=v["on_site"],
                                 consequence=a.get("consequence"), protein_change=a.get("protein_change"),
                                 annotated_frameshift="frameshift" in (a.get("consequence") or ""),
                                 length_change_modulo_three=(len(allele[3]) - len(allele[2])) % 3 if allele else None,
                                 vaccine=memberships[v["id"]], pipelines=v["pipelines"],
                                 corrections=a.get("corrections", []),
                                 allele_resolution=a.get("allele_resolution"),
                                 vaccine_peptides=source.get("vaccine_peptides", []),
                                 peptides=source.get("peptides", {}),
                                 validation=source.get("validation", {}))
        for row in source.get("vafs", []):
            source_summaries.append(dict(variant_id=v["id"], **row, evidence="published_depth_and_rounded_vaf",
                                         interpretation="Not raw allele counts; do not infer integer counts or pool products."))
    for control in controls:
        control["scope"] = "The tumor RNA products in the published allele-count matrix."
        control["nonzero_rna_summary_exceptions"] = [
            row for row in source_summaries if row["variant_id"] == control["variant_id"]
            and row["assay"] not in {"WGS", "WES"} and (row.get("vaf") or 0) > 0]
        if any(row["tissue"] == "tumor" for row in control["nonzero_rna_summary_exceptions"]):
            raise IntegrityError(f"Negative control has a nonzero tumor RNA summary: {control['variant_id']}")
    overlap_groups = inputs.get("evidence_overlap_registry", {}).get("groups", [])
    count_sources = {k: dict(v, **inputs["count_sources"].get(k, {})) for k, v in sources.items()}
    for source in count_sources.values():
        urls = {f["url"] for f in source.get("files", [])}
        source["evidence_overlap_groups"] = [g["id"] for g in overlap_groups
                                             if urls & {s["url"] for s in g["sources"]}]
        if source["evidence_overlap_groups"]:
            source["sample_identity_caution"] = "Published labels retained; biological specimen/timepoint unresolved."
    return dict(schema_version=1, snapshot=inputs["snapshot"], provenance=inputs["provenance"],
                suitable_for_abundance=False, variants=targets, vaccine_rows=vaccine_rows,
                count_sources=count_sources, evidence_overlaps=overlap_groups,
                allele_support=matrix, source_reported_support=source_summaries,
                negative_controls=controls, samples=inputs["samples"],
                historical_rna=historical, known_rna_disagreements=disagreements,
                structural_candidates=inputs["structural_candidates"],
                summary=dict(variants=len(targets), vaccine_variants=sum(m["included"] for m in memberships.values()),
                             vaccine_rows=len(vaccine_rows), count_sources=len(sources), support_cells=len(matrix),
                             rna_count_sources=sum(s["assay_type"] in RNA_ASSAYS for s in sources.values()),
                             annotated_frameshifts=sum(t["annotated_frameshift"] for t in targets.values()),
                             structural_candidates=len(inputs["structural_candidates"]["targets"]),
                             samples=len(inputs["samples"]), source_summary_rows=len(source_summaries),
                             source_summary_products=len({r["sample"] for r in source_summaries}),
                             statuses=dict(Counter(r["status"] for r in matrix)),
                             negative_controls=len(controls)))
