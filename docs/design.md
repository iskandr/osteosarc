# Shared osteosarc data API

The package owns data acquisition and source interpretation. Varcode owns
variant effects, Isovar owns RNA interpretation, Topiary owns annotation and
ranking, and Vaxrank owns vaccine design/reporting. Importing osteosarc does
not download data or import those projects.

## Contracts

* A `Dataset` is a named, immutable snapshot of source metadata. Creating a
  new snapshot is explicit. Opening an existing snapshot is offline.
* An `Asset` is one published object, not a biological replicate. Its stable
  ID derives from its complete URL. BAM processing products remain separate.
  Sample claims retain their source and conflicts; unknown is not inferred
  to mean Illumina, tumor, T0, GRCh38, or negative evidence.
* A shared `Cache` stores immutable downloaded bytes and SHA256 receipts.
  Reuse verifies bytes. Refresh is explicit; old snapshot references survive
  refresh. Downloads publish atomically under a per-URL lock.
* The full bucket listing remains queryable, including unclassified files.
  Curated BAM metadata enriches it, rather than defining the entire dataset.
  The website listing has its own date; it is not a live S3 inventory. An
  explicit paginated S3 listing can discover newer objects.
* Parsing is separate from downloading. Raw rows and original field names
  survive TSV/CSV/JSON parsing. VCF parsing delegates to pysam or varcode.
  Missing measurements remain missing, independently of zero measurements.
* Variant identities include assembly and literal alleles. Website entries
  without unique literal alleles remain visible with a status. Gene symbols
  and protein labels are never substitutes for exact alleles. Source-reported
  coordinates are not silently corrected or lifted over.
* `Region` uses zero-based, half-open coordinates and requires assembly.
  SAMtools conversion is one-based inclusive. Read extraction inspects BAM
  headers, rejects assembly conflicts/ambiguous contigs/missing indexes, and
  uses the indexed union. It retains original record multiplicity and tags.
  Defaults impose no allele, quality, duplicate, or template-count selection.
* Derived BAMs are cached by the complete request and source evidence, with
  checksums, tool version, filters and counts. An empty request is an error;
  an empty *result* is valid. No implicit full-BAM fallback.
* Fixture downsampling is a separate explicit operation. It cannot silently
  become the read source for VAF calculations.

## Deliverables

Python package, CLI, offline regression tests, reproducible live smoke recipe,
and migration guidance for the four existing repositories. Downstream
repositories are reviewed read-only during this initial extraction.
