# Evidence for catalogue issue #5

Retrieved 2026-09-20. `provenance.json` records the original URLs and SHA256s,
plus separate hashes for excerpts. VCF and FASTA coordinates are one-based;
UCSC mapping blocks in `osteosarc/data/allele_resolutions.json` are zero-based,
half-open. Source annotations and the public Natera catalogue rows retain
their original content; the JSON excerpts select fields and records only.

* **FAM157A:** Tempus Pindel GRCh37 3:197880130 G>G+42; GRCh38
  chr3:198153259. The complete 1,001-base RefSeq windows are identical.
  Equivalent anchored positions span 198153259–198153261; the first is used.
  NCBI's Gene record separately documents suppression of the protein model.
* **COL3A1:** Tempus Pindel GRCh37 2:189875615, 738-base REF>C;
  GRCh38 chr2:189010889. NM_000090.3:c.4254+1_4255-1del matches the
  catalogue. Both 1,738-base windows, including the entire REF, are identical.
  Equivalent anchored positions are 189010889 and 189010890; the first is used.
* **MUC3A:** The 102-base Tempus duplication's anchor falls in the source
  chain gap [100550653, 100550873). The corresponding GRCh38 gap is
  [100953013, 100958744), containing the catalogue position. The inserted
  unit does not occur verbatim in the checked GRCh38 window. This does not
  establish a GRCh38 allele. It replaces the older, insufficient point-mapping
  claim in the correction notes.
* **OTUD4:** The public Natera-derived catalogue reports only p.Ala153del;
  no transcript or genomic allele. The checked Tempus somatic/germline and
  CeGaT calls supplied no match. Resolution needs the original clinical call.
* **USH2A:** Tempus Freebayes GRCh37 1:215824094 C>A maps to the existing
  chr1:215650752 C>A entry; its RefSeq windows are identical. A transposed
  coordinate and matching protein label suggest a relationship to
  chr1:215560752 but cannot prove the original Natera call's identity.

`python -m pytest tests/test_allele_resolution.py` checks the evidence offline.
`python scripts/check_allele_sources.py --cache CACHE` verifies the whole
source checksums and chain mappings; add `--offline` to reuse cached originals.
The UCSC chain is not redistributed here. The script verifies the chain bytes
against the recorded checksum before interpreting its blocks.
