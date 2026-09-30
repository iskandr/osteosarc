# Regional Tempus overlap evidence

These two unmodified regional BAMs reproduce Osteosarc #100. Each contains reads
at the same 44 vaccine loci with 2,000-base flanks plus fetched mates. Original
acquisition receipts, remote identities, request regions and source URLs are in
`receipts.json`; the files are pinned in `osteosarc/data/evidence-overlaps.json`.
The original public source reads are CC0.

Run `python -m pytest tests/test_provenance.py` to verify the inputs and reproduce
the comparison without network access. It counts record occurrences and includes
base qualities, while excluding optional tags such as differing RG labels.
There are 40,684 identical core records, 12 different records per BAM and 18,971
shared read names. The files contain 40,696 records each.

The result applies to these regions, not whole libraries. Both original RG/SM
labels remain intact. The full published FASTQ catalogue is also pinned here: it
explicitly describes KCVBE1UI1P as a duplicate delivery of ALMY2X4KMV RNA. The
registry records matching lengths and multipart ETags for the corresponding mate
objects; these are server metadata, not whole-file SHA-256 hashes. The biological
specimen/timepoint remains disputed by the existing `tempus-file-labels` curation.
Specimen documentation is still needed before relabeling it.
