# Bounded BAM pipe lifetime (#106)

The Sid ONT retry reproduced the truncated BGZF-header and `Illegal seek`
destructor diagnostic. A timeout can kill samtools while its last header block
is still buffered during a sparse partner query. Pysam then attempts to clean up
a partially constructed AlignmentFile on a nonseekable pipe. Preserve the real
timeout, overflow, or producer error without this secondary exception.

Read the forward-only BAM stream with Python's gzip reader and copy its original
header and record blocks into a local BGZF file. Count complete record blocks
before the cap; never reconstruct their SEQ, QUAL, CIGAR, tags or floating-point
values. Use the BAM framing from SAMv1 section 4.2, checking lengths and truncated
blocks. Keep the subprocess deadline, kill/reap cleanup and partial-cache guard.

Regress successful large headers, a timeout midway through a BGZF header, a
truncated record, record overflow, nonzero producer exit with stderr, and a
subsequent successful call. Verify exact record multiplicities against the source.
This fixes pipe cleanup; it does not turn a timed-out remote query into negative
RNA evidence or establish that an earlier subprocess failure had the same cause.
