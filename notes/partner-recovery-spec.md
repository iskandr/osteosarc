# Bounded partner-query batches (#104)

Large regional RNA inputs can nominate thousands of mate/SA locations. A single
slow combined query currently times out before any of its available partners are
retained, even though the seed BAM is verified.

Add optional limits for regions per partner query, total partner queries and
per-partner-query seconds. Preserve default requests/cache identities and all
source/RG/segment/placement matching. When `on_timeout="incomplete"`, retain
completed batches, label only the failed batch's leads as timed out, and continue
with other batches. Every attempted query is charged to the query budget. Failed
queries never become visited regions or evidence of absence. Subsequent calls
reuse verified completed extractions and retry missing work.

Test late partners after an early timeout, exact record multiplicity, source
drift, exhausted budgets, input validation and offline cache reuse. Keep seed
acquisition semantics unchanged. This change supports Isovar's catalogue-wide
SV RNA audit without synthesizing observations from SAM pointers.
