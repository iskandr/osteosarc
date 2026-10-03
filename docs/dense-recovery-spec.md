# Dense seed context (#109)

Independent indexed counts found 1,044,612 reads in the failed 13-target 10x
batch and 1.35–5.67 million reads in individual dense targets. Retry the original
failures separately; neither a regional record cap nor a subprocess failure is
evidence against an expressed junction.

Stream seed records, retaining alignment objects only for templates with outgoing
mate/SA pointers or diagnostics. Include every seed record belonging to a
pointer-bearing name before resolving leads, and query those names in every round
so later leads can reuse earlier windows. Preserve all original seed records and
their multiplicities in a coordinate-ordered output stream; append only matched
partners. Verify the complete output record multiset. Keep per-record provenance,
unresolved pointers, caps and timeouts explicit. Version the changed recovery
request because the partner query-name set changes.

Regression checks cover ordinary context retention, seed partners without their
own pointers, recursive/cyclic pointers, duplicates, read groups and missing
QNAME/SEQ diagnostics. Measure memory and successful acquisition on the dense
pilot inputs before changing Isovar's configurable resource budgets. The pipe
cleanup failure remains a distinct issue unless reproduced and fixed.
