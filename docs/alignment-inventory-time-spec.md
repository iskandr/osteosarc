# Enforce snapshot modification times during alignment access

Issue: https://github.com/iskandr/osteosarc/issues/112

1. Share the whole-file downloader's numeric inventory timestamp comparison
   with remote alignment header inspection and indexed extraction.
   Compare UTC epochs using the existing one-second tolerance. HTTP
   Last-Modified has second precision; see RFC 9110 sections
   [5.6.7](https://www.rfc-editor.org/rfc/rfc9110.html#section-5.6.7) and
   [8.8.2](https://www.rfc-editor.org/rfc/rfc9110.html#section-8.8.2).
2. Reject older/newer remote dates before reading a header or downloading an
   index. Validate recorded identities when reopening cached headers and read
   derivatives, including old receipts, without network access or rewriting
   preserved provenance. Existing request identities and byte pins remain valid.
3. Retain the existing available-evidence policy: a nonnumeric/missing inventory
   time or missing server date cannot establish a date match. New alignment
   receipts explicitly record `inventory_modification` as `matched`,
   `inventory_timestamp_unavailable`, or `last_modified_unavailable`.
   Invalid supplied dates or nonfinite numeric inventory times raise
   IntegrityError. Local files continue to use their local identity checks.
4. Cover actual indexed BAM access via a local HTTP server, equal/fractional
   times, mismatch in both directions, missing evidence, malformed dates,
   extraction/recovery entry points and offline legacy cache validation.
   Confirm mismatches publish no completed acquisition. Run lint, full tests,
   package build and Python/consumer CI. Bump to 0.15.5, merge and publish.

No biological consequence or read-filter semantics change. This is the source
identity prerequisite for resuming the remaining selected Sid SV acquisitions.
