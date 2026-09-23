# Frozen additional SV research witnesses

These four original BAM subsets and their original receipts were handed off in
Osteosarc issue #15 on 2026-09-23. Public source URLs, source-call records, original
query intervals, direction and missing-quality evidence remain in the saved
JSON inputs. The data license is CC0-1.0. These are research candidates, not
independently validated biological or clinical interpretations.

The shared recipe is `osteosarc/data/additional_sv_recipe.json`. It preserves
all 251 archived source records in shared context members, plus exact per-event,
per-source witness membership. Reuse across members is not independent evidence.
Candidate reverse-direction matches remain labeled reverse; missing PacBio
qualities remain unknown. Original-query intervals are never transformed into
new witnesses. Source/RG/QNAME/segment selectors are pinned, not gene names.

Regenerate entirely offline from the repository root:

```sh
python -m osteosarc --offline fixtures generate osteosarc/data/additional_sv_recipe.json /tmp/additional-sv-bundle \
  --source T1-ONT-tagged=tests/data/additional_svs/T1-ONT-tagged.bam \
  --source T2-ONT-tagged=tests/data/additional_svs/T2-ONT-tagged.bam \
  --source T1-short=tests/data/additional_svs/T1-short.bam \
  --source T1-PacBio=tests/data/additional_svs/T1-PacBio.bam
python -m osteosarc fixtures verify /tmp/additional-sv-bundle
```

[VCF §5.4](https://samtools.github.io/hts-specs/VCFv4.3.pdf) breakend anchor coordinates are retained next to interbase boundaries. A
boundary after a retained-left anchor is POS; before a retained-right anchor it
is POS-1. Traversal orientation/inserted sequence stay unspecified where the
handoff does not independently establish them; original ALT strings remain
available. Consumers own reconstruction and interpretation assertions.
