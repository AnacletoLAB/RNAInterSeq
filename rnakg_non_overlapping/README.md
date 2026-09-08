# RNA-KG records that do **not** overlap RNAInterSeq

ncRNA–ncRNA interaction pairs present in RNA-KG and absent from RNAInterSeq,
intended as a leakage-free evaluation set for models trained on RNAInterSeq.

## Why this exists

**RNAInterSeq contains no train/test division.** There is no split, fold or
partition column in the released dataframe, and no splitting code anywhere in
the repository. Any model trained on RNAInterSeq therefore needs held-out data
from elsewhere, and RNA-KG is the natural candidate — but only after every
interaction RNAInterSeq already contains has been removed from it, otherwise
the evaluation leaks.

## Headline numbers

Source: `data/rna-kg_ncRNA_new.pkl` — 112,120 ncRNA–ncRNA records.

| | records | share |
|---|---:|---:|
| overlap RNAInterSeq (removed) | 90,716 | 80.91 % |
| **non-overlapping (kept)** | **21,404** | **19.09 %** |

At the entity level, of the 14,700 distinct sequences in the RNA-KG ncRNA set,
**2,885** never occur anywhere in RNAInterSeq.

## How the non-overlapping set is constructed

Run `python build_non_overlapping.py`. Steps:

1. **Reference set.** Load the released ncRNA–ncRNA dataset
   (`data/interactions_with_sequences_RNA_RNA_ncrna_ncrna_19-05-2026.pkl`,
   119,365 interactions). RNA-KG's ncRNA export is ncRNA–ncRNA only, so this is
   the matching scope; comparing instead against the full RNAInterSeq release
   gives exactly the same result (90,716 overlapping / 21,404 kept either way).

2. **Expand to transcript level.** RNAInterSeq stores interactions at the
   RNA-*entity* level and keeps every transcript isoform recovered for an
   entity; RNA-KG stores them at the transcript level. The two resources share
   no identifier, so **the join key is the nucleotide sequence itself** — the
   same choice made in the RNA-KG comparison cells of `mRNAInter_stats.ipynb`
   and `ncRNAInter_stats.ipynb`. Each interaction is expanded into all
   combinations of the transcripts of its two partners:
   692,913 ordered transcript pairs over 34,761 distinct sequences.

3. **Normalise both sides** — upper-case, `T`→`U`, strip. RNA-KG sequences are
   already clean RNA alphabet; RNAInterSeq mRNA sequences come raw from RefSeq
   in DNA alphabet, so this is what makes the two comparable.

4. **Canonicalise pair orientation.** `(A,B)` and `(B,A)` denote the same
   interaction, and the two resources do not agree on partner order. Pairs are
   sorted before comparison, giving 676,710 distinct canonical RNAInterSeq
   pairs.

5. **Subtract.** An RNA-KG record is *overlapping* if its canonical normalised
   sequence pair occurs in that set; everything else is kept.

6. **Entity level.** An RNA-KG sequence is non-overlapping if it appears
   nowhere in RNAInterSeq, as either partner, in any transcript.

Matching is **exact, not hashed**: every distinct RNAInterSeq sequence gets an
integer id and a canonical pair is encoded as one `int64`, so there is no
collision risk. An RNA-KG sequence absent from the RNAInterSeq vocabulary
cannot participate in any RNAInterSeq pair, so its record is non-overlapping
by construction.

## Two deliberate differences from the notebook comparison cells

Both are **conservative** — they can only find *more* overlap, i.e. remove more
records, never add unsafe ones.

| | notebook cells | here |
|---|---|---|
| sequence normalisation | mRNA only; the ncRNA cell compares raw strings | both subsets |
| pair orientation | ordered `(x,y)` | canonical (unordered) |
| comparison scope | matching subset only | matching subset (identical to full-release scope here) |

Measured effect on the ncRNA set:

| source | ordered (notebook method) | canonical (used here) | extra caught |
|---|---:|---:|---:|
| `rna-kg_ncRNA.pkl` (older) | 68,536 | 68,601 | 65 |
| `rna-kg_ncRNA_new.pkl` | 90,632 | 90,716 | 84 |

The ordered figure on the older file, **68,536**, reproduces the notebook
output exactly, which is what validates this pipeline against the existing
analysis. Restricting the comparison to RNAInterSeq's ncRNA–ncRNA subset
instead of the full dataset changes nothing here (90,716 either way): no RNA-KG
ncRNA pair occurs only inside the mRNA-associated subset.

## Which RNA-KG file

`data/` holds two ncRNA extractions and they are **not** nested:

* `rna-kg_ncRNA.pkl` — 84,345 pairs, the one the notebooks use;
* `rna-kg_ncRNA_new.pkl` — 112,120 pairs, newer; **used here**.

83,562 pairs are shared, 28,558 exist only in the newer file, and **783 exist
only in the older one** (mostly `snoRNA-snoRNA`). Those 783 are not silently
dropped: the **117** of them that are also absent from RNAInterSeq (103
`snoRNA-snoRNA`, 9 `miRNA-snoRNA`, 5 others) are written to
`ncrna_pairs_non_overlapping_oldfile_only.pkl`. Add them if you prefer to treat
the two extractions as a union, giving 21,521 pairs instead of 21,404.

## Files

| file | content |
|---|---|
| `data/ncrna_pairs_non_overlapping.pkl` | **main deliverable** — 21,404 RNA-KG ncRNA–ncRNA pairs absent from RNAInterSeq |
| `data/ncrna_entities_non_overlapping.pkl` | 2,885 RNA-KG ncRNA sequences absent from RNAInterSeq |
| `data/ncrna_pairs_non_overlapping_oldfile_only.pkl` | 117 non-overlapping pairs found only in the older extraction |
| `data/overlap_summary.json` | all counts, machine-readable |
| `build_non_overlapping.py` | the script that produces the above |

`ncrna_pairs_non_overlapping.pkl` keeps RNA-KG's own columns:
`RNA_sequence_x`, `RNA_sequence_y` (raw, *not* normalised — normalisation is
used for matching only), `Category_x`, `Category_y`, `IdX`, `IdY`,
`Interaction_Type`.

Composition of the 21,404 kept pairs — top classes: `lncRNA-miRNA` 6,175,
`miRNA-ncRNA` 6,033, `miRNA-pseudo` 5,416, `snoRNA-snoRNA` 673,
`lncRNA-ncRNA` 479. Verified clean: no duplicate rows, no reverse-orientation
duplicates, one self-pair `(A,A)`.

## Caveats before using this as a test set

* **Positives only.** Like RNAInterSeq, RNA-KG lists reported interactions and
  provides no experimentally established non-interacting pairs. Negatives still
  have to be constructed, with the usual care over degree-matched sampling.
* **Pair-level, not entity-level, disjointness.** These records contain no
  interaction RNAInterSeq has, but 11,815 of their 14,700 sequences *do* occur
  in RNAInterSeq. A model that memorises per-RNA propensity can still score
  well. For a stricter split, filter to pairs whose partners are both in
  `ncrna_entities_non_overlapping.pkl` — that leaves **647** pairs, a much
  harder and much smaller benchmark.
* **Exact sequence matching.** A transcript differing by one base — a different
  assembly patch, a SNP-bearing haplotype, a trimmed UTR — counts as absent
  from RNAInterSeq and is retained, so a small residue of near-duplicates is
  expected. An alignment-based filter (cf. `similarity_computation.py`) would
  remove it.
* **No species constraint.** Matching is on sequence alone; RNA-KG's ncRNA
  export carries no taxonomy column.
