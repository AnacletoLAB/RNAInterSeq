# Paper figures

Scripts that produce the figure panels used in the manuscript, from the current
release (`data/interactions_with_sequences_RNA_RNA_union_*.pkl`) and the entity
tables derived from it.

They exist as standalone scripts, rather than as cells of
`mRNAInter_stats.ipynb` / `ncRNAInter_stats.ipynb`, because the published panels
differ from the notebook versions in three respects:

1. **Relative frequencies instead of absolute counts.** Every panel that is a
   frequency now reports a fraction: interaction-biotype and species
   distributions, transcript length, GO domains, transcript-count bins and
   pairwise identity. Panels that are not frequencies — medians, boxplots,
   per-species entropy and the identity-versus-pairs scatter — keep their
   original scale.
2. **Two plotting defects fixed.** GO-domain bars were drawn from a `Counter`
   in insertion order against a fixed colour list, so labels and colours were
   not guaranteed to match BP/MF/CC; and the interaction counts in the GO
   boxplot labels were hardcoded and matched neither subset. Both are now
   derived from the data.
3. **Tight bounding box on export.** `predicted_evidence_type` called neither
   `tight_layout()` nor `bbox_inches="tight"`, which truncated the longest
   method labels on the left margin.

The notebooks are left as the record of the original analysis; these scripts are
the provenance of the panels that appear in the paper.

## Order

```
python 00_ncrna_similarity.py     # ncRNA transcript-similarity table (only if absent)
python 01_interaction_figures.py  # interaction-level panels, both subsets
python 02_entity_figures.py       # entity-level panels, both subsets
python 03_flow_diagram.py <out.png>
```

Steps 01–02 write into `data/mrna_stats_union/` and `data/ncrna_stats_union/`;
copy the PDFs into the manuscript's `fig/mrna_stats/` and `fig/ncrna_stats/`.
Step 03 takes its output path as an argument.

`00` reuses the per-entity results of a previous run for every entity whose
transcript set is unchanged and recomputes only the rest. It parallelises over
individual transcript *pairs* rather than whole entities: entity-level
parallelism serialises on single very large entities — TTN alone accounts for
6.0 of the 9.3 $\times 10^{12}$ alignment cells in a full mRNA run — and the two
are mathematically equivalent, since every pairwise identity is independent.

## Requirements

`pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `biopython`, and
`go-basic.obo` in the repository root (parsed directly, so `goatools` is not
needed).
