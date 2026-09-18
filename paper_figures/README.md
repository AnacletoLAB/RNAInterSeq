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
4. **Top and right frame lines removed.** Both scripts set
   `axes.spines.top/right = False`, but called `sns.set_style("white")`
   afterwards, which turns all four spines back on. `set_style` now runs first,
   so the panels match the notebook originals, which never called it.

`00b` replaces the entity tables built inside the stats notebooks
(`data/entities_union_*.pkl`, `df_step3_*_with_similarity.pkl`), which keyed
entities by (symbol, biotype, species *name*) and dropped the identifier
component (`make_entity_key` received the `Raw_ID` string, not a list). That
double-counted dog, listed by RNAInter under two names for taxid 9615, and
merged entities differing only by identifier: 176,444 "entities" for the
mRNA-associated subset of a dataset with 174,248 entities in total. The
by-key tables give 168,329 and 14,631 entities; per-entity alignments are
reused wherever the transcript set is unchanged. `01` also merges the two dog
names before counting species pairs.

`experimental_evidence_type.pdf` is produced by `06_evidence_panel.py` from the
evidence-method categories computed in `05_paper_numbers.py`, which ports the
hand-curated method-to-category map of `mRNAInter_stats.ipynb` (cells 24–26)
verbatim; the ncRNA–ncRNA panel and Table 4 of the manuscript reproduce the
notebook values exactly, which validates the port.

`05_paper_numbers.py` recomputes every number quoted in the manuscript from the
newest release (`data/paper_numbers_<date>.json`): interaction, entity, species
and sequence counts, transcript-level pairs, evidence categories, GO/Rfam
coverage, biotype/species-pair landscapes, per-entity length, multiplicity and
identity statistics. `03`/`04` read their counts from that file and from the
processed mRNA table, so they follow the pipeline instead of hardcoding values.

The notebooks are left as the record of the original analysis; these scripts are
the provenance of the panels that appear in the paper.

## Order

```
python 00_ncrna_similarity.py     # ncRNA transcript-similarity table (only if absent)
python 00b_entity_tables_by_key.py  # entity tables keyed by (symbol, biotype, taxid, ids) + similarity
python 01_interaction_figures.py  # interaction-level panels, both subsets
python 02_entity_figures.py       # entity-level panels, both subsets
python 05_paper_numbers.py        # data/paper_numbers_<date>.json (needed by 03 and 06)
python 06_evidence_panel.py       # experimental_evidence_type.pdf, both subsets
python 03_flow_diagram.py <out.png>   # fig/flow_RNA_RNA_inter.png
python 04_entity_diagram.py <out.png> # fig/RNA_entities_plot.png
```

`00b` and `05` run under the `rnacentral` env (pandas ≥ 2.3, numpy 2, Biopython);
the plotting scripts need matplotlib/seaborn with numpy 2 (the `vGenRNA` env), since
the release pickles are not readable with numpy 1.x.

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
