# GTEx tissue expression for the ncRNA-ncRNA RNAInterSEQ set

Side project that links the human ncRNAs in
`data/interactions_with_sequences_RNA_RNA_ncrna_ncrna_19-05-2026.pkl` (the
ncRNA-ncRNA RNAInterSEQ release) to their tissue expression profile in
[GTEx](https://gtexportal.org/) (Genotype-Tissue Expression). GTEx profiles
**postmortem human donor tissue only** -- there is no other-species axis in
the API (confirmed via its 54-tissue `tissueSiteDetailId` vocabulary), so
every non-human entity in the interaction set is out of scope by construction.

## Database scope

The source file contains **119,365 ncRNA-ncRNA interaction pairs** across
30 species-pair combinations:

| Species pair | Pairs | Share |
|---|---|---|
| Homo sapiens × Homo sapiens | 81,175 | 68.0 % |
| Mus musculus × Mus musculus | 35,236 | 29.5 % |
| Other / cross-species | 2,954 | 2.5 % |

Only the **81,175 human-human pairs** are annotated in this pipeline.

## Pipeline

Run scripts in this order (all with the `vGenRNA` conda env):

```
python 01_build_entities.py
python 02_fetch_gtex.py
python 05_filter_transcripts_by_sequence.py
python 03_compute_stats.py
python 06_annotate_interactions.py
python 07_tissue_figures.py
python 04_generate_report.py
```

Support modules (not run directly):

| Module | Purpose |
|---|---|
| `gtex_api.py` | GTEx Portal API v2 client (batched, cached, retry-aware) |
| `ensembl_api.py` | Ensembl REST API client for transcript sequences |
| `id_resolvers.py` | Converts NCBI/miRBase/circBase IDs → GTEx-queryable gene symbol |

### Step descriptions

1. **`01_build_entities.py`** — extracts the 10,068 unique human ncRNA
   entities and their ID buckets from the interaction table.
   → `data/ncrna_entities_human.pkl`

2. **`02_fetch_gtex.py`** — resolves each entity to one or more GENCODE
   genes (Ensembl ID, NCBI→Ensembl, symbol, miRBase precursor, circBase
   host gene, each tried at GENCODE v39/`gtex_v10` then v26/`gtex_v8`),
   then fetches median gene-level and transcript-level expression across
   all 54 GTEx tissues for every matched gene.
   → `data/gtex_gene_matches.pkl`, `data/gtex_gene_expression.pkl`,
      `data/gtex_transcript_expression.pkl`

3. **`05_filter_transcripts_by_sequence.py`** — fetches each candidate
   transcript's spliced cDNA from Ensembl and keeps only the ones whose
   sequence matches (exact or substring containment in either direction) a
   sequence RNAInterSEQ already has for that entity. Produces the
   sequence-confirmed transcript-expression table used downstream.
   → `data/ncrna_locus_transcript_matches.pkl`,
      `data/gtex_transcript_expression_sequence_confirmed.pkl`

4. **`03_compute_stats.py`** — coverage funnel, match-rate by ncRNA
   category, resolution-method/GENCODE-version breakdown, Yanai tau,
   top tissue-specific ncRNAs.
   → `data/ncrna_entity_locus_summary.pkl`, `figures/Fig1-Fig5.pdf`

5. **`06_annotate_interactions.py`** — joins all entity-level GTEx metadata
   onto the 81,175 human-human interaction rows. Adds 26 columns per
   interaction: per-entity confirmed transcripts, tau, peak tissue, peak
   TPM, 54-tissue TPM vector, biological system lists, plus shared-tissue
   and shared-system summary columns for each pair.
   → `data/human_ncrna_interactions_gtex_annotated.pkl`  *(final output)*

6. **`07_tissue_figures.py`** — interaction-level figures (Figs 6b–13).
   → `figures/Fig6_step5_funnel.pdf` through `figures/Fig13_tissue_system_coexpression.pdf`

7. **`04_generate_report.py`** — renders `gtex_ncrna_expression_report.tex`
   from all pickles and compiles it to PDF (uses `tectonic` from the
   `texreport` conda env if available).
   → `gtex_ncrna_expression_report.pdf`

## Output columns in the annotated interaction table

The final table (`human_ncrna_interactions_gtex_annotated.pkl`, 81,175 × 59)
extends the original 33 RNAInterSEQ columns with 26 GTEx annotation columns:

| Column | Description |
|---|---|
| `ent{1,2}_gtex_matched` | Entity resolved to ≥1 GENCODE gene in GTEx |
| `ent{1,2}_seq_confirmed` | ≥1 transcript sequence-confirmed against RNAInterSEQ |
| `ent{1,2}_confirmed_transcripts` | `{gencode_id: [transcript_ids]}` dict |
| `ent{1,2}_tau` | Yanai tissue-specificity index (0 = ubiquitous, 1 = single-tissue) |
| `ent{1,2}_peak_tissue` | GTEx tissue ID with highest confirmed median TPM |
| `ent{1,2}_peak_tpm` | Median TPM at peak tissue |
| `ent{1,2}_n_expressed_tissues` | Number of 54 tissues with max-TPM > 0 |
| `ent{1,2}_tissue_tpm` | `{tissue_id: max_median_tpm}` dict over all 54 tissues |
| `ent{1,2}_tissue_systems` | Sorted list of biological systems with expression > 0 |
| `both_gtex_matched` | Both partners have a GTEx match |
| `both_seq_confirmed` | Both partners have sequence-confirmed transcripts |
| `shared_expressed_tissues` | Sorted list of tissues where both TPM > 0 |
| `n_shared_expressed_tissues` | Count of shared expressed tissues |
| `shared_tissue_systems` | Sorted list of the biological systems to which the shared expressed tissues belong (non-empty only if the partners share a tissue) |
| `n_shared_tissue_systems` | Count of the above |
| `shared_systems` | Sorted list of biological systems in which both partners are detected (intersection of `ent1_tissue_systems` and `ent2_tissue_systems`; the partners may be detected in different tissues of the same system); sequence-confirmed pairs only |
| `n_shared_systems` | Count of shared biological systems |

### Tissue system grouping

The 54 GTEx tissues are collapsed to 8 biological systems:

| System | # tissues | Key tissues |
|---|---|---|
| CNS | 13 | All brain sub-regions + spinal cord |
| Peripheral nervous | 2 | Tibial nerve, pituitary |
| Cardiovascular | 5 | Heart (atrium, ventricle), arteries |
| GI | 10 | Colon, esophagus, liver, pancreas, stomach, small intestine |
| Reproductive | 9 | Testis, ovary, uterus, cervix, prostate, breast, vagina, fallopian tube |
| Immune | 3 | Whole blood, spleen, EBV-transformed lymphocytes |
| Muscle | 1 | Skeletal muscle |
| Other | 11 | Adipose, adrenal, bladder, fibroblasts, kidney, lung, skin, thyroid |

## Headline results (19-05-2026 release, GTEx v10 primary / v8 fallback)

**Entity-level (step 2)**
- 10,068 unique human ncRNA entities; **7,772 (77.2%)** resolved to ≥1
  GENCODE gene known to GTEx.
- Match rate by category: snoRNA 98.3%, lncRNA 97.7%, scaRNA 100%,
  circRNA 93.9%, pseudogene 89.1%, rRNA 85.7%, snRNA 81.0%, tRNA 80.0%,
  **miRNA 34.0%** (structural gap: bulk RNA-seq cannot capture short mature
  miRNA products).
- 753 entity-locus pairs (470 entities) only resolved at the v26/`gtex_v8`
  fallback — GENCODE v39 has since pruned these loci.
- Tissue specificity: median tau = 0.87 (most ncRNAs are highly
  tissue-selective); peak tissue skew is driven by Testis (most
  transcriptionally permissive GTEx tissue).

**Interaction-level (step 6)**
- 81,175 human-human ncRNA-ncRNA interaction pairs annotated.
- 54,746 (67.4%) have both partners GTEx-matched.
- 31,264 (38.5%) have both partners sequence-confirmed.
- 2,375 (2.9% of all pairs; 7.6% of the 31,264 sequence-confirmed pairs) share ≥1 expressed GTEx tissue (median TPM > 0 in both).
- By biological system: Reproductive leads (2,168 co-expressed pairs),
  driven by the broad transcriptional activity of testis; CNS is 4th
  (1,492 pairs), consistent with the large proportion of brain-enriched
  lncRNAs in the interaction set.

## Notes

- Run with `/home/mnicolini/miniconda3/envs/vGenRNA/bin/python <script>.py`
  (the only local conda env with `pandas`, `requests`, `matplotlib`).
- LaTeX report compiled with `tectonic` from a dedicated `texreport` conda
  env (`conda create -n texreport -c conda-forge tectonic`);
  `04_generate_report.py` finds it automatically.
- GTEx/mygene.info/Ensembl API responses are cached under `cache/`.

See `gtex_ncrna_expression_report.pdf` for the full formatted report with figures.
