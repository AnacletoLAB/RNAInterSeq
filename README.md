# RNAInterSeq — sequence-resolved RNA–RNA interactions

Pipeline that turns the RNA–RNA interaction records of
[RNAInter v4.0](http://www.rnainter.org/) into a **sequence-resolved** dataset:
every interaction whose two partners could both be mapped to at least one
nucleotide sequence, with identifiers, RNA biotypes, taxonomy, Gene Ontology
terms and Rfam families attached, and — for the human non-coding subset — GTEx
tissue-expression annotation.

This repository holds the **code only**. The dataset itself is distributed via
Zenodo (DOI to be added) and is deliberately excluded from version control.

## The dataset it produces

| | interactions | share |
|---|---:|---:|
| mRNA-associated | 7,329,973 | 98.4 % |
| ncRNA–ncRNA | 119,365 | 1.6 % |
| **total** | **7,449,338** | |

Built from 9,429,404 curated RNAInter interactions; those whose partners could
not both be sequence-resolved are dropped. The result covers **174,248 unique
RNA entities** (135,116 mRNA, 39,132 ncRNA) across 56 species (NCBI Taxonomy
identifiers; RNAInter lists dog under two names, *Canis familiaris* and
*Canis lupus familiaris*, both taxid 9615) and 14 RNA biotypes, carrying 795,593
distinct transcript sequences. Because entities keep
every recovered isoform rather than one representative, expanding interactions
to transcript level gives 91,559,557 pairs.

## Pipeline

Run the notebooks in order from the repository root.

| step | notebook / module | what it does |
|---|---|---|
| 1 | `1_RNAInter_sequence_map.ipynb` | filter RNAInter records, standardise species to NCBI taxids, parse identifiers into namespace buckets, emit the non-redundant RNA-entity table |
| 2a | `2_sequence_map_ncRNA.ipynb` | resolve ncRNA entities to transcripts via RNAcentral, miRBase, circBase, NCBI/RefSeq and Ensembl |
| 2b | `2_sequence_map_mRNA.ipynb` | resolve mRNA entities to RefSeq/Ensembl transcripts |
| 3 | `3_Interactions_map.ipynb` | project sequences and metadata back onto both partners; keep interactions where both are sequence-resolved |
| — | `similarity_computation.py` | pairwise global identity between the isoforms of an entity (batch job) |
| — | [`rnacentral-etl`](https://github.com/AnacletoLAB/rnacentral-etl) | submodule: DuckDB ETL over the RNAcentral FTP distribution |

Statistics and figures: `ncRNAInter_stats.ipynb`, `mRNAInter_stats.ipynb`.
The `*_union.ipynb` copies are the same notebooks pointed at the current
release, used to regenerate the figures after the mapping change described
below. `paper_figures/` holds the scripts that produce the panels as they
appear in the manuscript (run order in `paper_figures/README.md`): frequency
panels report relative frequencies rather than absolute counts, two plotting
defects present in the notebooks (GO domain bars drawn in `Counter` insertion
order against a fixed colour list, and hardcoded interaction counts in the GO
boxplot labels) are corrected there, and the entity-level panels are built
from `paper_figures/00b_entity_tables_by_key.py`, which keys RNA-entities by
the pipeline key (symbol, biotype, taxid, identifiers) — the notebooks keyed
them by species name without identifiers, which double-counted dog.

### Reproducing the pipeline

1. **Clone with the submodule**: `git clone --recurse-submodules …` (or
   `git submodule update --init`). The notebooks import `rnacentral.py` from
   `rnacentral-etl/` via `sys.path.append("rnacentral-etl")`.
2. **Inputs** (not distributed here), placed as the notebooks expect them:

   | file / folder | content | used by |
   |---|---|---|
   | `data/RNAInter_RNA_RNA.txt` | RNAInter v4.0 RNA–RNA download (tab-separated, 9,484,609 rows) | 1 |
   | `data/gene2refseq`, `data/gene2ensembl` | NCBI Gene FTP tables, uncompressed | 2a, 2b |
   | `data/miRNA.dat` | miRBase v22 EMBL-format release | 2a |
   | `data/circdata/` | circBase FASTA exports, one per taxid (`circBase_<taxid>.fa`), see step 4 | 2a |
   | `rnacentral_data/current_release/` | RNAcentral release 25 FTP directories `sequences/`, `id_mapping/`, `md5/`, `rfam/`, `go_annotations/` (rsync command in `rnacentral-etl/README.md`) | ETL |
   | `go-basic.obo` | Gene Ontology basic release, repository root | stats notebooks, `paper_figures` |

3. **Build the RNAcentral DuckDB once** (≈ tens of GB, several hours):
   `QueryRnaCentral` opens `rnacentral-etl/rnacentral.duckdb` read-only and
   assumes the tables exist. To create them, run once from the repository root
   with `read_only` temporarily set to `False` in `rnacentral.py`:

   ```python
   import sys; sys.path.append("rnacentral-etl")
   from rnacentral import QueryRnaCentral
   q = QueryRnaCentral(base_dir="rnacentral_data/current_release",
                       db_path="rnacentral-etl/rnacentral.duckdb")
   q._init_tables(); q.close()
   ```
4. **circBase is a manual step**: notebook 2a writes one identifier list per
   taxid to `data/circdata/<taxid>.txt`; submit each to the circBase batch
   search (https://www.circbase.org), save the returned FASTA as
   `data/circdata/circBase_<taxid>.fa`, then continue the notebook.
5. **External services** called live, with local caches: NCBI E-utilities
   (`Entrez.email` in 2a is the authors' — put your own; 2b calls `efetch`
   through `requests`, and `data/refseq_sequences_cache.json` caches the
   ≈945k RefSeq fetches, which take several hours), Ensembl REST, MyGene.info, and, for the GTEx
   sub-analysis, the GTEx Portal API v2. `ete3` downloads the NCBI taxonomy
   dump on first use.
6. **Run order**: `1_RNAInter_sequence_map.ipynb` → `2_sequence_map_ncRNA.ipynb`
   → `2_sequence_map_mRNA.ipynb` → `3_Interactions_map.ipynb` (writes the
   release, `data/interactions_with_sequences_RNA_RNA_<date>.pkl`; the figure
   scripts pick the newest `*_union_*.pkl`) → `mRNAInter_stats_union.ipynb` /
   `ncRNAInter_stats_union.ipynb` (evidence, GO/Rfam statistics; the mRNA
   notebook dumps its multi-transcript table to `data/df_step3_union.pkl` for
   `similarity_computation.py`, a batch job whose input/output paths are set at
   the top of the script, and reads the `*_with_similarity.pkl` result back) →
   the `paper_figures/` scripts → `gtex_ncrna_expression/` (own README).
7. **Resources**: the release table is a 5.8 GB pickle of Python objects and
   takes ≈ 15 GB of RAM to load (more for the per-row expansions in the stats
   notebooks); the pairwise-identity step is CPU-bound
   (≈ 10<sup>13</sup> alignment cells for the mRNA set, minutes on 48 cores
   with the pair-parallel scripts, hours with the entity-parallel batch job).

### Sub-analyses

* **`gtex_ncrna_expression/`** — links the human ncRNA–ncRNA interactions to
  GTEx tissue expression: entity→GENCODE locus resolution, sequence-level
  transcript confirmation against Ensembl, Yanai's τ, and interaction-level
  co-expression. See its own README for the step order.

## Identifier resolution

Interactors reach the pipeline as heterogeneous, often concatenated identifier
strings. They are parsed into per-source buckets (NCBI Gene, Ensembl, miRBase,
circBase, RNAcentral, GenBank, tRFdb) and **kept in their own namespace** rather
than collapsed into one universal id. Each namespace is then resolved
independently and the resulting transcript sets are **merged by union**, since
the routes are complementary rather than redundant. Every match is constrained
by the entity's taxonomy id and checked for RNA-biotype consistency.

Both the ncRNA and the mRNA pipeline now apply gene-symbol mapping to **all**
entities and merge by union. The mRNA side previously used it only as a
fallback for entities left unresolved; making the two symmetric resolves 2,844
additional entities and adds 8,092 transcript accessions.

## Requirements

Python ≥ 3.9 with `pandas`, `numpy`, `biopython`, `requests`, `tqdm`,
`matplotlib`, `seaborn`, `scipy`, `goatools`, `joblib`, `duckdb`, `ete3`,
`pyarrow`.

External inputs, none of them in this repository: the RNAInter v4.0 RNA–RNA
download, NCBI `gene2refseq` and `gene2ensembl`, miRBase v22 `miRNA.dat`,
circBase FASTA exports, an RNAcentral release (v25) for the ETL, and
`go-basic.obo`.

## Data availability

The dataset and its GTEx annotation are published on Zenodo (DOI to be
added). `.gitignore` keeps all `*.pkl`
tables, the `data/` tree, cached API responses and the assembled release folder
out of the repository.

## Citation

Nicolini M., Ferreyra S., Paccanaro A., Robinson P.N., Casiraghi E.,
Valentini G. *RNAInterSeq: RNA interaction sequence-based database.*
Manuscript in preparation.
