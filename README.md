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
| mRNA-associated | 7,332,037 | 98.4 % |
| ncRNA–ncRNA | 119,365 | 1.6 % |
| **total** | **7,451,402** | |

Built from 9,429,404 curated RNAInter interactions; those whose partners could
not both be sequence-resolved are dropped. The result covers **174,319 unique
RNA entities** (135,183 mRNA, 39,136 ncRNA) across 56 species (NCBI Taxonomy
identifiers; RNAInter lists dog under two names, *Canis familiaris* and
*Canis lupus familiaris*, both taxid 9615) and 14 RNA biotypes, carrying 796,477
distinct transcript sequences. Because entities keep
every recovered isoform rather than one representative, expanding interactions
to transcript level gives 91,665,313 pairs.

## Pipeline

Run the notebooks in order from the repository root.

| step | notebook / module | what it does |
|---|---|---|
| 1 | `1_RNAInter_sequence_map.ipynb` | filter RNAInter records, standardise species to NCBI taxids, parse identifiers into namespace buckets, emit the non-redundant RNA-entity table |
| 2a | `2_sequence_map_ncRNA.ipynb` | resolve ncRNA entities to transcripts via RNAcentral, miRBase, circBase, NCBI/RefSeq and Ensembl |
| 2b | `2_sequence_map_mRNA.ipynb` | resolve mRNA entities to RefSeq/Ensembl transcripts |
| 2c | `2c_refetch_refseq_failures.py` | re-fetch the RefSeq accessions lost to failed efetch batches in the original run and rebuild the mRNA entity table (see below) |
| 3 | `3_interactions_map.py` (script form of `3_Interactions_map.ipynb`) | project sequences and metadata back onto both partners; keep interactions where both are sequence-resolved; the transcript accessions of mRNA entities are retained |
| — | `similarity_computation.py` | pairwise global identity between the isoforms of an entity (batch job) |
| — | [`rnacentral-etl`](https://github.com/AnacletoLAB/rnacentral-etl) | submodule: DuckDB ETL over the RNAcentral FTP distribution |

Statistics and figures: `ncRNAInter_stats.ipynb`, `mRNAInter_stats.ipynb`
(the `*_union.ipynb` copies read the newest release);
`paper_figures/05_paper_numbers.py` recomputes every number quoted in the
manuscript from the current release (`data/paper_numbers_<date>.json`) and is
the authoritative source of those numbers, the notebooks being the record of
the original analysis; and
`zenodo_release/build_release.py` assembles the release folder (pickles, TSV,
FASTA, checksums).
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
   → `2_sequence_map_mRNA.ipynb` → `2c_refetch_refseq_failures.py` (only needed
   to repair a run whose efetch batches failed; the fetch cell of 2b now retries
   failed batches one accession at a time and reads every
   `data/refseq_sequences_cache*.json`) → `3_interactions_map.py` (writes the
   release, `data/interactions_with_sequences_RNA_RNA_union_<date>.pkl`; the figure
   scripts pick the newest `*_union_*.pkl`) → `mRNAInter_stats_union.ipynb` /
   `ncRNAInter_stats_union.ipynb` (evidence, GO/Rfam statistics; the mRNA
   notebook dumps its multi-transcript table to `data/df_step3_union.pkl` for
   `similarity_computation.py`, a batch job whose input/output paths are set at
   the top of the script, and reads the `*_with_similarity.pkl` result back) →
   the `paper_figures/` scripts (`00b`, `01`–`05`) → `gtex_ncrna_expression/` (own README)
   → `zenodo_release/build_release.py`.
7. **Resources**: the release table is a 6.0 GB pickle of Python objects and
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
fallback for entities left unresolved; making the two symmetric gives a symbol
match to 181,262 mRNA entities: 2,844 of them had no RefSeq accession from the
NCBI Gene step (152 were already resolved through Ensembl, 2,692 are resolved
only here) and 2,729 gain further accessions (8,092 accession slots added on top
of the NCBI Gene step). The three steps assign at least one accession to
194,559 + 318 + 2,692 = 197,569 mRNA entities, and all of them are
sequence-resolved in the current release. In the March 2026 run 29 batch
requests to NCBI efetch (1,450 accessions) had failed and were never retried, so
that 78 entities had no sequence and 491 lacked 1,352 isoforms (197,491
sequence-resolved entities); `2c_refetch_refseq_failures.py` re-fetched them on
18-09-2026 into `data/refseq_sequences_cache_refetch_18-09-2026.json`, rebuilt
the mRNA entity table and the release was regenerated (+78 entities, +2,064
interactions, +1,494 sequences).

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
