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
RNA entities** (135,116 mRNA, 39,132 ncRNA) across 57 species and 14 RNA
biotypes, carrying 795,593 distinct transcript sequences. Because entities keep
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
appear in the manuscript: frequency panels report relative frequencies rather
than absolute counts, and two plotting defects present in the notebooks (GO
domain bars drawn in `Counter` insertion order against a fixed colour list, and
hardcoded interaction counts in the GO boxplot labels) are corrected there.

### Sub-analyses

* **`gtex_ncrna_expression/`** — links the human ncRNA–ncRNA interactions to
  GTEx tissue expression: entity→GENCODE locus resolution, sequence-level
  transcript confirmation against Ensembl, Yanai's τ, and interaction-level
  co-expression. See its own README for the step order.
* **`rnakg_non_overlapping/`** — builds the subset of
  [RNA-KG](https://github.com/AnacletoLAB/RNA-KG) that RNAInterSeq does *not*
  contain, for use as a leakage-free evaluation set. See its README.

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
`matplotlib`, `seaborn`, `goatools`, `joblib`, `duckdb`, `ete3`.

External inputs, none of them in this repository: the RNAInter v4.0 RNA–RNA
download, NCBI `gene2refseq` and `gene2ensembl`, miRBase v22 `miRNA.dat`,
circBase FASTA exports, an RNAcentral release (v25) for the ETL, and
`go-basic.obo`.

## Data availability

The dataset, its GTEx annotation and the RNA-KG non-overlapping evaluation set
are published on Zenodo (DOI to be added). `.gitignore` keeps all `*.pkl`
tables, the `data/` tree, cached API responses and the assembled release folder
out of the repository.

## Citation

Nicolini M., Ferreyra S., Paccanaro A., Robinson P.N., Casiraghi E.,
Valentini G. *RNAInterSeq: RNA interaction sequence-based database.*
Manuscript in preparation.
