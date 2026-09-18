# RNAInterSeq — a sequence-resolved dataset of RNA–RNA interactions

**7,451,402 RNA–RNA interactions in which both partners carry nucleotide
sequences**, derived from RNAInter v4.0 and covering **both mRNA-associated and
ncRNA–ncRNA interactions**, with identifiers, RNA biotypes, taxonomy, Gene
Ontology terms and Rfam families attached — plus GTEx tissue-expression
annotation for the human non-coding subset.

Code: <https://github.com/AnacletoLAB/RNAInter_sequenceParser>

## Contents

| file | rows | size | content |
|---|---:|---:|---|
| `data/RNAInterSeq_interactions_full.pkl` | 7,451,402 | 6.0 GB | **the dataset** — both interaction classes, all columns (pandas pickle) |
| `data/RNAInterSeq_human_ncRNA_GTEx.pkl` | 81,175 | 97 MB | human ncRNA–ncRNA subset + 26 GTEx columns (pandas pickle) |
| `data/RNAInterSeq_interactions.tsv.gz` | 7,451,402 | 133 MB | one row per interaction: RNAInter fields, evidence, the two entity ids |
| `data/RNAInterSeq_entities.tsv.gz` | 174,319 | 10 MB | one row per RNA entity: key, biotype, taxonomy, identifiers, transcript accessions, GO, Rfam |
| `data/RNAInterSeq_transcripts.tsv.gz` | 915,541 | 7 MB | one row per (entity, transcript): index, accession, source, length |
| `data/RNAInterSeq_transcripts.fasta.gz` | 915,541 | 346 MB | the transcript sequences, header `entity_id\|index\|accession` |
| `data/RNAInterSeq_human_ncRNA_GTEx.tsv.gz` | 81,175 | 9 MB | the GTEx annotation columns, flattened (dict-valued fields as JSON) |

The two pickles are pandas DataFrames (`pandas.read_pickle(path)`, pandas ≥ 2.3,
numpy ≥ 2) and carry the complete, nested content. The TSV/FASTA files carry the
same content in a normalised, language-neutral form: `entity_id` (`E` + six
digits) is the join key between interactions, entities and transcripts;
list-valued fields are `;`-separated, dict-valued fields are JSON strings, and
a missing gene symbol (NaN in the pickle) is an empty field. MD5 checksums of
every file are in `CHECKSUMS.md5` / `CHECKSUMS.json`.

## Composition of the main file

| | interactions | share |
|---|---:|---:|
| **mRNA-associated** | 7,332,037 | 98.4 % |
| **ncRNA–ncRNA** | 119,365 | 1.6 % |
| **total** | **7,451,402** | |

Select either class with
`df[df.Category1.eq("mRNA") | df.Category2.eq("mRNA")]` and its negation.

**174,319 unique RNA entities** — 135,183 mRNA and 39,136 non-coding
(piRNA 14,713 · lncRNA 11,098 · miRNA 10,941 · snoRNA 882 · pseudo 881 ·
snRNA 260 · circRNA 189 · rRNA 111 · scaRNA 22 · tRNA 20 · ncRNA 17 ·
ribozyme 1 · processed_transcript 1) — across **56 species** (NCBI Taxonomy
identifiers; 57 organism names, dog being listed as both *Canis familiaris* and
*Canis lupus familiaris*, taxid 9615) and 14 RNA
biotypes, carrying **796,477 distinct transcript sequences**
(length: min 11, median 3,130, mean 3,949, max 205,012 nt).

Entities keep **every isoform recovered**, not one representative, so expanding
interactions to transcript level yields **91,665,313** transcript pairs
(90,972,400 mRNA-associated, 692,913 ncRNA–ncRNA).

Species coverage is uneven and biased toward well-studied organisms; within the
ncRNA–ncRNA class 68.0 % of pairs are human–human and 29.5 % mouse–mouse.

## Columns (37)

**Interaction and partners**

| column | content |
|---|---|
| `RNAInterID` | RNAInter v4.0 record id — *not unique*, 37 rows repeat |
| `Interactor1.Symbol`, `Interactor2.Symbol` | gene symbol per partner |
| `Category1`, `Category2` | RNA biotype per partner |
| `Species1`, `Species2` | organism name |
| `Species1_TaxID`, `Species2_TaxID` | NCBI Taxonomy id |
| `Raw_ID1`, `Raw_ID2` | original RNAInter identifier strings |
| `key1`, `key2` | composite entity key `(symbol, biotype, taxid, ids)` |

**Sequences and identifiers** — `ent1_*` / `ent2_*`, one set per partner

| column | content |
|---|---|
| `ent*_sequences` | **list** of transcript sequences, all isoforms recovered |
| `ent*_refseq_accessions` | mRNA entities: RefSeq accessions of the transcripts, in the order of the RefSeq sequences in `ent*_sequences` (the Ensembl-derived sequences, if any, follow them); NaN for ncRNA entities, whose RNAcentral accessions sit in `ent*_ID_buckets` |
| `ent*_ensembl_ids` | mRNA entities: Ensembl gene ids expanded to protein-coding transcripts through the Ensembl REST API; NaN for ncRNA entities |
| `ent*_ID_buckets` | identifiers per source: Ensembl, NCBI, miRBase, circBase, RNAcentral, GenBank, tRFdb |
| `ent*_Raw_key` | normalised identifier tuple |
| `ent*_typecheck` | `False` where the ncRNA biotype could not be confirmed |

Sequence alphabet differs by provenance: ncRNA sequences are RNA (`AUCG`),
mRNA sequences come raw from RefSeq in DNA alphabet. Normalise with
`s.upper().replace("T","U")` before comparing across the two.

**Functional annotation**

| column | content |
|---|---|
| `ent*_GO_RNAcentral` | GO term ids from RNAcentral |
| `ent*_GO_raw` | GO terms from MyGene.info, split BP / MF / CC |
| `ent*_mygene_raw` | full MyGene.info record |
| `ent*_Rfam` | Rfam family accessions (e.g. `RF00272`) |

**Evidence** — inherited verbatim from RNAInter

| column | content | rows |
|---|---|---:|
| `score` | RNAInter confidence score | |
| `strong` | low-throughput, interaction-specific assays | 409,774 |
| `weak` | transcriptome-wide screens | 1,335,630 |
| `predict` | computational support | 6,421,117 |

The three are **not mutually exclusive**; a record may carry several.
*Validated* means `strong ∪ weak` = 1,665,399 — not a fourth independent class.
Per subset: mRNA-associated 378,590 / 1,310,988 / 6,353,227 (validated
1,610,445); ncRNA–ncRNA 31,184 / 24,642 / 67,890 (validated 54,954).

## GTEx annotation

`RNAInterSeq_human_ncRNA_GTEx.pkl` adds 26 columns to the **81,175 human–human
ncRNA–ncRNA** interactions: per partner the sequence-confirmed GENCODE
transcripts, Yanai's τ, peak tissue and TPM, the full 54-tissue TPM vector and
the biological systems in which it is detected, plus the tissues shared by both
partners (`shared_expressed_tissues`), the systems those shared tissues belong to
(`shared_tissue_systems`) and the systems in which both partners are detected,
possibly in different tissues (`shared_systems`). The subset was cut from the
previous release and carries its 33 base columns (no accession columns).

It covers only that subset by construction — GTEx profiles postmortem *Homo
sapiens* tissue and its API has no species axis. Coverage inside the subset is
partial: both partners matched to a GTEx gene in 54,746 pairs (67.4 %), both
sequence-confirmed in 31,264 (38.5 %), sharing ≥1 expressed tissue in 2,375
(2.9 %) and ≥1 biological system in 2,499 (3.1 %). That last figure is a **conservative lower bound**: in 91.9 % of
sequence-confirmed pairs at least one partner is undetected across all 54
tissues, chiefly mature miRNAs and snoRNAs that poly(A)-selected bulk RNA-seq
does not capture — only 0.5 % are cases of both partners detected in
non-overlapping tissues.

Join back on `["RNAInterID", "key1", "key2"]`, not on `RNAInterID` alone.

One mRNA entity (NCBI Gene 39571, *Drosophila melanogaster*, `E000001` in the
TSV files) has no gene symbol in RNAInter: its key in the pickles is
`(nan, 'mRNA', '7227', ('NCBI:39571',))`, which cannot be matched by value
(`NaN != NaN`); select it with `df.key1.apply(lambda k: k[3] == ('NCBI:39571',))`
or through `entity_id` in the TSV files, where the symbol is empty.

## How it was built

RNAInter v4.0 → filtering of ambiguous, non-RNA and unmappable records →
species standardisation to NCBI taxids → identifier parsing into per-source
namespaces → non-redundant RNA-entity table → transcript-sequence mapping
(ncRNA via RNAcentral, miRBase, circBase, NCBI/RefSeq, Ensembl; mRNA via
RefSeq `gene2refseq`, gene symbol and the Ensembl REST API) → GO and Rfam
annotation → projection back onto both partners, keeping only interactions
where both are sequence-resolved (9,429,404 → 7,451,402).

Every namespace is resolved independently under the entity's taxonomy id and
checked for biotype consistency; the resulting transcript sets are merged by
**union**, in both the ncRNA and the mRNA pipeline.

## Limitations

* **Positives only.** Reported interactions; no experimentally established
  non-interacting pairs. Negatives must be constructed deliberately —
  random pairing produces false negatives, degree-matched sampling introduces
  its own bias.
* **Gene-level records, transcript-level sequences.** RNAInter interactions are
  generally not resolved to specific transcript pairs, so isoform-specific
  predictions are hypotheses rather than observations.
* **Strong class imbalance** (98.4 % / 1.6 %) and uneven species and biotype
  coverage. Absence of an interaction is not evidence of absence.
* **Evidence classes overlap** and originate from assays differing in
  directness and reliability; treat them as provenance, not ground-truth labels.
* **Annotation is incomplete and uneven**; some GO terms attach to the gene
  rather than to a specific transcript.
* **Inherits RNAInter v4.0's temporal coverage**; interactions reported later
  are absent.
* **No train/test division.** None is provided; splits are left to the user.
  They should separate shared entities, highly similar transcripts and Rfam
  families — random interaction-level splits alone give optimistic estimates.

## License and citation

Released under an open license, subject to the terms of the source resources
(RNAInter, RNAcentral, miRBase, circBase, RefSeq, Ensembl, GTEx, Rfam, GO).

Nicolini M., Ferreyra S., Paccanaro A., Robinson P.N., Casiraghi E.,
Valentini G. *RNAInterSeq: RNA interaction sequence-based database.*
Manuscript in preparation.
