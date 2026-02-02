# RNAInter_sequenceParser

**Description**
`RNAInter_sequenceParser` is a preprocessing toolkit for building sequence-level ncRNA–ncRNA RNAInter dataset. It standardizes RNA identifiers across NCBI, Ensembl, and RNAcentral, resolves many-to-many mappings, filters for non-coding RNA types, and outputs deduplicated, model-ready sequences. The pipeline is species-aware and designed for reproducible construction of large-scale interaction datasets.

**Features**

* Translate NCBI Gene IDs to Ensembl Gene/Transcript IDs
* Map transcripts to RNAcentral URS sequences
* Filter and select non-coding RNA types
* Deduplicate identical sequences for ML-ready datasets
* Species-aware processing

