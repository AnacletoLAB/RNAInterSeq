"""Script form of 3_Interactions_map.ipynb: project the sequence-resolved RNA
entities back onto both partners of every curated interaction and keep the
interactions whose two partners are both sequence-resolved.

Differences with respect to the notebook:
  * the mRNA entity table is the rebuilt one written by
    2c_refetch_refseq_failures.py (the notebook reads the pre-union file);
  * the transcript accessions of mRNA-entities (ent*_refseq, the RefSeq
    accessions in the same order as the RefSeq sequences of ent*_sequences, and
    ent*_ensembl_ids) are RETAINED in the release; the notebook dropped them,
    so the released mRNA sequences carried no accession-level provenance;
  * the output file name carries the run date.
Everything else -- the entity key, the merge, the filter, the column set --
is the notebook's.

Usage: python 3_interactions_map.py [mrna_entities.pkl] [out.pkl]
"""
import sys, time
import pandas as pd

DATE = "18-09-2026"
INTER = "data/interactions_processed_noSeq.pkl"
NC = "data/nc_rna_entities_processed.pkl"
MR = sys.argv[1] if len(sys.argv) > 1 else f"data/m_rna_entities_processed_union_{DATE}.pkl"
OUT = sys.argv[2] if len(sys.argv) > 2 else f"data/interactions_with_sequences_RNA_RNA_union_{DATE}.pkl"

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

df_filt = pd.read_pickle(INTER); log("curated interactions:", len(df_filt))
rna_nc = pd.read_pickle(NC); rna_mr = pd.read_pickle(MR)
log("ncRNA entities:", len(rna_nc), "| mRNA entities:", len(rna_mr), "from", MR)

def normalize_raw_ids(x):
    return tuple(x) if isinstance(x, list) else x

df_filt["Raw_ID1_key"] = df_filt["Raw_ID1_list"].apply(normalize_raw_ids)
df_filt["Raw_ID2_key"] = df_filt["Raw_ID2_list"].apply(normalize_raw_ids)
df_filt["key1"] = list(zip(df_filt["Interactor1.Symbol"], df_filt["Category1"], df_filt["Species1_TaxID"], df_filt["Raw_ID1_key"]))
df_filt["key2"] = list(zip(df_filt["Interactor2.Symbol"], df_filt["Category2"], df_filt["Species2_TaxID"], df_filt["Raw_ID2_key"]))

def prepare_entities(df):
    df = df.copy()
    df["Raw_key"] = df["RNA_IDs"].apply(normalize_raw_ids)
    df["entity_key"] = list(zip(df["Symbol"], df["Category"], df["Species"], df["Raw_key"]))
    return df

entities_all = pd.concat([prepare_entities(rna_nc), prepare_entities(rna_mr)], ignore_index=True)
log("entity entries:", len(entities_all))
entities_all = entities_all.drop_duplicates("entity_key")
log("unique entity keys:", len(entities_all))

entities1 = entities_all.add_prefix("ent1_").rename(columns={"ent1_entity_key": "key1"})
entities2 = entities_all.add_prefix("ent2_").rename(columns={"ent2_entity_key": "key2"})
df_merged = df_filt.merge(entities1, on="key1", how="left").merge(entities2, on="key2", how="left")
log("merged:", df_merged.shape)

has_seq1 = df_merged["ent1_sequences"].apply(lambda x: isinstance(x, list) and len(x) > 0)
has_seq2 = df_merged["ent2_sequences"].apply(lambda x: isinstance(x, list) and len(x) > 0)
log("lacking sequence: Interactor1", int((~has_seq1).sum()), "| Interactor2", int((~has_seq2).sum()),
    "| both", int((~has_seq1 & ~has_seq2).sum()))
df_filtered = df_merged[has_seq1 & has_seq2].copy()
log("interactions before filtering:", len(df_merged), "| after:", len(df_filtered))

cols_to_drop = [
    "Category1_std", "Category2_std",
    "ent1_Category", "ent2_Category",
    "ent1_RNA_IDs", "ent2_RNA_IDs",
    "ent1_Raw_ID_key", "ent2_Raw_ID_key",
    "ent1_Species", "ent2_Species",
    "ent1_Symbol", "ent2_Symbol",
    "ent1_GO_MF_count", "ent2_GO_MF_count",
    "ent1_GO_CC_count", "ent2_GO_CC_count",
    "ent1_GO_BP_count", "ent2_GO_BP_count",
    "ent1_refseq_seq", "ent2_refseq_seq",        # redundant with ent*_sequences
    "ent1_ensembl_seq", "ent2_ensembl_seq",      # idem
    "Raw_ID1_list", "Raw_ID2_list",
    "Raw_ID1_key", "Raw_ID2_key",
    # ent*_refseq and ent*_ensembl_ids are kept (transcript accessions of mRNA-entities)
]
df_filtered = df_filtered.drop(columns=[c for c in cols_to_drop if c in df_filtered.columns])
df_filtered = df_filtered.rename(columns={"ent1_GO": "ent1_GO_RNAcentral", "ent2_GO": "ent2_GO_RNAcentral",
                                          "ent1_refseq": "ent1_refseq_accessions", "ent2_refseq": "ent2_refseq_accessions",
                                          "ent1_ensembl_ids": "ent1_ensembl_ids", "ent2_ensembl_ids": "ent2_ensembl_ids"})
# put the accession columns next to the sequences
order = []
for c in df_filtered.columns:
    if c in ("ent1_refseq_accessions", "ent1_ensembl_ids", "ent2_refseq_accessions", "ent2_ensembl_ids"):
        continue
    order.append(c)
    if c == "ent1_sequences": order += ["ent1_refseq_accessions", "ent1_ensembl_ids"]
    if c == "ent2_sequences": order += ["ent2_refseq_accessions", "ent2_ensembl_ids"]
df_filtered = df_filtered[order]
log("final shape:", df_filtered.shape); log("columns:", list(df_filtered.columns))

is_m = df_filtered["Category1"].eq("mRNA") | df_filtered["Category2"].eq("mRNA")
log("mRNA-associated:", int(is_m.sum()), "| ncRNA-ncRNA:", int((~is_m).sum()))
keys = pd.concat([df_filtered["key1"], df_filtered["key2"]]).drop_duplicates()
log("unique RNA entities retained:", len(keys),
    "| mRNA:", sum(1 for k in keys if k[1] == "mRNA"), "| ncRNA:", sum(1 for k in keys if k[1] != "mRNA"))
df_filtered.to_pickle(OUT)
log("saved", OUT)
