"""Assemble the Zenodo release folder from the current pipeline outputs.

  data/RNAInterSeq_interactions_full.pkl   the release table (pandas pickle, all columns)
  data/RNAInterSeq_human_ncRNA_GTEx.pkl     human ncRNA-ncRNA subset + GTEx annotation (pickle)
  data/RNAInterSeq_interactions.tsv.gz      one row per interaction: RNAInter fields, evidence,
                                            and the two entity identifiers
  data/RNAInterSeq_entities.tsv.gz          one row per RNA entity: key, biotype, taxonomy,
                                            identifiers, transcript accessions, GO, Rfam
  data/RNAInterSeq_transcripts.tsv.gz       one row per (entity, transcript): accession, length
  data/RNAInterSeq_transcripts.fasta.gz     the transcript sequences, header entity_id|index|accession
  data/RNAInterSeq_human_ncRNA_GTEx.tsv.gz  the GTEx annotation columns, flattened
  CHECKSUMS.json / CHECKSUMS.md5

The TSV/FASTA files carry the same information as the pickles in a normalised,
language-neutral form; entity_id (E + 6 digits) is assigned to the unique
entity keys of the release in sorted order and is the join key between the
three tables. List-valued fields are ';'-separated, dict-valued fields are JSON;
a missing gene symbol (NaN in the pickle) is written as an empty field.
"""
import glob, gzip, hashlib, json, os, re, sys, time
import pandas as pd

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)
ROOT = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(ROOT)
os.chdir(REPO)
SRC = sorted(glob.glob("data/interactions_with_sequences_RNA_RNA_union_*.pkl"))[-1]
GTEX = "gtex_ncrna_expression/data/human_ncrna_interactions_gtex_annotated.pkl"
D = os.path.join(ROOT, "data"); os.makedirs(D, exist_ok=True)

def nv(x): return [s for s in x if isinstance(s, str) and len(s) > 0] if isinstance(x, list) else []
import numpy as np
def lst(x):
    if isinstance(x, (np.ndarray, pd.Series)): x = list(x)
    return ";".join(str(v) for v in x if not (isinstance(v, float) and np.isnan(v))) if isinstance(x, (list, tuple, set)) and len(x) else ""
def js(x): return json.dumps(x, default=str) if isinstance(x, (dict, list)) else ""

df = pd.read_pickle(SRC); log("release:", SRC, df.shape); SHAPE = {"RNAInterSeq_interactions_full.pkl": df.shape}
# canonical, hashable form of the entity key: a missing symbol is NaN in the pickle (NaN != NaN,
# so such keys cannot be looked up by value) and becomes '' here
def ck(k): return tuple("" if (isinstance(x, float) and x != x) else x for x in k)
df["key1"] = df["key1"].apply(ck); df["key2"] = df["key2"].apply(ck)
import shutil
shutil.copyfile(SRC, os.path.join(D, "RNAInterSeq_interactions_full.pkl")); log("copied pickle")

# ---------------------------------------------------------------- entities
ent = {}
for side in ("1", "2"):
    cols = {f"ent{side}_ID_buckets": "ID_buckets", f"ent{side}_sequences": "sequences", f"ent{side}_typecheck": "typecheck",
            f"ent{side}_GO_RNAcentral": "GO_RNAcentral", f"ent{side}_Rfam": "Rfam", f"ent{side}_GO_raw": "GO_raw",
            f"ent{side}_mygene_raw": "mygene_raw", f"ent{side}_refseq_accessions": "refseq_accessions",
            f"ent{side}_ensembl_ids": "ensembl_ids", f"Species{side}": "species_name"}
    sub = df[[f"key{side}"] + list(cols)].rename(columns={f"key{side}": "key", **cols}).drop_duplicates("key")
    for r in sub.itertuples(index=False):
        if r.key not in ent: ent[r.key] = r
keys = sorted(ent, key=lambda k: (str(k[0]), str(k[1]), str(k[2]), str(k[3])))
eid = {k: f"E{i+1:06d}" for i, k in enumerate(keys)}
log("unique entities:", len(keys))

def go_ids(raw, dom):
    if not isinstance(raw, dict): return ""
    return ";".join(sorted({x["id"] for x in raw.get(dom, []) if isinstance(x, dict) and x.get("id")}))
rows = []
for k in keys:
    r = ent[k]
    rows.append(dict(entity_id=eid[k], key=repr(k), symbol=k[0], biotype=k[1], taxid=k[2], species_name=r.species_name,
                     rnainter_ids=";".join(k[3]) if isinstance(k[3], tuple) else str(k[3]),
                     id_buckets=js(r.ID_buckets), typecheck="" if pd.isna(r.typecheck) else str(r.typecheck),
                     n_transcripts=len(nv(r.sequences)),
                     refseq_accessions=lst(r.refseq_accessions), ensembl_gene_ids=lst(r.ensembl_ids),
                     GO_RNAcentral=lst(r.GO_RNAcentral),
                     GO_BP=go_ids(r.GO_raw, "BP"), GO_MF=go_ids(r.GO_raw, "MF"), GO_CC=go_ids(r.GO_raw, "CC"),
                     mygene_type_of_gene=(r.mygene_raw or {}).get("type_of_gene", "") if isinstance(r.mygene_raw, dict) else "",
                     Rfam=lst(r.Rfam)))
E = pd.DataFrame(rows)
E.to_csv(os.path.join(D, "RNAInterSeq_entities.tsv.gz"), sep="\t", index=False, compression="gzip")
log("entities.tsv.gz written:", E.shape)

# ---------------------------------------------------------------- transcripts (tsv + fasta)
ntx = 0
with gzip.open(os.path.join(D, "RNAInterSeq_transcripts.fasta.gz"), "wt") as fa, \
     gzip.open(os.path.join(D, "RNAInterSeq_transcripts.tsv.gz"), "wt") as tsv:
    tsv.write("entity_id\ttranscript_index\taccession\tsource\tlength\n")
    for k in keys:
        r = ent[k]; seqs = nv(r.sequences)
        acc = list(r.refseq_accessions) if isinstance(r.refseq_accessions, list) else []
        for i, s in enumerate(seqs):
            if k[1] == "mRNA":
                a, src = (acc[i], "RefSeq") if i < len(acc) else ("", "Ensembl")
            else:
                a, src = "", "ncRNA pipeline"
            fa.write(f">{eid[k]}|{i}|{a or 'NA'}\n")
            for j in range(0, len(s), 80): fa.write(s[j:j+80] + "\n")
            tsv.write(f"{eid[k]}\t{i}\t{a}\t{src}\t{len(s)}\n"); ntx += 1
log("transcripts written:", ntx)

# ---------------------------------------------------------------- interactions
I = pd.DataFrame({
    "RNAInterID": df.RNAInterID, "entity1_id": [eid[k] for k in df.key1], "entity2_id": [eid[k] for k in df.key2],
    "Interactor1.Symbol": df["Interactor1.Symbol"], "Category1": df.Category1, "Species1": df.Species1, "Species1_TaxID": df.Species1_TaxID, "Raw_ID1": df.Raw_ID1,
    "Interactor2.Symbol": df["Interactor2.Symbol"], "Category2": df.Category2, "Species2": df.Species2, "Species2_TaxID": df.Species2_TaxID, "Raw_ID2": df.Raw_ID2,
    "score": df.score, "strong": df.strong, "weak": df.weak, "predict": df.predict})
I.to_csv(os.path.join(D, "RNAInterSeq_interactions.tsv.gz"), sep="\t", index=False, compression="gzip")
log("interactions.tsv.gz written:", I.shape)
del df, ent, I

# ---------------------------------------------------------------- GTEx subset
g = pd.read_pickle(GTEX); log("GTEx table:", g.shape); SHAPE["RNAInterSeq_human_ncRNA_GTEx.pkl"] = g.shape
shutil.copyfile(GTEX, os.path.join(D, "RNAInterSeq_human_ncRNA_GTEx.pkl"))
G = pd.DataFrame({"RNAInterID": g.RNAInterID, "entity1_id": [eid[k] for k in g.key1], "entity2_id": [eid[k] for k in g.key2]})
for side in ("1", "2"):
    p = f"ent{side}_"
    G[p + "gtex_matched"] = g[p + "gtex_matched"]; G[p + "seq_confirmed"] = g[p + "seq_confirmed"]
    G[p + "confirmed_transcripts"] = g[p + "confirmed_transcripts"].apply(js)
    for c in ("tau", "peak_tissue", "peak_tpm", "n_expressed_tissues"): G[p + c] = g[p + c]
    G[p + "tissue_tpm"] = g[p + "tissue_tpm"].apply(js); G[p + "tissue_systems"] = g[p + "tissue_systems"].apply(lst)
for c in ("both_gtex_matched", "both_seq_confirmed", "n_shared_expressed_tissues", "n_shared_tissue_systems", "n_shared_systems"): G[c] = g[c]
for c in ("shared_expressed_tissues", "shared_tissue_systems", "shared_systems"): G[c] = g[c].apply(lst)
G.to_csv(os.path.join(D, "RNAInterSeq_human_ncRNA_GTEx.tsv.gz"), sep="\t", index=False, compression="gzip")
log("GTEx tsv written:", G.shape)

# ---------------------------------------------------------------- checksums
ck = {}
for f in sorted(os.listdir(D)):
    p = os.path.join(D, f); h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""): h.update(chunk)
    ck[f] = {"md5": h.hexdigest(), "bytes": os.path.getsize(p)}
for f, (r, c) in SHAPE.items(): ck[f].update(rows=int(r), cols=int(c))
json.dump(ck, open(os.path.join(ROOT, "CHECKSUMS.json"), "w"), indent=2)
with open(os.path.join(ROOT, "CHECKSUMS.md5"), "w") as fh:
    for f, v in ck.items(): fh.write(f"{v['md5']}  data/{f}\n")
log("checksums written"); print(json.dumps(ck, indent=2))
