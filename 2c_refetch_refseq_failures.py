"""Re-fetch the RefSeq accessions that the original efetch run never downloaded
and rebuild the mRNA entity table.

The March 2026 run of the fetch cell of 2_sequence_map_mRNA.ipynb lost 29 batch
requests (50 accessions each, 1,450 in total: "Failed batch:" in the notebook
output); the September 2026 union run only fetched the accessions that were new
with respect to that run, so the 1,450 were never retried. As a result 78
mRNA-entities carrying only failed accessions had no sequence at all and 491
resolved entities lacked 1,352 isoforms. The accessions are valid and are still
served by NCBI.

This script changes nothing but the fetch:
  * the RefSeq accession lists are recomputed from data/gene2refseq exactly as
    the notebook does (NCBI Gene step, then gene symbol step merged by union)
    and checked against the lists stored in the processed table;
  * the missing accessions are downloaded into a NEW cache file (the 3.8 GB
    main cache is left untouched);
  * refseq_seq / sequences are rebuilt for every entity, Ensembl sequences and
    GO / MyGene.info annotations are copied from the processed table, and a
    single MyGene.info request annotates the entities that were previously
    unresolved (the notebook fetched GO only for resolved entities).
The result is written to a new dated file; 3_interactions_map.py reads it.
"""
import json, os, sys, time
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

DATE = "18-09-2026"
TO_BE = "data/m_rna_entities_to_be_processed.pkl"
UNION = "data/m_rna_entities_processed_union.pkl"
CACHES = ["data/refseq_sequences_cache.json", "data/refseq_sequences_cache_union_delta.json"]
NEW_CACHE = f"data/refseq_sequences_cache_refetch_{DATE}.json"
OUT = f"data/m_rna_entities_processed_union_{DATE}.pkl"

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH_SIZE, MAX_WORKERS, RETRY, DELAY = 50, 3, 5, 0.35   # as in the notebook

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

# ---------------------------------------------------------------- 1. accession lists
ent = pd.read_pickle(TO_BE)
log("to-be-processed entities:", len(ent))
g2r = pd.read_csv("data/gene2refseq", sep="\t", comment="#", header=None, usecols=[0, 1, 2, 3, 15],
                  names=["taxid", "GeneID", "status", "rna_acc", "Symbol"], dtype=str, low_memory=False)
g2r = g2r[(g2r["status"] != "SUPPRESSED") & (g2r["rna_acc"] != "-")
          & (g2r["rna_acc"].str.startswith(("NM_", "XM_"), na=False))].copy()
log("gene2refseq mRNA rows:", len(g2r))
refseq_dict = {k: sorted(set(v)) for k, v in g2r.groupby(["taxid", "GeneID"])["rna_acc"]}
symbol_dict = {k: sorted(set(v)) for k, v in g2r.groupby(["taxid", "Symbol"])["rna_acc"]}
del g2r

def get_refseq_transcripts(bucket, species):            # cell 6
    if not isinstance(bucket, dict):
        return None
    tr = []
    for geneid in bucket.get("NCBI", []):
        tr.extend(refseq_dict.get((str(species), str(geneid)), []))
    return sorted(set(tr)) if tr else None

def get_refseq_by_symbol(species, symbol):              # cell 11
    if not symbol or not isinstance(symbol, str):
        return set()
    return set(symbol_dict.get((str(species), symbol), []))

ent["refseq"] = [get_refseq_transcripts(b, s) for b, s in zip(ent["ID_buckets"], ent["Species"])]
log("NCBI Gene step: entities with >=1 accession:", ent["refseq"].notna().sum())
stats = dict(matched_symbol=0, bucket_filled_empty=0, bucket_extended=0, bucket_no_new_transcripts=0)
merged = []
for species, symbol, prev in zip(ent["Species"], ent["Symbol"], ent["refseq"]):
    A = set(prev) if isinstance(prev, list) else set()
    N = get_refseq_by_symbol(species, symbol)
    if N: stats["matched_symbol"] += 1
    if not A:
        if N: stats["bucket_filled_empty"] += 1
        m = N
    else:
        if N:
            if N <= A: stats["bucket_no_new_transcripts"] += 1
            else: stats["bucket_extended"] += 1
        m = A | N
    merged.append(sorted(m) if m else None)
ent["refseq"] = merged
log("symbol step:", stats, "| entities with >=1 accession:", ent["refseq"].notna().sum())

# ---------------------------------------------------------------- 2. check against the processed table
proc = pd.read_pickle(UNION)
log("processed (union) table rows:", len(proc))
key = lambda df: list(zip(df["Symbol"].astype(str), df["Species"].astype(str), df["Raw_ID_key"].astype(str)))
ent["k"] = key(ent); proc["k"] = key(proc)
assert ent["k"].is_unique and proc["k"].is_unique
pk = dict(zip(proc["k"], proc["refseq"]))
common = ent[ent["k"].isin(pk)]
mism = sum(1 for k, r in zip(common["k"], common["refseq"]) if (pk[k] or []) != (r or []))
log("entities in common:", len(common), "| accession lists differing from the processed table:", mism)
assert mism == 0, "recomputed accession lists do not match the processed table"

# ---------------------------------------------------------------- 3. fetch what is missing
cache = {}
for f in CACHES:
    cache.update(json.load(open(f))); log("cache", f, "->", len(cache), "sequences")
if os.path.exists(NEW_CACHE):
    cache.update(json.load(open(NEW_CACHE))); log("existing", NEW_CACHE, "->", len(cache))
all_acc = sorted(set(a for L in ent["refseq"].dropna() for a in L))
to_fetch = [a for a in all_acc if a not in cache]
log("distinct accessions:", len(all_acc), "| to fetch:", len(to_fetch))

def parse_fasta(text):
    seqs, sid, seq = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if sid: seqs[sid] = "".join(seq)
            sid = line.split()[0][1:]; seq = []
        else:
            seq.append(line.strip())
    if sid: seqs[sid] = "".join(seq)
    return seqs

def fetch_batch(batch):
    params = {"db": "nuccore", "id": ",".join(batch), "rettype": "fasta", "retmode": "text"}
    s = requests.Session()
    for attempt in range(RETRY):
        try:
            r = s.get(EFETCH_URL, params=params, timeout=120)
            if r.status_code == 200:
                return parse_fasta(r.text)
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(2 ** attempt); continue
            log("  HTTP", r.status_code, "for", batch[:3]); time.sleep(2 ** attempt)
        except requests.exceptions.RequestException as e:
            log("  request error:", e); time.sleep(2 ** attempt)
    return None                                            # a failed batch is reported, never swallowed

new = {}
if to_fetch:
    batches = [to_fetch[i:i + BATCH_SIZE] for i in range(0, len(to_fetch), BATCH_SIZE)]
    failed = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = {ex.submit(fetch_batch, b): b for b in batches}
        for fu in as_completed(futs):
            res = fu.result()
            if res is None: failed.append(futs[fu])
            else: new.update(res)
            time.sleep(DELAY)
    # second pass, one accession at a time, for the batches that failed as a whole
    for b in failed:
        for a in b:
            res = fetch_batch([a])
            if res: new.update(res)
            time.sleep(DELAY)
    json.dump(new, open(NEW_CACHE, "w"))
    cache.update(new)
    log("fetched:", len(new), "->", NEW_CACHE, "| still missing:", sum(a not in cache for a in to_fetch))
missing = [a for a in all_acc if a not in cache]
log("accessions without a sequence after the re-fetch:", len(missing), missing[:10])

# ---------------------------------------------------------------- 4. rebuild the table
def extract_ensembl_ids(bucket):                          # cell 9
    if not isinstance(bucket, dict): return []
    return [x for x in bucket.get("Ensembl", []) if isinstance(x, str) and x.startswith("ENS")]
ens_seq = {k: v for k, v in zip(proc["k"], proc["ensembl_seq"]) if isinstance(v, list) and v}
ent["ensembl_ids"] = ent["ID_buckets"].apply(extract_ensembl_ids)
ent["ensembl_seq"] = [ens_seq.get(k) for k in ent["k"]]
log("entities with Ensembl sequences (copied):", sum(1 for k in ent["k"] if k in ens_seq), "(processed table:", len(ens_seq), ")")

def map_sequences(ids):                                    # cell 13
    if not ids: return None
    seqs = [cache.get(i) for i in ids if i in cache]
    return seqs if seqs else None
ent["refseq_seq"] = ent["refseq"].apply(map_sequences)

def merge_seqs(r, e):                                      # cell 14
    seqs = list(r) if r else []
    if e: seqs.extend(e)
    return seqs if seqs else None
ent["sequences"] = [merge_seqs(r, e) for r, e in zip(ent["refseq_seq"], ent["ensembl_seq"])]
n_acc = ent["refseq"].apply(lambda x: len(x) if isinstance(x, list) else 0)
n_seq = ent["refseq_seq"].apply(lambda x: len(x) if isinstance(x, list) else 0)
log("entities with accessions but fewer RefSeq sequences:", int(((n_acc > n_seq) & (n_acc > 0)).sum()),
    "| missing accession slots:", int((n_acc - n_seq)[n_acc > n_seq].sum()))
res = ent[ent["sequences"].apply(lambda x: x is not None and len(x) > 0)].copy()
log("sequence-resolved mRNA-entities:", len(res), "of", len(ent),
    "| previously:", len(proc), "| new:", sum(1 for k in res["k"] if k not in pk))

# ---------------------------------------------------------------- 5. GO / MyGene.info annotations
go_cols = ["GO_raw", "GO_BP_count", "GO_MF_count", "GO_CC_count", "mygene_raw"]
have = res["k"].isin(pk)
for c in go_cols:                                          # copied by key; object columns hold dicts/lists
    lut = dict(zip(proc["k"], proc[c]))
    res[c] = [lut.get(k) for k in res["k"]]
newrows = res[~have]
log("entities needing a MyGene.info query:", len(newrows))

def normalize_go(block):                                   # cell 16
    if not block: return []
    if isinstance(block, dict): block = [block]
    return [{"id": x.get("id"), "term": x.get("term"), "evidence": x.get("evidence"),
             "qualifier": x.get("qualifier"), "pubmed": x.get("pubmed")}
            for x in block if x.get("id") and x.get("term")]
go_cache, raw_cache = {}, {}
def fetch_go_batch(pairs, source="NCBI"):
    if not pairs: return
    ids = [str(g) for g, _ in pairs]
    data = {"ids": ",".join(ids), "fields": "go,taxid,symbol,name,type_of_gene,genomic_pos,summary"}
    if source == "Ensembl": data["scope"] = "ensembl.gene"
    r = requests.post("https://mygene.info/v3/gene", data=data,
                      headers={"Content-Type": "application/x-www-form-urlencoded"}, timeout=120)
    r.raise_for_status()
    tax_map = {str(g): str(t) for g, t in pairs}
    for rec in r.json():
        gid = str(rec.get("query")); tx = str(rec.get("taxid")) if rec.get("taxid") else None
        if gid not in tax_map or tx != tax_map[gid]: continue
        go = rec.get("go") or {}
        go_cache[(gid, tax_map[gid])] = {"BP": normalize_go(go.get("BP")), "MF": normalize_go(go.get("MF")),
                                         "CC": normalize_go(go.get("CC"))}
        raw_cache[(gid, tax_map[gid])] = rec
def id_pairs(row):
    out = []; taxid = str(int(row["Species"])); b = row["ID_buckets"]
    if not isinstance(b, dict): return out
    for src in ["NCBI", "Ensembl"]:
        ids = b.get(src, []); ids = [ids] if isinstance(ids, (int, str)) else ids
        out += [(str(g), taxid) for g in ids if pd.notna(g)]
    return out
if len(newrows):
    pairs = sorted({p for _, row in newrows.iterrows() for p in id_pairs(row)})
    fetch_go_batch([p for p in pairs if p[0].isdigit()], "NCBI")
    fetch_go_batch([p for p in pairs if not p[0].isdigit()], "Ensembl")
    log("MyGene.info records:", len(raw_cache), "of", len(pairs), "id/taxid pairs")
    def map_go(row):
        out = {"BP": [], "MF": [], "CC": []}
        for k in id_pairs(row):
            if k in go_cache:
                for c in out: out[c].extend(go_cache[k][c])
        for c in out:
            seen = set(); out[c] = [x for x in out[c] if x["id"] not in seen and not seen.add(x["id"])]
        return out
    def map_raw(row):
        for k in id_pairs(row):
            if k in raw_cache: return raw_cache[k]
        return None
    go_new = {k: map_go(r) for k, (_, r) in zip(newrows["k"], newrows.iterrows())}
    raw_new = {k: map_raw(r) for k, (_, r) in zip(newrows["k"], newrows.iterrows())}
    res["GO_raw"] = [go_new.get(k, v) for k, v in zip(res["k"], res["GO_raw"])]
    res["mygene_raw"] = [raw_new.get(k, v) for k, v in zip(res["k"], res["mygene_raw"])]
    for c in ("BP", "MF", "CC"):
        res[f"GO_{c}_count"] = [float(len(go_new[k][c])) if k in go_new else v for k, v in zip(res["k"], res[f"GO_{c}_count"])]
    log("annotated new entities with GO terms:", sum(1 for d in go_new.values() if any(d.values())), "of", len(go_new))
for c in ("GO_BP_count", "GO_MF_count", "GO_CC_count"):
    res[c] = res[c].astype(float)

# ---------------------------------------------------------------- 6. save, same columns as the processed table
res = res[[c for c in proc.columns if c != "k"]]
res.to_pickle(OUT)
nseq = res["sequences"].apply(len).sum(); ndist = len(set(s for L in res["sequences"] for s in L))
log("saved", OUT, "| rows", len(res), "| sequences", nseq, "| distinct", ndist,
    "| GO BP/MF/CC sums", int(res.GO_BP_count.sum()), int(res.GO_MF_count.sum()), int(res.GO_CC_count.sum()))
