"""
Construct the subset of RNA-KG that does NOT overlap RNAInterSeq.

Rationale
---------
RNAInterSeq carries no train/test division. A model trained on RNAInterSeq
therefore needs an evaluation set drawn from elsewhere, and RNA-KG is the
obvious candidate -- but only after every interaction RNAInterSeq already
contains has been removed from it, otherwise the evaluation leaks.

Matching key
------------
RNAInterSeq stores interactions at the RNA-*entity* level and keeps every
transcript isoform recovered for an entity, whereas RNA-KG stores them at the
transcript level. There is no identifier shared by both resources, so the join
key is the nucleotide sequence itself, exactly as in the RNA-KG comparison
cells of mRNAInter_stats.ipynb / ncRNAInter_stats.ipynb. Each RNAInterSeq
interaction is expanded into all combinations of the transcripts of its two
partners before matching.

Two corrections relative to those notebook cells (see README):
  * sequences are normalised on BOTH sides (upper-case, T->U, strip);
  * interaction pairs are canonicalised as UNORDERED, since (A,B) and (B,A)
    denote the same interaction and the two resources do not agree on
    partner order.

Both are conservative: they can only find MORE overlap, i.e. remove more
records from the non-overlapping output.

Exactness
---------
No hashing is used. Every distinct RNAInterSeq sequence is assigned an integer
id; a canonical pair is encoded as a single int64 from the two ids. An RNA-KG
sequence absent from the RNAInterSeq vocabulary cannot take part in any
RNAInterSeq pair, so its record is non-overlapping by construction.
"""
import pandas as pd, numpy as np, glob, json, time, re, os
from itertools import product

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT  = os.path.join(ROOT, "data")
os.makedirs(OUT, exist_ok=True)

NONACGU = re.compile(r"[^AUCG]")
def norm(s):
    if not isinstance(s, str): return None
    return s.strip().upper().replace("T", "U")

# ----------------------------------------------------------------- RNAInterSeq
# Reference set: the released ncRNA-ncRNA dataset. RNA-KG's ncRNA export is
# ncRNA-ncRNA only, so this is the matching scope; verified to give exactly the
# same result as comparing against the full RNAInterSeq release (90,716
# overlapping / 21,404 kept either way).
RIS = os.path.join(ROOT, "..", "data",
                   "interactions_with_sequences_RNA_RNA_ncrna_ncrna_19-05-2026.pkl")
log("RNAInterSeq reference:", os.path.basename(RIS))
ris = pd.read_pickle(RIS)
log("  interactions:", len(ris))

seq2id, id_ct = {}, 0
# pickle memoises identical strings, so the same sequence object recurs across
# many interaction rows. Caching on object identity avoids re-normalising
# multi-kb sequences (TTN's transcripts are ~100 kb) once per row they occur in.
_by_obj, _keep = {}, []
def sid_cached(x):
    global id_ct
    k = id(x)
    j = _by_obj.get(k, -1)
    if j >= 0: return j
    s = norm(x)
    i = seq2id.get(s, -1)
    if i < 0:
        i = id_ct; seq2id[s] = i; id_ct += 1
    _by_obj[k] = i
    return i

def to_ids(col):
    out = []
    for l in col:
        if isinstance(l, list):
            _keep.append(l)          # keep refs alive so id() stays valid
            out.append([sid_cached(x) for x in l])
        else:
            out.append([])
    return out

ent1 = to_ids(ris["ent1_sequences"])
ent2 = to_ids(ris["ent2_sequences"])
log("  distinct normalised RNAInterSeq sequences:", len(seq2id))
SHIFT = 1
while (1 << SHIFT) <= len(seq2id): SHIFT += 1
log("  pair-key shift:", SHIFT)

keys = np.empty(sum(len(a)*len(b) for a, b in zip(ent1, ent2)), dtype=np.int64)
k = 0
for a, b in zip(ent1, ent2):
    for x, y in product(a, b):
        lo, hi = (x, y) if x <= y else (y, x)
        keys[k] = (lo << SHIFT) | hi; k += 1
log("  transcript-level ordered pairs expanded:", k)
del ent1, ent2, ris, _by_obj, _keep
ris_pairs = np.unique(keys); del keys
log("  distinct canonical (unordered) RNAInterSeq pairs:", len(ris_pairs))

summary = {"rnainterseq_release": os.path.basename(RIS),
           "rnainterseq_distinct_sequences": len(seq2id),
           "rnainterseq_canonical_pairs": int(len(ris_pairs))}

def in_ris(sx, sy):
    """Vectorised membership of canonical pairs in the RNAInterSeq pair set."""
    ix = np.fromiter((seq2id.get(s, -1) for s in sx), dtype=np.int64, count=len(sx))
    iy = np.fromiter((seq2id.get(s, -1) for s in sy), dtype=np.int64, count=len(sy))
    known = (ix >= 0) & (iy >= 0)
    lo = np.minimum(ix, iy); hi = np.maximum(ix, iy)
    kk = (lo << SHIFT) | hi
    pos = np.searchsorted(ris_pairs, kk)
    pos[pos >= len(ris_pairs)] = 0
    return known & (ris_pairs[pos] == kk), ix, iy

# ---------------------------------------------------------------------- RNA-KG
def run(tag, path, cx, cy, extra):
    log(f"--- {tag}: {os.path.basename(path)}")
    kg = pd.read_pickle(path)
    log("  rows:", len(kg))
    sx = [norm(s) for s in kg[cx]]
    sy = [norm(s) for s in kg[cy]]
    bad = sum(1 for s in sx[:200000] if s and NONACGU.search(s))
    log(f"  non-AUCG sequences in first 200k of {cx}: {bad}")

    ov, ix, iy = in_ris(sx, sy)
    n_ov = int(ov.sum())
    log(f"  overlapping records: {n_ov} ({100*n_ov/len(kg):.2f}%)")
    log(f"  NON-overlapping records: {len(kg)-n_ov} ({100*(1-n_ov/len(kg)):.2f}%)")

    out = kg.loc[~ov].copy()
    out.to_pickle(os.path.join(OUT, f"{tag}_pairs_non_overlapping.pkl"))
    log(f"  saved {tag}_pairs_non_overlapping.pkl ({len(out)} rows)")

    # ---- entity level: a sequence is non-overlapping if RNAInterSeq never has it
    ents = {}
    for s, i, idv, ty in [(sx, ix, kg[extra[0]], kg[extra[1]]),
                          (sy, iy, kg[extra[2]], kg[extra[3]])]:
        for seq, known, ident, typ in zip(s, i, idv, ty):
            if seq is None or known >= 0: continue          # present in RNAInterSeq
            if seq in ents: continue
            ents[seq] = (ident, typ if not isinstance(typ, list) else "|".join(map(str, typ)))
    allseq = {s for s in sx if s} | {s for s in sy if s}
    ent_df = pd.DataFrame([{"sequence": s, "rnakg_id": v[0], "rnakg_type": v[1],
                            "length": len(s)} for s, v in ents.items()])
    ent_df.to_pickle(os.path.join(OUT, f"{tag}_entities_non_overlapping.pkl"))
    log(f"  distinct RNA-KG sequences: {len(allseq)} | absent from RNAInterSeq: {len(ents)}")
    log(f"  saved {tag}_entities_non_overlapping.pkl")

    summary[tag] = {"source_file": os.path.basename(path),
                    "rnakg_records": int(len(kg)),
                    "overlapping_records": n_ov,
                    "non_overlapping_records": int(len(kg) - n_ov),
                    "rnakg_distinct_sequences": len(allseq),
                    "non_overlapping_sequences": len(ents)}
    del kg, sx, sy, out
    return

run("ncrna", os.path.join(ROOT, "..", "data", "rna-kg_ncRNA_new.pkl"),
    "RNA_sequence_x", "RNA_sequence_y", ("IdX", "Category_x", "IdY", "Category_y"))
# The RNA-KG mRNA export is out of scope for this deliverable.

json.dump(summary, open(os.path.join(OUT, "overlap_summary.json"), "w"), indent=2)
log("wrote data/overlap_summary.json")
log(json.dumps(summary, indent=2))
