"""Entity tables for the transcript-multiplicity / similarity analysis, keyed by the
pipeline entity key (symbol, biotype, NCBI taxid, identifiers) = key1/key2 of the release.

The tables previously used by 02_entity_figures.py (data/entities_union_{mrna,ncrna}.pkl)
were keyed by (symbol, biotype, species *name*) with the identifier component dropped:
`make_entity_key` received the Raw_ID string instead of a list, so `raw_ids` was always ().
Two consequences: entities that differ only by identifier were merged, and dog -- listed
by RNAInter both as "Canis familiaris" and "Canis lupus familiaris", taxid 9615 -- was
counted twice (8,798 entities in the mRNA-associated subset). That is why the manuscript
quoted 176,444 entities for a subset of a dataset that has 174,248 in total.

Per-entity pairwise-identity statistics are reused from the previous similarity run for
every entity whose (length-filtered) transcript set is unchanged, and computed pair-parallel
only for the rest, exactly as 00_ncrna_similarity.py does. Same measure: global
Needleman-Wunsch identity, match=1, mismatch=0, no gap penalty, normalised by the longer
sequence.

Outputs: data/entities_bykey_{mrna,ncrna}.pkl, data/df_step3_bykey_{mrna,ncrna}_with_similarity.pkl
"""
import pandas as pd, numpy as np, time, sys, os, glob, multiprocessing as mp
from itertools import combinations
from Bio.Align import PairwiseAligner

t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)
MAX_LEN = {"mrna": 100_000, "ncrna": 1_000_000}   # as in mRNAInter_stats_union.ipynb / 00_ncrna_similarity.py
NCBI_NAME = {"9615": "Canis lupus familiaris"}     # canonical NCBI Taxonomy name for the double-named taxid

aligner = PairwiseAligner()
aligner.mode = "global"; aligner.match_score = 1; aligner.mismatch_score = 0
aligner.open_gap_score = 0; aligner.extend_gap_score = 0
SEQS = []
def one_pair(t):
    e, i, j = t
    a = SEQS[e][i]; b = SEQS[e][j]
    return aligner.score(a, b) / max(len(a), len(b))

def nvalid(x): return [s for s in x if isinstance(s, str) and len(s) > 0] if isinstance(x, list) else []

def build_entity_table(sub):
    rows = []
    for side in [1, 2]:
        rows.append(pd.DataFrame({
            "entity_key": sub[f"key{side}"],
            "symbol": sub[f"Interactor{side}.Symbol"], "category": sub[f"Category{side}"],
            "taxid": sub[f"Species{side}_TaxID"].astype(str), "species_name": sub[f"Species{side}"],
            "seq_list": sub[f"ent{side}_sequences"]}))
    e = pd.concat(rows, ignore_index=True)
    e["key_str"] = e.entity_key.astype(str)
    e = e.drop_duplicates("key_str").reset_index(drop=True)
    e["species"] = [NCBI_NAME.get(t, n) for t, n in zip(e.taxid, e.species_name)]
    e["valid"] = e.seq_list.apply(nvalid)
    e["n_transcripts"] = e.valid.apply(len)
    e["mean_len"] = e.valid.apply(lambda v: float(np.mean([len(s) for s in v])) if v else np.nan)
    return e

def main():
    mp.set_start_method("fork", force=True)
    global SEQS
    PATH = sorted(glob.glob("data/interactions_with_sequences_RNA_RNA_union_*.pkl"))[-1]
    log("loading", PATH)
    df = pd.read_pickle(PATH)
    is_m = df.Category1.eq("mRNA") | df.Category2.eq("mRNA")
    # previous per-entity results, reused where the transcript set is unchanged: the notebook
    # tables and, when present, the by-key tables of an earlier run of this script
    for tag, sub, prev_fs in [("mrna", df[is_m], ["data/df_step3_union_with_similarity.pkl", "data/df_step3_bykey_mrna_with_similarity.pkl"]),
                              ("ncrna", df[~is_m], ["data/df_step3_ncrna_with_similarity.pkl", "data/df_step3_bykey_ncrna_with_similarity.pkl"])]:
        log(f"===== {tag}: {len(sub)} interactions")
        e = build_entity_table(sub)
        log(" entities by pipeline key:", len(e), "| taxids:", e.taxid.nunique(), "| species names:", e.species.nunique())
        e.drop(columns=["valid"]).to_pickle(f"data/entities_bykey_{tag}.pkl")
        # ---- similarity: reuse where the filtered transcript set is unchanged ----
        prev_stats = {}
        for prev_f in prev_fs:
            if not os.path.exists(prev_f): continue
            prev = pd.read_pickle(prev_f)
            for sl, npairs, ms, md, mn, mx, vs in zip(prev.seq_list_filtered, prev.n_pairs, prev.mean_sim,
                                                      prev.median_sim, prev.min_sim, prev.max_sim, prev.var_sim):
                sl = sl if isinstance(sl, list) else eval(sl)
                prev_stats[tuple(sorted(sl))] = dict(n_pairs=int(npairs), mean_sim=float(ms), median_sim=float(md),
                                                     min_sim=float(mn), max_sim=float(mx), var_sim=float(vs))
            log(f" previous results from {prev_f}: {len(prev_stats)} transcript sets"); del prev
        e["seq_list_filtered"] = e.valid.apply(lambda v: [s for s in v if len(s) <= MAX_LEN[tag]])
        e["n_transcripts_filtered"] = e.seq_list_filtered.apply(len)
        st = e[e.n_transcripts_filtered >= 2].copy().reset_index(drop=True)
        log(" multi-transcript entities:", len(st))
        keys = [tuple(sorted(l)) for l in st.seq_list_filtered]
        reuse = [prev_stats.get(k) for k in keys]
        todo = [i for i, r in enumerate(reuse) if r is None]
        log(f" reused: {len(st)-len(todo)} | to compute: {len(todo)}")
        if todo:
            SEQS = [st.seq_list_filtered.iloc[i] for i in todo]
            tasks = [(k, a, b) for k, L in enumerate(SEQS) for a, b in combinations(range(len(L)), 2)]
            tasks.sort(key=lambda t: -len(SEQS[t[0]][t[1]]) * len(SEQS[t[0]][t[2]]))
            cells = sum(len(SEQS[k][a]) * len(SEQS[k][b]) for k, a, b in tasks)
            log(f" pairs: {len(tasks)} | cells: {cells:.3e}")
            sims = {k: [] for k in range(len(SEQS))}
            with mp.Pool(mp.cpu_count()) as pool:
                for n, (t, v) in enumerate(zip(tasks, pool.imap(one_pair, tasks, chunksize=1)), 1):
                    sims[t[0]].append(v)
                    if n % 50000 == 0: log(f"   {n}/{len(tasks)}")
            for k, i in enumerate(todo):
                a = np.array(sims[k])
                reuse[i] = dict(n_pairs=len(a), mean_sim=a.mean(), median_sim=np.median(a),
                                min_sim=a.min(), max_sim=a.max(), var_sim=a.var())
            log(" alignments done")
        res = pd.DataFrame(reuse)
        st = pd.concat([st.drop(columns=["valid", "seq_list"]), res], axis=1)
        st.to_pickle(f"data/df_step3_bykey_{tag}_with_similarity.pkl")
        log(f" saved data/df_step3_bykey_{tag}_with_similarity.pkl | mean {st.mean_sim.mean():.4f} median {st.mean_sim.median():.4f} "
            f"min {st.mean_sim.min():.4f} max {st.mean_sim.max():.4f} var mean {st.var_sim.mean():.4f} var median {st.var_sim.median():.4f}")
        del e, st, prev_stats
    log("DONE")

if __name__ == "__main__":
    main()
