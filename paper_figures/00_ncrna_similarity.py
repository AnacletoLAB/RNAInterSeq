"""ncRNA transcript-similarity stats, pair-parallel (same measure as
similarity_computation.py: global Needleman-Wunsch identity, match=1,
mismatch=0, no gap penalty, normalised by the longer sequence)."""
import pandas as pd, numpy as np, time, multiprocessing as mp
from itertools import combinations
from Bio.Align import PairwiseAligner

t0=time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)
MAX_LEN=1_000_000
SEQS=[]
aligner=PairwiseAligner()
aligner.mode='global'; aligner.match_score=1; aligner.mismatch_score=0
aligner.open_gap_score=0; aligner.extend_gap_score=0

def one_pair(t):
    e,i,j=t
    a=SEQS[e][i]; b=SEQS[e][j]
    return aligner.score(a,b)/max(len(a),len(b))

def main():
    global SEQS
    mp.set_start_method("fork", force=True)
    e=pd.read_pickle("data/entities_union_ncrna.pkl")
    e["f"]=e["seq_list"].apply(lambda x:[s for s in x if isinstance(s,str) and len(s)<=MAX_LEN] if isinstance(x,list) else [])
    e["n"]=e["f"].apply(len)
    new=e[e["n"]>=2].copy().reset_index(drop=True)
    log("ncRNA entities with >=2 transcripts:", len(new))
    SEQS=list(new["f"])
    tasks=[(i,a,b) for i,L in enumerate(SEQS) for a,b in combinations(range(len(L)),2)]
    tasks.sort(key=lambda t: -len(SEQS[t[0]][t[1]])*len(SEQS[t[0]][t[2]]))
    cells=sum(len(SEQS[i][a])*len(SEQS[i][b]) for i,a,b in tasks)
    log(f"pairs: {len(tasks)} | cells: {cells:.3e}")
    sims={i:[] for i in range(len(SEQS))}
    with mp.Pool(mp.cpu_count()) as pool:
        for k,(t,v) in enumerate(zip(tasks, pool.imap(one_pair, tasks, chunksize=1)),1):
            sims[t[0]].append(v)
            if k%50000==0: log(f"  {k}/{len(tasks)}")
    log("alignments done")
    rows=[]
    for i,vals in sims.items():
        a=np.array(vals)
        rows.append({"n_pairs":len(a),"mean_sim":a.mean(),"median_sim":np.median(a),
                     "min_sim":a.min(),"max_sim":a.max(),"var_sim":a.var()})
    st=pd.concat([new.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    st=st.rename(columns={"n":"n_transcripts_filtered","f":"seq_list_filtered"})
    st.drop(columns=["seq_list"]).to_pickle("data/df_step3_ncrna_with_similarity.pkl")
    log("saved data/df_step3_ncrna_with_similarity.pkl")
    log(f"mean {st.mean_sim.mean():.6f} median {st.mean_sim.median():.6f} "
        f"min {st.mean_sim.min():.6f} max {st.mean_sim.max():.6f}")
    log(f"var mean {st.var_sim.mean():.6f} median {st.var_sim.median():.6f}")

if __name__=="__main__":
    main()
