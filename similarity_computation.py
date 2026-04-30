import pandas as pd
import numpy as np
from itertools import combinations
from joblib import Parallel, delayed
import multiprocessing
from Bio.Align import PairwiseAligner
import pickle
# -----------------------------
# load step 1 data
# -----------------------------
# from # save df_step3 to pickle file now
#with open("df_step3.pkl", "wb") as f:
#    pickle.dump(df_step3, f)
with open("df_step3.pkl", "rb") as f:
    df_step3 = pickle.load(f)

# -----------------------------
# global aligner (same config)
# -----------------------------
aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = 0
aligner.extend_gap_score = 0


def global_identity(seq1, seq2):
    score = aligner.score(seq1, seq2)
    max_len = max(len(seq1), len(seq2))
    return score / max_len


def transcript_similarity_stats(seq_list):
    seqs = [s for s in seq_list if isinstance(s, str) and len(s) > 0]

    if len(seqs) < 2:
        return {
            "n_pairs": 0,
            "mean_sim": np.nan,
            "median_sim": np.nan,
            "min_sim": np.nan,
            "max_sim": np.nan,
            "var_sim": np.nan
        }

    sims = [global_identity(s1, s2) for s1, s2 in combinations(seqs, 2)]
    sims = np.array(sims)

    return {
        "n_pairs": len(sims),
        "mean_sim": sims.mean(),
        "median_sim": np.median(sims),
        "min_sim": sims.min(),
        "max_sim": sims.max(),
        "var_sim": sims.var()
    }


# -----------------------------
# parallel execution (same logic)
# -----------------------------
n_cores = multiprocessing.cpu_count()

results = Parallel(n_jobs=n_cores)(
    delayed(transcript_similarity_stats)(seqs)
    for seqs in df_step3["seq_list_filtered"]
)

results_df = pd.DataFrame(results)

df_step3 = pd.concat(
    [df_step3.reset_index(drop=True), results_df],
    axis=1
)

# -----------------------------
# save output in pickle 
# -----------------------------
with open("df_step3_with_similarity.pkl", "wb") as f:
    pickle.dump(df_step3, f)

print("Similarity computation complete.")

