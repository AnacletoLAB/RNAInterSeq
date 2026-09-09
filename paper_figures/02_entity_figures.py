"""Entity-level panels: transcript multiplicity and pairwise identity.
Frequency panels report relative frequencies (EC/GV); boxplots and the scatter
keep their original scale."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib as mpl
import seaborn as sns, pandas as pd, numpy as np, time, os

t0=time.time()
def log(*a): print(f"[{time.time()-t0:7.1f}s]", *a, flush=True)
MD="data/mrna_stats_union/"; ND="data/ncrna_stats_union/"
mpl.rcParams.update({"figure.figsize":(7,5),"figure.dpi":300,"axes.labelsize":18,
    "axes.titlesize":20,"xtick.labelsize":16,"ytick.labelsize":16,
    "axes.spines.top":False,"axes.spines.right":False,"font.family":"sans-serif"})
sns.set_style("white")
PASTEL_COLOR="#A6CEE3"; PASTEL_SIM="#B2DF8A"; PASTEL_BOX="#FDBF6F"
bin_labels=["1","2","3–5","6–10",">10"]
def b(n): return "1" if n==1 else "2" if n==2 else "3–5" if 3<=n<=5 else "6–10" if 6<=n<=10 else ">10"

for tag, OUT, ents_f, sim_f, cnt_names in [
    ("mRNA",  MD, "data/entities_union_mrna.pkl",  "data/df_step3_union_with_similarity.pkl",
     ["transcripts_counts.pdf","transcripts-count-m.pdf"]),
    ("ncRNA", ND, "data/entities_union_ncrna.pkl", "data/df_step3_ncrna_with_similarity.pkl",
     ["transcripts_counts.pdf","transcripts-count-nc.pdf"]),
]:
    log(f"===== {tag}")
    e=pd.read_pickle(ents_f)
    e["valid"]=e["seq_list"].apply(lambda x:[s for s in x if isinstance(s,str)] if isinstance(x,list) else [])
    e["n"]=e["valid"].apply(len)
    counts=e["n"].apply(b).value_counts().reindex(bin_labels, fill_value=0)
    frac=counts/counts.sum()
    log(" entities:", len(e), "| bins:", counts.to_dict())
    for name in cnt_names:
        plt.figure(figsize=(7,5), dpi=300)
        sns.barplot(x=list(frac.index), y=frac.values, hue=list(frac.index),
                    palette=[PASTEL_COLOR]*len(frac), legend=False)
        plt.xlabel("Transcripts per RNA entity"); plt.ylabel("Fraction of entities")
        plt.tight_layout(); plt.savefig(OUT+name, bbox_inches="tight"); plt.close()
    log(" transcript-count panels done (proportion)")
    del e

    st=pd.read_pickle(sim_f)
    st["bin"]=st["n_transcripts_filtered"].apply(b)
    log(" multi-transcript entities:", len(st))

    plt.figure(figsize=(7,5), dpi=300)
    sns.histplot(st["mean_sim"].dropna(), bins=30, color=PASTEL_SIM, stat="probability")
    plt.xlabel("Mean pairwise transcript identity"); plt.ylabel("Fraction of RNA entities")
    plt.tight_layout(); plt.savefig(OUT+"transcript_pairwise_id.pdf", bbox_inches="tight"); plt.close()
    log(" transcript_pairwise_id done (proportion)")

    bo=["2","3–5","6–10",">10"]
    plt.figure(figsize=(8,5), dpi=300)
    sns.boxplot(x="bin", y="mean_sim", data=st, order=bo, color=PASTEL_BOX)
    plt.xlabel("Transcript multiplicity bin"); plt.ylabel("Mean pairwise identity")
    plt.tight_layout(); plt.savefig(OUT+"transcript_multiplicity_bin.pdf", bbox_inches="tight"); plt.close()
    log(" transcript_multiplicity_bin done")

    st["log_n_pairs"]=np.log10(st["n_pairs"]+1)
    plt.figure(figsize=(7,5), dpi=300)
    sns.scatterplot(x="log_n_pairs", y="mean_sim", hue="bin", hue_order=bo,
                    palette="Set2", data=st, alpha=0.7, s=12, edgecolor=None)
    plt.xlabel("log10(number of transcript pairs + 1)"); plt.ylabel("Mean pairwise identity")
    plt.legend(title="Transcript bin", fontsize=9, title_fontsize=10)
    plt.tight_layout(); plt.savefig(OUT+"transcript_id_vs_pairs.pdf", bbox_inches="tight"); plt.close()
    log(" transcript_id_vs_pairs done")
    del st
log("ALL ENTITY FIGURES DONE")
