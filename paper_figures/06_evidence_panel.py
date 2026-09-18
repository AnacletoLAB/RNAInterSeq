"""experimental_evidence_type.pdf for both subsets, from the method categories of
data/paper_numbers_<date>.json (05_paper_numbers.py), which reproduce the
hand-curated method -> category map of mRNAInter_stats.ipynb cells 24-26.
Until this script existed the panel was the notebook output and was not
regenerated with the release; the plot itself is the notebook's (percentage of
interactions per evidence category, counts next to the bars)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib as mpl
import seaborn as sns, pandas as pd, glob, json, os

sns.set_style("white")
mpl.rcParams.update({"figure.figsize": (7, 5), "figure.dpi": 300, "axes.labelsize": 18, "axes.titlesize": 20,
                     "xtick.labelsize": 14, "ytick.labelsize": 14, "axes.spines.top": False, "axes.spines.right": False,
                     "font.family": "sans-serif"})
N = json.load(open(sorted(glob.glob("data/paper_numbers_*.json"))[-1]))
print("release:", N["release"])
for tag, key, out in [("mRNA", "evidence_mrna", "data/mrna_stats_union/"), ("ncRNA", "evidence_ncrna", "data/ncrna_stats_union/")]:
    cs = pd.DataFrame(N[key]["method_categories"], columns=["Category", "Count", "Percentage"]).sort_values("Count", ascending=True)
    plt.figure(figsize=(10, 6))
    sns.barplot(x="Percentage", y="Category", data=cs, hue="Category", palette="viridis", legend=False)
    for i, count in enumerate(cs["Count"]):
        plt.text(cs["Percentage"].iloc[i] + 0.5, i, f"{count}", va="center")
    plt.title("Distribution of Experimental Evidence by Category", fontsize=14)
    plt.xlabel("Percentage of interactions"); plt.ylabel("")
    plt.xlim(0, cs["Percentage"].max() + 5)
    plt.tight_layout(); os.makedirs(out, exist_ok=True)
    plt.savefig(out + "experimental_evidence_type.pdf", bbox_inches="tight"); plt.close()
    print(tag, "->", out + "experimental_evidence_type.pdf", cs.to_dict("records"))
