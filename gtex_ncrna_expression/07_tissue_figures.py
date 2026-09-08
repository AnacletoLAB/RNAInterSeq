"""
Step 7: Figures on tissue-expression patterns in the human ncRNA-ncRNA interaction set.

Requires data/human_ncrna_interactions_gtex_annotated.pkl (run 06 first).

Figures produced:
  Fig6_step5_funnel.pdf           -- genes/transcripts before vs after seq filter
  Fig7_tissue_coexpression.pdf    -- per-tissue count of co-expressed interaction pairs
  Fig8_shared_tissue_histogram.pdf -- distribution of n_shared_expressed_tissues
  Fig9_category_pair_heatmap.pdf   -- fraction co-expressed by ncRNA category combination
  Fig10_tau_by_category.pdf        -- tissue-specificity (tau) violin per ncRNA category
  Fig11_peak_tissue_treemap.pdf    -- which peak tissues dominate the interaction set
"""
from pathlib import Path
import warnings

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as mticker

warnings.filterwarnings("ignore")

ROOT    = Path(__file__).parent
DATA    = ROOT / "data"
FIGURES = ROOT / "figures"
FIGURES.mkdir(exist_ok=True)

ANNOTATED_PKL  = DATA / "human_ncrna_interactions_gtex_annotated.pkl"
LOCUS_TX_PKL   = DATA / "ncrna_locus_transcript_matches.pkl"

# ── pastel palette ────────────────────────────────────────────────────────────
PASTEL = [
    "#A8C5DA",  # pastel blue
    "#F4A7A3",  # pastel rose
    "#B5D7A8",  # pastel green
    "#F9CEAA",  # pastel orange
    "#C8B4D8",  # pastel purple
    "#F9E4A3",  # pastel yellow
    "#B4D4D4",  # pastel teal
    "#F4C6D0",  # pastel pink
    "#C5D5B5",  # pastel sage
    "#D4C0A0",  # pastel sand
]

CATEGORY_ORDER = ["lncRNA", "miRNA", "pseudo", "snoRNA", "circRNA",
                  "snRNA", "scaRNA", "rRNA", "tRNA", "ncRNA"]

SYSTEM_ORDER = ["CNS", "Peripheral nervous", "Cardiovascular", "GI",
                "Reproductive", "Immune", "Muscle", "Other"]
CAT_COLORS = {c: PASTEL[i % len(PASTEL)] for i, c in enumerate(CATEGORY_ORDER)}

_TISSUE_LABELS = {
    "Adipose_Subcutaneous": "Adipose (subcut.)",
    "Adipose_Visceral_Omentum": "Adipose (visceral)",
    "Adrenal_Gland": "Adrenal gland",
    "Artery_Aorta": "Artery (aorta)",
    "Artery_Coronary": "Artery (coronary)",
    "Artery_Tibial": "Artery (tibial)",
    "Bladder": "Bladder",
    "Brain_Amygdala": "Brain – amygdala",
    "Brain_Anterior_cingulate_cortex_BA24": "Brain – ACC",
    "Brain_Caudate_basal_ganglia": "Brain – caudate",
    "Brain_Cerebellar_Hemisphere": "Brain – cerebellum (H)",
    "Brain_Cerebellum": "Brain – cerebellum",
    "Brain_Cortex": "Brain – cortex",
    "Brain_Frontal_Cortex_BA9": "Brain – frontal ctx",
    "Brain_Hippocampus": "Brain – hippocampus",
    "Brain_Hypothalamus": "Brain – hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia": "Brain – nucleus acc.",
    "Brain_Putamen_basal_ganglia": "Brain – putamen",
    "Brain_Spinal_cord_cervical_c-1": "Brain – spinal cord",
    "Brain_Substantia_nigra": "Brain – substantia nigra",
    "Breast_Mammary_Tissue": "Breast",
    "Cells_Cultured_fibroblasts": "Cells – fibroblasts",
    "Cells_EBV-transformed_lymphocytes": "Cells – lymphocytes",
    "Cervix_Ectocervix": "Cervix (ecto)",
    "Cervix_Endocervix": "Cervix (endo)",
    "Colon_Sigmoid": "Colon (sigmoid)",
    "Colon_Transverse": "Colon (transverse)",
    "Esophagus_Gastroesophageal_Junction": "Esophagus (GEJ)",
    "Esophagus_Mucosa": "Esophagus (mucosa)",
    "Esophagus_Muscularis": "Esophagus (musc.)",
    "Fallopian_Tube": "Fallopian tube",
    "Heart_Atrial_Appendage": "Heart (atrium)",
    "Heart_Left_Ventricle": "Heart (ventricle)",
    "Kidney_Cortex": "Kidney (cortex)",
    "Kidney_Medulla": "Kidney (medulla)",
    "Liver": "Liver",
    "Lung": "Lung",
    "Minor_Salivary_Gland": "Salivary gland",
    "Muscle_Skeletal": "Muscle",
    "Nerve_Tibial": "Nerve (tibial)",
    "Ovary": "Ovary",
    "Pancreas": "Pancreas",
    "Pituitary": "Pituitary",
    "Prostate": "Prostate",
    "Skin_Not_Sun_Exposed_Suprapubic": "Skin (unexposed)",
    "Skin_Sun_Exposed_Lower_leg": "Skin (exposed)",
    "Small_Intestine_Terminal_Ileum": "Small intestine",
    "Spleen": "Spleen",
    "Stomach": "Stomach",
    "Testis": "Testis",
    "Thyroid": "Thyroid",
    "Uterus": "Uterus",
    "Vagina": "Vagina",
    "Whole_Blood": "Whole blood",
}

def _short(tissue):
    return _TISSUE_LABELS.get(tissue, tissue.replace("_", " "))


def _save(fig, name, tight=True):
    path = FIGURES / name
    if tight:
        fig.savefig(path, bbox_inches="tight", dpi=150)
    else:
        fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  saved {path.name}")


# ── Fig 6: step 5 funnel (genes & transcripts before/after seq filter) ────────

def fig_step5_funnel(locus_tx):
    confirmed = locus_tx[locus_tx["sequence_confirmed"]]
    cats = [c for c in CATEGORY_ORDER if c in locus_tx["category"].unique()]

    genes_before = locus_tx.groupby("category")["gencode_id"].nunique().reindex(cats, fill_value=0)
    genes_after  = confirmed.groupby("category")["gencode_id"].nunique().reindex(cats, fill_value=0)
    tx_before    = locus_tx.groupby("category")["transcript_id"].nunique().reindex(cats, fill_value=0)
    tx_after     = confirmed.groupby("category")["transcript_id"].nunique().reindex(cats, fill_value=0)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    x = np.arange(len(cats))
    w = 0.38

    for ax, before, after, label in [
        (axes[0], genes_before, genes_after, "GENCODE genes"),
        (axes[1], tx_before,    tx_after,    "GENCODE transcripts"),
    ]:
        ax.bar(x - w/2, before.values, width=w, color="#A8C5DA", label="before filter", zorder=3)
        ax.bar(x + w/2, after.values,  width=w, color="#B5D7A8", label="after filter",  zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels(cats, rotation=35, ha="right", fontsize=9)
        ax.set_ylabel(f"Unique {label}", fontsize=9)
        ax.set_title(f"{label}: before vs after\nsequence-confirmed filter", fontsize=10)
        ax.legend(fontsize=8, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
        ax.grid(axis="y", color="#e0e0e0", zorder=0)

    fig.suptitle("Step 5 – sequence-based transcript filtering", fontsize=11, y=1.01)
    fig.tight_layout()
    _save(fig, "Fig6_step5_funnel.pdf")


# ── Fig 7: per-tissue co-expression count ─────────────────────────────────────

def fig_tissue_coexpression(ann):
    both_conf = ann[ann["both_seq_confirmed"]].copy()
    if both_conf.empty:
        print("  [skip Fig7] no both-confirmed rows")
        return

    tissue_counts = {}
    all_tissues = list(_TISSUE_LABELS.keys())
    for t in all_tissues:
        mask = both_conf["shared_expressed_tissues"].apply(lambda lst: t in lst)
        tissue_counts[t] = mask.sum()

    counts = pd.Series(tissue_counts).sort_values(ascending=False)

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(range(len(counts)), counts.values, color=PASTEL[0], zorder=3)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels([_short(t) for t in counts.index],
                       rotation=55, ha="right", fontsize=7)
    ax.set_ylabel("Interaction pairs\n(both partners detected)", fontsize=9)
    ax.set_title("Co-expressed ncRNA-ncRNA interaction pairs per GTEx tissue\n"
                 "(sequence-confirmed; median TPM > 0 in both partners)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e0e0e0", zorder=0)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    fig.tight_layout()
    _save(fig, "Fig7_tissue_coexpression.pdf")


# ── Fig 8: distribution of n_shared_expressed_tissues ─────────────────────────

def fig_shared_tissue_histogram(ann):
    both_conf = ann[ann["both_seq_confirmed"]]["n_shared_expressed_tissues"]
    if both_conf.empty:
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    bins = range(0, 56)
    ax.hist(both_conf, bins=bins, color=PASTEL[0], edgecolor="white", linewidth=0.4, zorder=3)
    ax.axvline(both_conf.median(), color="#888", linestyle="--", linewidth=1.2,
               label=f"Median = {both_conf.median():.0f}")
    ax.set_xlabel("Number of shared expressed tissues (out of 54)", fontsize=9)
    ax.set_ylabel("Interaction pairs", fontsize=9)
    ax.set_title("Tissue overlap between interacting ncRNA partners\n"
                 "(sequence-confirmed pairs only)", fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e0e0e0", zorder=0)
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    fig.tight_layout()
    _save(fig, "Fig8_shared_tissue_histogram.pdf")


# ── Fig 9: category-pair heatmap (fraction co-expressed) ─────────────────────

def fig_category_pair_heatmap(ann):
    both_m = ann[ann["both_gtex_matched"]].copy()
    if both_m.empty:
        return

    cats = [c for c in CATEGORY_ORDER
            if c in both_m["Category1"].values or c in both_m["Category2"].values]
    cat_set = set(cats)
    both_m = both_m[both_m["Category1"].isin(cat_set) & both_m["Category2"].isin(cat_set)]

    # Sort the category pair so (c1, c2) and (c2, c1) map to the same cell.
    # This gives a symmetric matrix without iterrows.
    both_m = both_m.copy()
    both_m["_cA"] = both_m[["Category1", "Category2"]].min(axis=1)
    both_m["_cB"] = both_m[["Category1", "Category2"]].max(axis=1)

    count_raw = pd.crosstab(both_m["_cA"], both_m["_cB"]).reindex(index=cats, columns=cats, fill_value=0)
    # Make fully symmetric: add upper and lower triangles, avoid doubling the diagonal.
    count_mat = count_raw.copy()
    for c1 in cats:
        for c2 in cats:
            v = int(count_raw.loc[c1, c2]) if c1 in count_raw.index and c2 in count_raw.columns else 0
            v2 = int(count_raw.loc[c2, c1]) if c2 in count_raw.index and c1 in count_raw.columns else 0
            count_mat.loc[c1, c2] = v + (v2 if c1 != c2 else 0)

    both_c = both_m[both_m["both_seq_confirmed"] & (both_m["n_shared_expressed_tissues"] > 0)]
    coex_raw = pd.crosstab(both_c["_cA"], both_c["_cB"]).reindex(index=cats, columns=cats, fill_value=0)
    coex_mat = coex_raw.copy()
    for c1 in cats:
        for c2 in cats:
            v = int(coex_raw.loc[c1, c2]) if c1 in coex_raw.index and c2 in coex_raw.columns else 0
            v2 = int(coex_raw.loc[c2, c1]) if c2 in coex_raw.index and c1 in coex_raw.columns else 0
            coex_mat.loc[c1, c2] = v + (v2 if c1 != c2 else 0)

    mat_n = count_mat.astype(float)
    frac = (coex_mat.astype(float) / mat_n.replace(0, np.nan)).fillna(np.nan)

    cmap = LinearSegmentedColormap.from_list(
        "pastel_heat", ["#f7f3ef", "#A8C5DA", "#3a7ca5"], N=256)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5),
                             gridspec_kw={"width_ratios": [1, 1]})

    for ax, data, title, fmt in [
        (axes[0], mat_n,  "Interaction count\n(both GTEx-matched)", ".0f"),
        (axes[1], frac,   "Fraction co-expressed\n(≥1 shared tissue, seq-confirmed)", ".2f"),
    ]:
        im = ax.imshow(data.values, cmap=cmap, aspect="auto",
                       vmin=0, vmax=(data.max().max() if not data.isnull().all().all() else 1))
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels(cats, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(cats)))
        ax.set_yticklabels(cats, fontsize=8)
        ax.set_title(title, fontsize=9, pad=8)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        for i in range(len(cats)):
            for j in range(len(cats)):
                v = data.iloc[i, j]
                if not np.isnan(v):
                    text = f"{v:{fmt}}" if v > 0 else ""
                    ax.text(j, i, text, ha="center", va="center",
                            fontsize=6.5, color="#444")

    fig.suptitle("ncRNA category-pair co-expression landscape", fontsize=11)
    fig.tight_layout()
    _save(fig, "Fig9_category_pair_heatmap.pdf")


# ── Fig 10: tau violin by ncRNA category ──────────────────────────────────────

def fig_tau_violin(ann):
    # pool unique entities (deduplicate by entity key to avoid counting the same
    # entity once per interaction row it participates in)
    rows = []
    for key_col, tau_col, cat_col in [
        ("key1", "ent1_tau", "Category1"),
        ("key2", "ent2_tau", "Category2"),
    ]:
        sub = ann[[key_col, tau_col, cat_col]].copy()
        sub.columns = ["entity_key", "tau", "category"]
        rows.append(sub)
    pool = (pd.concat(rows)
              .dropna(subset=["tau"])
              .drop_duplicates(subset="entity_key"))

    cats = [c for c in CATEGORY_ORDER if c in pool["category"].unique()]
    data = [pool[pool["category"] == c]["tau"].values for c in cats]
    colors = [CAT_COLORS[c] for c in cats]

    fig, ax = plt.subplots(figsize=(9, 4.5))
    parts = ax.violinplot(data, positions=range(len(cats)),
                          showmedians=True, showextrema=False)
    for i, (body, col) in enumerate(zip(parts["bodies"], colors)):
        body.set_facecolor(col)
        body.set_edgecolor("#aaaaaa")
        body.set_alpha(0.85)
    parts["cmedians"].set_color("#555")
    parts["cmedians"].set_linewidth(1.5)

    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, rotation=30, ha="right", fontsize=9)
    ax.set_ylabel("Yanai tissue-specificity (tau)", fontsize=9)
    ax.set_ylim(0, 1.05)
    ax.set_title("Tissue specificity of interacting ncRNAs by category\n"
                 "(tau = 0 ubiquitous, tau = 1 single-tissue; one point per entity per interaction)",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#e0e0e0", zorder=0)
    fig.tight_layout()
    _save(fig, "Fig10_tau_by_category.pdf")


# ── Fig 11: peak tissue distribution ─────────────────────────────────────────

def fig_peak_tissue(ann):
    # pool unique entities (deduplicate by entity key)
    rows = []
    for key_col, pt_col in [("key1", "ent1_peak_tissue"), ("key2", "ent2_peak_tissue")]:
        sub = ann[[key_col, pt_col]].copy()
        sub.columns = ["entity_key", "peak_tissue"]
        rows.append(sub)
    pool = (pd.concat(rows)
              .drop_duplicates(subset="entity_key")
              .dropna(subset=["peak_tissue"]))

    top_tissues = pool["peak_tissue"].value_counts().head(20)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.barh([_short(t) for t in top_tissues.index],
            top_tissues.values,
            color=PASTEL[0], edgecolor="white", linewidth=0.5, zorder=3)
    ax.invert_yaxis()
    ax.set_xlabel("Unique ncRNA entities\n(both interaction partners pooled)", fontsize=9)
    ax.set_title("Top 20 peak expression tissues for interacting ncRNAs", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#e0e0e0", zorder=0)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    fig.tight_layout()
    _save(fig, "Fig11_peak_tissue_distribution.pdf")


# ── Fig 12: tau correlation between interaction partners ──────────────────────

def fig_tau_scatter(ann):
    both = ann[ann["both_seq_confirmed"] &
               ann["ent1_tau"].notna() & ann["ent2_tau"].notna()].copy()
    if len(both) < 50:
        return

    cat_pairs = both["Category1"] + "–" + both["Category2"]
    common_pairs = cat_pairs.value_counts().head(5).index.tolist()

    fig, ax = plt.subplots(figsize=(6, 6))
    # background scatter (all pairs, light grey)
    ax.scatter(both["ent1_tau"], both["ent2_tau"],
               alpha=0.05, s=8, color="#cccccc", rasterized=True, zorder=1)

    for i, pair in enumerate(common_pairs):
        mask = cat_pairs == pair
        ax.scatter(both.loc[mask, "ent1_tau"], both.loc[mask, "ent2_tau"],
                   alpha=0.3, s=12, color=PASTEL[i % len(PASTEL)],
                   label=pair, rasterized=True, zorder=2)

    # diagonal reference
    ax.plot([0, 1], [0, 1], color="#aaaaaa", linewidth=0.8, linestyle="--", zorder=3)

    r = both[["ent1_tau", "ent2_tau"]].corr().iloc[0, 1]
    ax.text(0.05, 0.93, f"Pearson r = {r:.2f}", transform=ax.transAxes,
            fontsize=9, color="#444")
    ax.set_xlabel("Tau – partner 1", fontsize=9)
    ax.set_ylabel("Tau – partner 2", fontsize=9)
    ax.set_title("Tissue-specificity correlation between interacting partners\n"
                 "(sequence-confirmed pairs with GTEx data)", fontsize=10)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7.5, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    _save(fig, "Fig12_tau_partner_scatter.pdf")


# ── Fig 13: tissue system co-expression ──────────────────────────────────────

def fig_tissue_system_coexpression(ann):
    both_conf = ann[ann["both_seq_confirmed"]].copy()
    if both_conf.empty:
        print("  [skip Fig13] no both-confirmed rows")
        return

    systems = SYSTEM_ORDER

    # pairs per system
    sys_counts = {s: 0 for s in systems}
    for lst in both_conf["shared_tissue_systems"]:
        for s in lst:
            if s in sys_counts:
                sys_counts[s] += 1
    counts = pd.Series(sys_counts).reindex(systems, fill_value=0)

    # system-pair co-occurrence matrix (symmetric: pair shares system A AND B)
    cooc = pd.DataFrame(0, index=systems, columns=systems, dtype=int)
    for lst in both_conf["shared_tissue_systems"]:
        sl = [s for s in lst if s in systems]
        for s1 in sl:
            for s2 in sl:
                cooc.loc[s1, s2] += 1

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5),
                             gridspec_kw={"width_ratios": [1, 1.4]})

    # ── left: bar chart per system ────────────────────────────────────────
    ax = axes[0]
    sys_rev = systems[::-1]
    colors_rev = [PASTEL[i % len(PASTEL)] for i in range(len(systems))][::-1]
    ax.barh(sys_rev, counts[sys_rev].values, color=colors_rev,
            edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_xlabel("Co-expressed interaction pairs", fontsize=9)
    ax.set_title("Co-expressed pairs\nper biological system", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#e0e0e0", zorder=0)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    for i, (s, v) in enumerate(zip(sys_rev, counts[sys_rev].values)):
        ax.text(v + counts.max() * 0.01, i, f"{v:,}", va="center", fontsize=8, color="#444")

    # ── right: co-occurrence heatmap ──────────────────────────────────────
    ax = axes[1]
    cmap = LinearSegmentedColormap.from_list(
        "pastel_heat", ["#f7f3ef", "#A8C5DA", "#3a7ca5"], N=256)
    vmax = cooc.values.max() or 1
    im = ax.imshow(cooc.values, cmap=cmap, aspect="auto", vmin=0, vmax=vmax)
    ax.set_xticks(range(len(systems)))
    ax.set_xticklabels(systems, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(systems)))
    ax.set_yticklabels(systems, fontsize=8)
    ax.set_title("System-pair co-occurrence\n(both partners expressed in A and B)", fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Interaction pairs")
    for i in range(len(systems)):
        for j in range(len(systems)):
            v = int(cooc.iloc[i, j])
            if v > 0:
                ax.text(j, i, f"{v:,}", ha="center", va="center",
                        fontsize=6.5, color="#222" if v < vmax * 0.6 else "#fff")

    fig.suptitle("Biological system grouping of co-expressed ncRNA-ncRNA pairs\n"
                 "(sequence-confirmed; 54 GTEx tissues collapsed to 8 systems)",
                 fontsize=11)
    fig.tight_layout()
    _save(fig, "Fig13_tissue_system_coexpression.pdf")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading annotated interaction table …")
    ann      = pd.read_pickle(ANNOTATED_PKL)
    locus_tx = pd.read_pickle(LOCUS_TX_PKL)
    print(f"  {len(ann):,} human-human interaction rows")

    print("\nFig 6 – step 5 funnel …")
    fig_step5_funnel(locus_tx)

    print("Fig 7 – tissue co-expression bar chart …")
    fig_tissue_coexpression(ann)

    print("Fig 8 – shared-tissue histogram …")
    fig_shared_tissue_histogram(ann)

    print("Fig 9 – category-pair heatmap …")
    fig_category_pair_heatmap(ann)

    print("Fig 10 – tau violin by category …")
    fig_tau_violin(ann)

    print("Fig 11 – peak tissue distribution …")
    fig_peak_tissue(ann)

    print("Fig 12 – tau partner scatter …")
    fig_tau_scatter(ann)

    print("Fig 13 – tissue system co-expression …")
    fig_tissue_system_coexpression(ann)

    print("\nDone. All figures saved to", FIGURES)


if __name__ == "__main__":
    main()
