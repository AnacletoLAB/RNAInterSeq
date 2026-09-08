"""
Step 3: turn the matched GTEx expression data into a coverage report and a
set of summary figures.

The matches table is long (one row per entity-locus pair: a handful of
entities, mostly miRNAs with paralogous precursor loci, resolve to more
than one GENCODE gene), so every entity-level count below is computed by
grouping on entity_key first rather than counting rows directly.

Tissue-specificity is summarized with Yanai's tau (0 = uniformly expressed
across all GTEx tissues, 1 = expressed in essentially one tissue only):
    tau = sum_i(1 - x_i/x_max) / (n_tissues - 1)
computed only over genes with at least one tissue showing nonzero expression.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DATA = ROOT / "data"
FIGS = ROOT / "figures"

ENTITIES_PKL = DATA / "ncrna_entities_human.pkl"
MATCHES_PKL = DATA / "gtex_gene_matches.pkl"
GENE_EXPR_PKL = DATA / "gtex_gene_expression.pkl"
TX_EXPR_PKL = DATA / "gtex_transcript_expression.pkl"
LOCUS_TX_MATCH_PKL = DATA / "ncrna_locus_transcript_matches.pkl"
REPORT_TXT = DATA / "gtex_coverage_report.txt"
LOCUS_SUMMARY_PKL = DATA / "ncrna_entity_locus_summary.pkl"
STATS_JSON = DATA / "gtex_stats_for_report.json"


def tau(values):
    x = np.asarray(values, dtype=float)
    xmax = x.max()
    if xmax <= 0:
        return np.nan
    return float(np.sum(1 - x / xmax) / (len(x) - 1))


def main():
    lines, stats = [], {}

    def report(msg=""):
        print(msg)
        lines.append(msg)

    ents = pd.read_pickle(ENTITIES_PKL)
    matches = pd.read_pickle(MATCHES_PKL)
    gene_expr = pd.read_pickle(GENE_EXPR_PKL)
    tx_expr = pd.read_pickle(TX_EXPR_PKL)
    locus_tx = pd.read_pickle(LOCUS_TX_MATCH_PKL)
    n_tissues = gene_expr["tissueSiteDetailId"].nunique()
    matched_loci = matches[matches["matched"]].copy()
    locus_confirmed = locus_tx.groupby(["entity_key", "gencode_id"])["sequence_confirmed"].any()

    report("=" * 70)
    report("GTEx tissue-expression coverage report")
    report("RNAInterSEQ ncRNA-ncRNA dataset (interactions_with_sequences_RNA_RNA_"
           "ncrna_ncrna_19-05-2026), GTEx v10/v39 primary + v8/v26 fallback")
    report("=" * 70)

    # ---- coverage funnel (entity-level) -----------------------------------
    per_entity = matches.groupby("entity_key").agg(
        category=("category", "first"), matched=("matched", "any"),
        n_loci=("gencode_id", lambda s: s.notna().sum()))
    n_total = len(per_entity)
    n_matched = int(per_entity["matched"].sum())
    n_multi_locus = int((per_entity["n_loci"] > 1).sum())

    gene_expr_by_gene = gene_expr.groupby("gencodeId")["median"].max()
    expressed_gencode_ids = set(gene_expr_by_gene[gene_expr_by_gene > 0].index)
    entity_expressed = matched_loci.groupby("entity_key")["gencode_id"].apply(
        lambda s: bool(set(s) & expressed_gencode_ids))
    n_expressed = int(entity_expressed.sum())

    tx_per_gene = tx_expr.groupby("gencodeId")["transcriptId"].nunique()
    n_single_tx = int((tx_per_gene == 1).sum())
    n_multi_tx = int((tx_per_gene > 1).sum())

    report("\n--- Coverage funnel ---")
    report(f"{n_total:>7,}  unique human ncRNA entities in the 19-05-2026 dataset")
    report(f"{n_matched:>7,}  resolved to >=1 GENCODE gene known to GTEx "
           f"({n_matched/n_total:.1%})")
    report(f"{n_multi_locus:>7,}  of those resolved to multiple genomic loci "
           f"(paralogous precursors, mostly miRNA)")
    report(f"{n_expressed:>7,}  detected (median TPM > 0) in >=1 of {n_tissues} GTEx "
           f"tissues ({n_expressed/max(n_matched,1):.1%} of matched)")
    report(f"{n_single_tx:>7,}  matched genes have exactly one GENCODE transcript")
    report(f"{n_multi_tx:>7,}  matched genes have multiple annotated transcripts "
           f"(isoform-level expression available)")

    stats["funnel"] = dict(n_total=n_total, n_matched=n_matched,
                            n_multi_locus=n_multi_locus, n_expressed=n_expressed,
                            n_single_tx=n_single_tx, n_multi_tx=n_multi_tx,
                            n_tissues=int(n_tissues))

    # ---- match rate by category --------------------------------------------
    report("\n--- Match rate by RNAInter category (entity-level) ---")
    rate = per_entity.groupby("category")["matched"].agg(n_total="count", n_matched="sum")
    rate["match_rate"] = rate["n_matched"] / rate["n_total"]
    rate = rate.sort_values("n_total", ascending=False)
    report(rate.to_string())
    stats["match_rate_by_category"] = rate.reset_index().to_dict(orient="records")

    # ---- which identifier system resolved each locus ----------------------
    report("\n--- Resolution method (entity-locus pairs) ---")
    method_counts = matched_loci["match_method"].value_counts()
    report(method_counts.to_string())
    stats["match_method_counts"] = method_counts.to_dict()

    report("\n--- GENCODE version used (entity-locus pairs) ---")
    version_counts = matched_loci["gencode_version"].value_counts()
    report(version_counts.to_string())
    n_v26_only_entities = int(
        matched_loci.groupby("entity_key")["gencode_version"]
        .apply(lambda s: set(s) == {"v26"}).sum()
    )
    report(f"\n{n_v26_only_entities:,} entities matched ONLY at the v26/gtex_v8 fallback "
           f"(no v39/gtex_v10 record exists for them -- GENCODE pruned these loci out "
           f"of the current annotation, e.g. several miRNA host genes).")
    stats["gencode_version_counts"] = version_counts.to_dict()
    stats["n_v26_only_entities"] = n_v26_only_entities

    # ---- sequence-confirmed transcripts ------------------------------------
    # Matching to a GENCODE gene (above) says "this entity's symbol/ID points
    # at this locus" -- it does not say which of that gene's transcripts is
    # actually the molecule RNAInterSEQ has on record. RNAInterSEQ's sequence
    # for an entity is usually sourced from miRBase/NCBI RefSeq/RNAcentral,
    # not Ensembl, so there's no shared ID to join on; 05_filter_transcripts_
    # by_sequence.py instead fetched every candidate GENCODE transcript's
    # cDNA sequence from Ensembl and tested it for containment against the
    # entity's own known sequence(s).
    n_loci_total = len(matched_loci)
    n_loci_confirmed = int(locus_confirmed.reindex(
        pd.MultiIndex.from_frame(matched_loci[["entity_key", "gencode_id"]]), fill_value=False
    ).sum())
    report(f"\n--- Sequence-confirmed transcripts ---")
    report(f"{n_loci_confirmed:,}/{n_loci_total:,} matched entity-locus pairs "
           f"({n_loci_confirmed/n_loci_total:.1%}) have >=1 GENCODE transcript whose "
           f"Ensembl-fetched sequence matches (exact or containment) a sequence "
           f"RNAInterSEQ already has on record for that entity.")
    stats["n_loci_confirmed"] = n_loci_confirmed
    stats["n_loci_total_for_confirmation"] = n_loci_total

    confirm_by_cat = (
        locus_tx.groupby(["entity_key", "gencode_id", "category"])["sequence_confirmed"]
        .any().reset_index()
        .groupby("category")["sequence_confirmed"].agg(n="count", confirmed="sum")
    )
    confirm_by_cat["rate"] = confirm_by_cat["confirmed"] / confirm_by_cat["n"]
    confirm_by_cat = confirm_by_cat.sort_values("n", ascending=False)
    report("\nSequence-confirmation rate by category (entity-locus level):")
    report(confirm_by_cat.to_string())
    stats["confirm_by_category"] = confirm_by_cat.reset_index().to_dict(orient="records")

    tx_per_gene_total = locus_tx.groupby("gencode_id")["transcript_id"].nunique()
    tx_per_gene_confirmed = (locus_tx[locus_tx["sequence_confirmed"]]
                              .groupby("gencode_id")["transcript_id"].nunique())
    narrowed = (tx_per_gene_confirmed / tx_per_gene_total.reindex(tx_per_gene_confirmed.index)).dropna()
    report(f"\nFor genes with >=1 confirmed transcript, the confirmed set is on average "
           f"{narrowed.mean():.1%} of that gene's total annotated transcripts "
           f"(median {narrowed.median():.1%}) -- sequence confirmation narrows "
           f"multi-isoform genes down to the specific isoform(s) RNAInterSEQ has.")
    stats["transcript_narrowing_mean"] = float(narrowed.mean())
    stats["transcript_narrowing_median"] = float(narrowed.median())

    # ---- category vs GTEx biotype agreement -------------------------------
    primary_locus = matched_loci.drop_duplicates("entity_key")
    agree = primary_locus.apply(
        lambda r: (r["category"].lower() in str(r["gtex_gene_type"]).lower()) or
                  (str(r["gtex_gene_type"]).lower() in r["category"].lower()),
        axis=1,
    )
    report(f"\n--- RNAInter category vs GTEx/GENCODE biotype agreement "
           f"(primary locus per entity) ---")
    report(f"{agree.mean():.1%} of matched entities have a GTEx biotype consistent "
           f"with their RNAInter category label.")
    stats["category_biotype_agreement"] = float(agree.mean())

    # ---- tissue specificity (tau) at gene level ----------------------------
    gene_tau = gene_expr.groupby("gencodeId")["median"].apply(tau).dropna()
    report(f"\n--- Tissue specificity (Yanai's tau, gene level, n={len(gene_tau):,}) ---")
    report(gene_tau.describe().to_string())
    stats["tau_describe"] = gene_tau.describe().to_dict()

    top_specific = (
        gene_expr.loc[gene_expr.groupby("gencodeId")["median"].idxmax()]
        .set_index("gencodeId")[["geneSymbol", "tissueSiteDetailId", "median"]]
        .join(gene_tau.rename("tau"))
        .sort_values("tau", ascending=False)
        .head(15)
    )
    report("\nTop 15 most tissue-specific matched ncRNAs (by tau) and their peak tissue:")
    report(top_specific.to_string())
    stats["top_specific"] = top_specific.reset_index().to_dict(orient="records")

    # ---- save entity-locus summary table -----------------------------------
    peak = gene_expr.loc[gene_expr.groupby("gencodeId")["median"].idxmax(),
                          ["gencodeId", "tissueSiteDetailId", "median"]]
    peak.columns = ["gencode_id", "peak_tissue", "peak_median_tpm"]
    summary = matched_loci.merge(
        gene_tau.rename("tau"), left_on="gencode_id", right_index=True, how="left"
    ).merge(
        tx_per_gene.rename("n_gtex_transcripts"), left_on="gencode_id", right_index=True, how="left"
    ).merge(peak, on="gencode_id", how="left")
    summary = summary.merge(
        locus_confirmed.rename("sequence_confirmed"),
        left_on=["entity_key", "gencode_id"], right_index=True, how="left")
    summary["sequence_confirmed"] = summary["sequence_confirmed"].fillna(False).infer_objects(copy=False)
    summary.to_pickle(LOCUS_SUMMARY_PKL)
    report(f"\nSaved per entity-locus GTEx summary table ({len(summary):,} rows) -> "
           f"{LOCUS_SUMMARY_PKL}")

    # ---- figures ------------------------------------------------------------
    FIGS.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 4))
    rate["match_rate"].plot.bar(ax=ax, color="#3b6ea5")
    ax.set_ylabel("Fraction matched to a GTEx gene")
    ax.set_title("GTEx match rate by RNAInter ncRNA category")
    plt.tight_layout()
    fig.savefig(FIGS / "Fig1_match_rate_by_category.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(gene_tau.values, bins=40, color="#5a9367")
    ax.set_xlabel("tau (tissue specificity)")
    ax.set_ylabel("number of ncRNA genes")
    ax.set_title("Tissue specificity of matched ncRNAs across GTEx tissues")
    plt.tight_layout()
    fig.savefig(FIGS / "Fig2_tau_distribution.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    capped = tx_per_gene.clip(upper=10)
    capped.value_counts().sort_index().plot.bar(ax=ax, color="#a55a5a")
    ax.set_xlabel("annotated transcripts per gene (10 = 10 or more)")
    ax.set_ylabel("number of ncRNA genes")
    ax.set_title("Transcript multiplicity of matched ncRNA genes")
    plt.tight_layout()
    fig.savefig(FIGS / "Fig3_transcripts_per_gene.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    method_counts.sort_values().plot.barh(ax=ax, color="#8a6bbf")
    ax.set_xlabel("entity-locus pairs")
    ax.set_title("Which identifier system resolved each GTEx match")
    plt.tight_layout()
    fig.savefig(FIGS / "Fig4_match_method_breakdown.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    version_counts.plot.bar(ax=ax, color="#c98a3b")
    ax.set_ylabel("entity-locus pairs")
    ax.set_xlabel("GENCODE version (paired GTEx dataset)")
    ax.set_title("v39/gtex_v10 vs v26/gtex_v8 fallback usage")
    plt.xticks(rotation=0)
    plt.tight_layout()
    fig.savefig(FIGS / "Fig5_gencode_version_usage.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    confirm_by_cat["rate"].plot.bar(ax=ax, color="#4a9b9b")
    ax.set_ylabel("Fraction with >=1 sequence-confirmed transcript")
    ax.set_title("Sequence-confirmation rate by RNAInter ncRNA category")
    plt.tight_layout()
    fig.savefig(FIGS / "Fig6_sequence_confirmation_rate.pdf")
    plt.close(fig)

    report(f"\nFigures saved to {FIGS}/")

    REPORT_TXT.write_text("\n".join(lines))
    STATS_JSON.write_text(pd.Series(stats).to_json())
    report(f"Full report written to {REPORT_TXT}")
    report(f"Machine-readable stats written to {STATS_JSON}")


if __name__ == "__main__":
    main()
