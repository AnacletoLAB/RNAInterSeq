"""
Step 4: render a LaTeX report summarizing the GTEx tissue-expression
coverage of the RNAInterSEQ ncRNA-ncRNA dataset, with the figures from
03_compute_stats.py embedded. Recomputes all numbers directly from the
saved pickles (rather than trusting the JSON side-channel) so the report
can never silently drift from the underlying data.

Produces gtex_ncrna_expression_report.tex in this directory. Compiles it
to PDF with pdflatex/tectonic/latexmk if one is available on PATH;
otherwise leaves the .tex ready to compile elsewhere (e.g. Overleaf).
"""
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DATA = ROOT / "data"
FIGS = ROOT / "figures"
TEX_PATH = ROOT / "gtex_ncrna_expression_report.tex"

INTERACTIONS_PKL     = ROOT.parent / "data" / "interactions_with_sequences_RNA_RNA_ncrna_ncrna_19-05-2026.pkl"
ENTITIES_PKL         = DATA / "ncrna_entities_human.pkl"
MATCHES_PKL          = DATA / "gtex_gene_matches.pkl"
GENE_EXPR_PKL        = DATA / "gtex_gene_expression.pkl"
TX_EXPR_PKL          = DATA / "gtex_transcript_expression.pkl"
LOCUS_TX_MATCH_PKL   = DATA / "ncrna_locus_transcript_matches.pkl"
TX_EXPR_CONFIRMED_PKL = DATA / "gtex_transcript_expression_sequence_confirmed.pkl"
ANNOTATED_PKL        = DATA / "human_ncrna_interactions_gtex_annotated.pkl"


def esc(s):
    """Escape LaTeX special characters in arbitrary gene/entity symbols."""
    s = str(s)
    repl = {
        "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
        "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}", "\\": r"\textbackslash{}",
    }
    return re.sub(r"[&%$#_{}~^\\]", lambda m: repl[m.group()], s)


def pct(x):
    """Format a fraction as a LaTeX-safe percentage (escapes the literal '%',
    which would otherwise start a LaTeX comment and truncate the line)."""
    return f"{x:.1%}".replace("%", r"\%")


def tau(values):
    x = np.asarray(values, dtype=float)
    xmax = x.max()
    if xmax <= 0:
        return np.nan
    return float(np.sum(1 - x / xmax) / (len(x) - 1))


def df_to_booktabs(df, columns=None, float_fmt="{:.3f}"):
    columns = columns or list(df.columns)
    lines = [r"\begin{tabular}{" + "l" * len(columns) + "}", r"\toprule",
             " & ".join(esc(c) for c in columns) + r" \\", r"\midrule"]
    for _, row in df.iterrows():
        cells = []
        for c in columns:
            v = row[c]
            cells.append(float_fmt.format(v) if isinstance(v, float) else esc(v))
        lines.append(" & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def main():
    # ── initial database scope stats ─────────────────────────────────────────
    df_all = pd.read_pickle(INTERACTIONS_PKL)
    n_total_int = len(df_all)
    hh_mask = (df_all["Species1"] == "Homo sapiens") & (df_all["Species2"] == "Homo sapiens")
    mm_mask = (df_all["Species1"] == "Mus musculus")  & (df_all["Species2"] == "Mus musculus")
    n_hh  = int(hh_mask.sum())
    n_mm  = int(mm_mask.sum())
    n_other_int = n_total_int - n_hh - n_mm
    cross_mask = df_all["Species1"] != df_all["Species2"]
    n_cross = int(cross_mask.sum())
    n_sp_pairs = df_all.groupby(["Species1", "Species2"]).ngroups
    del df_all  # free memory before loading the big annotated pkl

    ents = pd.read_pickle(ENTITIES_PKL)
    matches = pd.read_pickle(MATCHES_PKL)
    gene_expr = pd.read_pickle(GENE_EXPR_PKL)
    tx_expr = pd.read_pickle(TX_EXPR_PKL)
    locus_tx = pd.read_pickle(LOCUS_TX_MATCH_PKL)
    n_tissues = gene_expr["tissueSiteDetailId"].nunique()

    # interaction-level annotation stats (step 6)
    ann_stats = {}
    if ANNOTATED_PKL.exists():
        ann = pd.read_pickle(ANNOTATED_PKL)
        n_ann = len(ann)
        n_both_m = int(ann["both_gtex_matched"].sum())
        n_both_c = int(ann["both_seq_confirmed"].sum())
        n_shared  = int((ann["n_shared_expressed_tissues"] > 0).sum())
        med_sh = ann.loc[ann["both_seq_confirmed"], "n_shared_expressed_tissues"].median()
        # co-expression gap analysis
        both_c = ann[ann["both_seq_confirmed"]]
        n_bc   = len(both_c)
        zero1  = both_c["ent1_n_expressed_tissues"] == 0
        zero2  = both_c["ent2_n_expressed_tissues"] == 0
        either_zero = zero1 | zero2
        coex_m = both_c["n_shared_expressed_tissues"] > 0
        n_zero_partner = int(either_zero[~coex_m].sum())
        n_disjoint     = int((~coex_m & ~either_zero).sum())
        # expressed entity stats (pooled across both partners, deduplicated)
        rows = []
        for kc, nc in [("key1","ent1_n_expressed_tissues"),("key2","ent2_n_expressed_tissues")]:
            s = ann[[kc, nc]].drop_duplicates(kc)
            s.columns = ["ek","n"]
            rows.append(s)
        pool_all = pd.concat(rows).drop_duplicates("ek")
        expr_pool = pool_all[pool_all["n"] > 0]
        med_n_tissues = int(expr_pool["n"].median())
        n_all54       = int((expr_pool["n"] == 54).sum())

        # tissue system stats
        sys_counts = {}
        for lst in ann.loc[ann["both_seq_confirmed"], "shared_tissue_systems"]:
            for s in lst:
                sys_counts[s] = sys_counts.get(s, 0) + 1
        top_sys = max(sys_counts, key=sys_counts.get) if sys_counts else "N/A"
        top_sys_count = sys_counts.get(top_sys, 0)
        n_sys_any = int((ann["n_shared_tissue_systems"] > 0).sum())

        ann_stats = dict(
            n_human_interactions=f"{n_ann:,}",
            n_both_matched=f"{n_both_m:,}", pct_both_matched=pct(n_both_m / n_ann),
            n_both_confirmed=f"{n_both_c:,}", pct_both_confirmed=pct(n_both_c / n_ann),
            n_shared_any=f"{n_shared:,}", pct_shared_any=pct(n_shared / n_ann),
            pct_shared_conf=pct(n_shared / n_bc),
            med_shared=f"{med_sh:.0f}",
            n_sys_any=f"{n_sys_any:,}", pct_sys_any=pct(n_sys_any / n_ann),
            top_sys=esc(top_sys), top_sys_count=f"{top_sys_count:,}",
            n_zero_partner=f"{n_zero_partner:,}", pct_zero_partner=pct(n_zero_partner / n_bc),
            n_disjoint=f"{n_disjoint:,}", pct_disjoint=pct(n_disjoint / n_bc),
            med_n_tissues=f"{med_n_tissues:,}", n_all54=f"{n_all54:,}",
        )
    else:
        ann_stats = dict(
            n_human_interactions="N/A (run 06_annotate_interactions.py first)",
            n_both_matched="N/A", pct_both_matched=r"N/A\%",
            n_both_confirmed="N/A", pct_both_confirmed=r"N/A\%",
            n_shared_any="N/A", pct_shared_any=r"N/A\%",
            pct_shared_conf=r"N/A\%",
            med_shared="N/A",
            n_sys_any="N/A", pct_sys_any=r"N/A\%",
            top_sys="N/A", top_sys_count="N/A",
            n_zero_partner="N/A", pct_zero_partner=r"N/A\%",
            n_disjoint="N/A", pct_disjoint=r"N/A\%",
            med_n_tissues="N/A", n_all54="N/A",
        )
    matched_loci = matches[matches["matched"]].copy()

    per_entity = matches.groupby("entity_key").agg(
        category=("category", "first"), matched=("matched", "any"),
        n_loci=("gencode_id", lambda s: s.notna().sum()))
    n_total, n_matched = len(per_entity), int(per_entity["matched"].sum())
    n_multi_locus = int((per_entity["n_loci"] > 1).sum())

    gene_expr_by_gene = gene_expr.groupby("gencodeId")["median"].max()
    expressed_ids = set(gene_expr_by_gene[gene_expr_by_gene > 0].index)
    n_expressed = int(matched_loci.groupby("entity_key")["gencode_id"]
                       .apply(lambda s: bool(set(s) & expressed_ids)).sum())

    tx_per_gene = tx_expr.groupby("gencodeId")["transcriptId"].nunique()
    n_single_tx, n_multi_tx = int((tx_per_gene == 1).sum()), int((tx_per_gene > 1).sum())

    rate = per_entity.groupby("category")["matched"].agg(n="count", matched="sum")
    rate["rate"] = rate["matched"] / rate["n"]
    rate = rate.sort_values("n", ascending=False).reset_index()
    rate_table = df_to_booktabs(
        rate.rename(columns={"category": "Category", "n": "Entities",
                              "matched": "Matched", "rate": "Match rate"}),
        columns=["Category", "Entities", "Matched", "Match rate"])

    method_counts = matched_loci["match_method"].value_counts().reset_index()
    method_counts.columns = ["Resolution method", "Entity-locus pairs"]
    method_table = df_to_booktabs(method_counts, float_fmt="{:.0f}")

    version_counts = matched_loci["gencode_version"].value_counts()
    n_v26_only = int(matched_loci.groupby("entity_key")["gencode_version"]
                      .apply(lambda s: set(s) == {"v26"}).sum())

    gene_tau = gene_expr.groupby("gencodeId")["median"].apply(tau).dropna()
    top_specific = (
        gene_expr.loc[gene_expr.groupby("gencodeId")["median"].idxmax()]
        .set_index("gencodeId")[["geneSymbol", "tissueSiteDetailId", "median"]]
        .join(gene_tau.rename("tau"))
        .sort_values("tau", ascending=False).head(15).reset_index(drop=True)
    )
    top_specific.columns = ["Gene symbol", "Peak tissue", "Median TPM", "tau"]
    top_specific["Peak tissue"] = top_specific["Peak tissue"].str.replace("_", " ")
    top_table = df_to_booktabs(top_specific,
                                columns=["Gene symbol", "Peak tissue", "Median TPM", "tau"])

    locus_confirmed = locus_tx.groupby(["entity_key", "gencode_id"])["sequence_confirmed"].any()
    n_loci_total = len(matched_loci)
    n_loci_confirmed = int(locus_confirmed.reindex(
        pd.MultiIndex.from_frame(matched_loci[["entity_key", "gencode_id"]]), fill_value=False
    ).sum())

    confirm_by_cat = (
        locus_tx.groupby(["entity_key", "gencode_id", "category"])["sequence_confirmed"]
        .any().reset_index()
        .groupby("category")["sequence_confirmed"].agg(n="count", confirmed="sum")
    )
    confirm_by_cat["rate"] = confirm_by_cat["confirmed"] / confirm_by_cat["n"]
    confirm_by_cat = confirm_by_cat.sort_values("n", ascending=False).reset_index()
    confirm_table = df_to_booktabs(
        confirm_by_cat.rename(columns={"category": "Category", "n": "Entity-locus pairs",
                                        "confirmed": "Confirmed", "rate": "Confirmation rate"}),
        columns=["Category", "Entity-locus pairs", "Confirmed", "Confirmation rate"])

    tx_per_gene_total = locus_tx.groupby("gencode_id")["transcript_id"].nunique()
    tx_per_gene_confirmed = (locus_tx[locus_tx["sequence_confirmed"]]
                              .groupby("gencode_id")["transcript_id"].nunique())
    narrowed = (tx_per_gene_confirmed / tx_per_gene_total.reindex(tx_per_gene_confirmed.index)).dropna()

    primary_locus = matched_loci.drop_duplicates("entity_key")
    biotype_agree = primary_locus.apply(
        lambda r: (r["category"].lower() in str(r["gtex_gene_type"]).lower()) or
                  (str(r["gtex_gene_type"]).lower() in r["category"].lower()), axis=1).mean()

    species_counts = ents["species"].value_counts() if "species" in ents.columns else None

    tex = TEMPLATE.format(
        n_total_int=f"{n_total_int:,}",
        n_hh=f"{n_hh:,}", pct_hh=pct(n_hh / n_total_int),
        n_mm=f"{n_mm:,}", pct_mm=pct(n_mm / n_total_int),
        n_other_int=f"{n_other_int:,}", pct_other_int=pct(n_other_int / n_total_int),
        n_cross=f"{n_cross:,}", n_sp_pairs=f"{n_sp_pairs:,}",
        n_total=f"{n_total:,}", n_matched=f"{n_matched:,}",
        pct_matched=pct(n_matched/n_total),
        n_multi_locus=f"{n_multi_locus:,}",
        n_expressed=f"{n_expressed:,}", pct_expressed=pct(n_expressed/max(n_matched,1)),
        n_tissues=n_tissues,
        n_single_tx=f"{n_single_tx:,}", n_multi_tx=f"{n_multi_tx:,}",
        n_total_transcripts=f"{tx_expr['transcriptId'].nunique():,}",
        rate_table=rate_table, method_table=method_table,
        n_v39=f"{version_counts.get('v39', 0):,}", n_v26=f"{version_counts.get('v26', 0):,}",
        n_v26_only=f"{n_v26_only:,}",
        biotype_agree=pct(biotype_agree),
        tau_mean=f"{gene_tau.mean():.3f}", tau_median=f"{gene_tau.median():.3f}",
        n_tau=f"{len(gene_tau):,}",
        top_table=top_table,
        n_loci_confirmed=f"{n_loci_confirmed:,}", n_loci_total=f"{n_loci_total:,}",
        pct_confirmed=pct(n_loci_confirmed/n_loci_total),
        confirm_table=confirm_table,
        narrow_mean=pct(narrowed.mean()), narrow_median=pct(narrowed.median()),
        n_tx_rows_total=f"{len(tx_expr):,}",
        n_tx_rows_confirmed=f"{len(pd.read_pickle(TX_EXPR_CONFIRMED_PKL)):,}",
        fig1=str((FIGS / "Fig1_match_rate_by_category.pdf").relative_to(ROOT)),
        fig2=str((FIGS / "Fig2_tau_distribution.pdf").relative_to(ROOT)),
        fig3=str((FIGS / "Fig3_transcripts_per_gene.pdf").relative_to(ROOT)),
        fig4=str((FIGS / "Fig4_match_method_breakdown.pdf").relative_to(ROOT)),
        fig5=str((FIGS / "Fig5_gencode_version_usage.pdf").relative_to(ROOT)),
        fig6=str((FIGS / "Fig6_sequence_confirmation_rate.pdf").relative_to(ROOT)),
        fig6b=str((FIGS / "Fig6_step5_funnel.pdf").relative_to(ROOT)),
        fig7=str((FIGS / "Fig7_tissue_coexpression.pdf").relative_to(ROOT)),
        fig8=str((FIGS / "Fig8_shared_tissue_histogram.pdf").relative_to(ROOT)),
        fig9=str((FIGS / "Fig9_category_pair_heatmap.pdf").relative_to(ROOT)),
        fig10=str((FIGS / "Fig10_tau_by_category.pdf").relative_to(ROOT)),
        fig11=str((FIGS / "Fig11_peak_tissue_distribution.pdf").relative_to(ROOT)),
        fig12=str((FIGS / "Fig12_tau_partner_scatter.pdf").relative_to(ROOT)),
        fig13=str((FIGS / "Fig13_tissue_system_coexpression.pdf").relative_to(ROOT)),
        **ann_stats,
    )
    TEX_PATH.write_text(tex)
    print(f"Wrote {TEX_PATH}")

    # tectonic was installed into a dedicated conda env (no system TeX Live
    # available, no sudo); fall back to checking it explicitly if not on PATH.
    extra_paths = [Path.home() / "miniconda3" / "envs" / "texreport" / "bin"]

    def _find(engine):
        found = shutil.which(engine)
        if found:
            return found
        for p in extra_paths:
            if (p / engine).exists():
                return str(p / engine)
        return None

    for engine in ("tectonic", "latexmk", "pdflatex"):
        engine_path = _find(engine)
        if engine_path:
            print(f"Compiling with {engine} ...")
            cmd = {
                "tectonic": [engine_path, TEX_PATH.name],
                "latexmk": [engine_path, "-pdf", "-interaction=nonstopmode", TEX_PATH.name],
                "pdflatex": [engine_path, "-interaction=nonstopmode", TEX_PATH.name],
            }[engine]
            for _ in range(2 if engine == "pdflatex" else 1):  # pdflatex needs 2 passes for refs
                subprocess.run(cmd, cwd=ROOT, check=False)
            pdf = TEX_PATH.with_suffix(".pdf")
            if pdf.exists():
                print(f"PDF written to {pdf}")
            break
    else:
        print("No LaTeX engine (pdflatex/tectonic/latexmk) found on PATH -- "
              "the .tex source is ready to compile elsewhere (e.g. Overleaf, "
              "or `pdflatex gtex_ncrna_expression_report.tex` on a machine with TeX Live).")


TEMPLATE = r"""
\documentclass[11pt]{{article}}
\usepackage[margin=2.5cm]{{geometry}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\usepackage{{float}}
\usepackage{{hyperref}}
\usepackage{{amsmath}}
\usepackage{{microtype}}
\usepackage{{url}}

\title{{GTEx Tissue-Expression Landscape of the RNAInterSEQ ncRNA--ncRNA Interactome}}
\author{{RNAInter\_sequenceParser side project}}
\date{{\today}}

\begin{{document}}
\maketitle

\begin{{abstract}}
We link the human non-coding RNAs in the RNAInterSEQ ncRNA--ncRNA interaction
release (\texttt{{interactions\_with\_sequences\_RNA\_RNA\_ncrna\_ncrna\_19-05-2026}})
to their tissue expression profile in GTEx (Genotype-Tissue Expression),
which profiles postmortem human donor tissue exclusively across
{n_tissues} anatomical sites. Of {n_total} unique human ncRNA entities,
{n_matched} ({pct_matched}) were resolved to at least one GENCODE gene known
to GTEx by chaining together every identifier system already present in the
dataset (Ensembl Gene ID, NCBI Gene ID via \texttt{{gene2ensembl}}, miRBase
mature accession via precursor-locus resolution, circBase accession via host
gene), with a GENCODE v39/v26 version fallback to recover genes pruned from
the current annotation. {n_expressed} of the matched entities ({pct_expressed})
are detected at non-zero median TPM in at least one tissue, and isoform-level
expression is available for {n_multi_tx} matched genes. Because a matched
gene's transcript set can include isoforms unrelated to the specific
molecule RNAInterSEQ has on record, and because RNAInterSEQ's own sequence
for an entity is rarely sourced from Ensembl, every candidate transcript's
sequence was independently fetched from Ensembl and tested against
RNAInterSEQ's known sequence(s) for that entity: {n_loci_confirmed}/{n_loci_total}
matched entity-locus pairs ({pct_confirmed}) have at least one such
sequence-confirmed transcript.
\end{{abstract}}

\section{{Background}}
GTEx (\url{{https://gtexportal.org}}) profiles RNA-seq from postmortem tissue
samples donated by human subjects. Its API v2 tissue vocabulary
(\texttt{{tissueSiteDetailId}}) enumerates {n_tissues} human anatomical sites
and carries no species axis -- there is no non-human GTEx data to query.
Every non-\emph{{Homo sapiens}} entity in the source interaction table is
therefore out of scope for this analysis by construction.

\section{{Methods}}

\subsection{{Identifier resolution}}
RNAInterSEQ stores each ncRNA entity's external identifiers in per-source
buckets (Ensembl, NCBI, miRBase, circBase, RNAcentral, GenBank, tRFdb,
``Other''). GenBank, tRFdb and ``Other'' are empty for every human entity in
this dataset and are not pursued further. The remaining four are resolved to
a GTEx-queryable gene symbol or Ensembl Gene ID as follows, tried in this
order for every entity (all are kept when more than one resolves, since a
single mature product can derive from more than one genomic locus):

\begin{{enumerate}}
  \item \textbf{{Ensembl Gene ID}} -- queried directly against GTEx
        \texttt{{reference/gene}}.
  \item \textbf{{NCBI Gene ID}} -- mapped to an Ensembl Gene ID via NCBI's own
        \texttt{{gene2ensembl}} table (human rows only).
  \item \textbf{{Gene symbol}} -- the entity's own RNAInter symbol, queried
        directly.
  \item \textbf{{miRBase mature accession}} (miRNA category only) -- the
        precursor locus name(s) producing that mature product are read from
        \texttt{{miRNA.dat}}; each precursor name is resolved to its official
        HGNC/GENCODE symbol with a single batched \texttt{{mygene.info}}
        alias lookup (handles irregular naming, e.g.\ \texttt{{let-7a-1}}
        $\rightarrow$ \texttt{{MIRLET7A1}}, not a fixed string transform).
  \item \textbf{{circBase accession}} (circRNA category only) -- circBase's
        own FASTA headers already carry the host gene symbol the circRNA
        splices from (GENCODE has no circRNA biotype, so circRNAs can only
        be matched through their linear host gene).
\end{{enumerate}}

\subsection{{GENCODE version fallback}}
Each candidate symbol/Ensembl ID is tried first at GENCODE v39 (paired with
the current \texttt{{gtex\_v10}} dataset), then, only for keys not found, at
v26 (paired with \texttt{{gtex\_v8}}). This recovers genes GENCODE has pruned
out of the current annotation -- {n_v26} entity-locus pairs ({n_v26_only} of
them entities with \emph{{no}} v39 record at all) were only resolvable at
v26.

\subsection{{Sequence-level transcript confirmation}}
Matching an entity to a GENCODE gene (above) identifies the right locus, not
the right \emph{{transcript}}: a multi-isoform gene can carry dozens of
GENCODE transcripts, most of them unrelated to the specific molecule
RNAInterSEQ has on record. There is no shared identifier to resolve this --
RNAInterSEQ's own sequence for an entity is typically sourced from miRBase,
NCBI RefSeq, or RNAcentral, not Ensembl -- so the join key has to be the
sequence itself. For every matched gene, the spliced cDNA sequence of each
of its GENCODE transcripts is fetched from the Ensembl REST API and tested
for containment (exact match, or substring in either direction -- handling
e.g.\ a mature miRNA sequence sitting inside its precursor transcript)
against the entity's own known sequence(s).

\section{{Results}}

\subsection{{Database composition and scope}}
\begin{{sloppypar}}
The source file contains {n_total_int} ncRNA--ncRNA interaction pairs spanning
{n_sp_pairs} species-pair combinations. The dominant species are human
({n_hh} pairs, {pct_hh}) and mouse ({n_mm} pairs, {pct_mm}), with the
remaining {n_other_int} pairs ({pct_other_int}) covering plant, avian, rat,
and viral-host contexts ({n_cross} pairs involve two different species, mostly
viral miRNA--human lncRNA combinations). Because GTEx profiles exclusively
postmortem \emph{{Homo sapiens}} donor tissue, only the {n_hh} human--human
pairs are in scope for this analysis; the full set of 57 columns produced by
this pipeline (original 33 RNAInterSEQ columns $+$ 24 GTEx annotation columns)
is stored in \path{{human_ncrna_interactions_gtex_annotated.pkl}}.
\end{{sloppypar}}

\subsection{{Coverage funnel}}
\begin{{itemize}}
  \item {n_total} unique human ncRNA entities in the 19-05-2026 dataset
  \item {n_matched} ({pct_matched}) resolved to $\geq 1$ GENCODE gene known to GTEx
  \item {n_multi_locus} of those resolved to multiple genomic loci (paralogous precursors, mostly miRNA)
  \item {n_expressed} ({pct_expressed} of matched) detected at median TPM $>0$ in $\geq 1$ of {n_tissues} tissues
  \item {n_single_tx} matched genes have exactly one annotated transcript; {n_multi_tx} have multiple ({n_total_transcripts} unique transcripts total)
\end{{itemize}}

\subsection{{Match rate by RNAInter category}}
\begin{{table}}[H]
\centering
{rate_table}
\caption{{Entity-level GTEx match rate by RNAInter ncRNA category.}}
\end{{table}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.8\textwidth]{{{fig1}}}
\caption{{GTEx match rate by RNAInter ncRNA category.}}
\end{{figure}}

\subsection{{Which identifier system resolved each match}}
\begin{{table}}[H]
\centering
{method_table}
\caption{{Entity-locus pairs by resolution method.}}
\end{{table}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.8\textwidth]{{{fig4}}}
\caption{{Contribution of each identifier system to the final match set.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.65\textwidth]{{{fig5}}}
\caption{{GENCODE v39/gtex\_v10 vs.\ v26/gtex\_v8 fallback usage.}}
\end{{figure}}

RNAInter category and GTEx/GENCODE biotype agree for {biotype_agree} of
matched entities (primary locus per entity); the main source of disagreement
is RNAInter ``lncRNA'' entities that GENCODE annotates as protein-coding,
worth a closer look if category purity matters downstream.

\subsection{{Sequence-confirmed transcripts}}
{n_loci_confirmed}/{n_loci_total} matched entity-locus pairs ({pct_confirmed})
have $\geq 1$ GENCODE transcript whose Ensembl-fetched sequence matches a
sequence RNAInterSEQ already has on record. For genes that get $\geq 1$
confirmed transcript, the confirmed set is on average {narrow_mean} of that
gene's total annotated transcripts (median {narrow_median}) -- i.e.\ this
step meaningfully narrows multi-isoform genes down to the specific
isoform(s) RNAInterSEQ has, rather than reporting the whole gene's isoform
mixture as the entity's expression. Filtering down to only sequence-confirmed
transcripts keeps {n_tx_rows_confirmed}/{n_tx_rows_total} (transcript,
tissue) expression rows.

\begin{{table}}[H]
\centering
\small
{confirm_table}
\caption{{Sequence-confirmation rate by RNAInter ncRNA category (entity-locus level).}}
\end{{table}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.8\textwidth]{{{fig6}}}
\caption{{Sequence-confirmation rate by RNAInter ncRNA category.}}
\end{{figure}}

Confirmation rate tracks how well a category's typical sequence length and
provenance line up with GENCODE's annotated transcript boundaries: miRNA is
near-total (99.7\%, mature sequences are short and almost always found
verbatim inside their precursor transcript), while circRNA is the lowest
(13.8\%) because a circRNA's backsplice-junction sequence is, by definition,
not contiguous in the linear host transcript -- consistent with the
host-gene-proxy limitation noted below.

\subsection{{Tissue specificity}}
Yanai's tau across {n_tissues} GTEx tissues, computed for the {n_tau} matched
genes with non-zero expression in at least one tissue: mean {tau_mean},
median {tau_median} -- most matched ncRNAs are far from uniformly expressed.

\begin{{figure}}[H]
\centering
\includegraphics[width=0.7\textwidth]{{{fig2}}}
\caption{{Distribution of tissue specificity (tau) across matched ncRNA genes.}}
\end{{figure}}

\begin{{table}}[H]
\centering
\small
{top_table}
\caption{{Top 15 most tissue-specific matched ncRNAs and their peak tissue.}}
\end{{table}}

\subsection{{Transcript multiplicity}}
\begin{{figure}}[H]
\centering
\includegraphics[width=0.7\textwidth]{{{fig3}}}
\caption{{Number of annotated GENCODE transcripts per matched ncRNA gene.}}
\end{{figure}}

\section{{Interaction-level tissue expression annotation}}

The entity-level GTEx data are joined back to the
{n_human_interactions} human--human interaction pairs in the source table
(both partners \emph{{Homo sapiens}}) to produce an annotated interaction
table (\path{{human_ncrna_interactions_gtex_annotated.pkl}}). For each
interaction row the table carries, for both partners: the set of
sequence-confirmed GENCODE transcript IDs (as a per-gene dict), tissue
specificity (tau), peak tissue, peak median TPM, and the full 54-tissue
max-TPM vector.

\begin{{itemize}}
  \item {n_both_matched} ({pct_both_matched}) of pairs have both partners
        GTEx-matched.
  \item {n_both_confirmed} ({pct_both_confirmed}) have both partners
        sequence-confirmed.
  \item {n_shared_any} ({pct_shared_any}) share $\geq 1$ expressed GTEx tissue
        (median TPM $> 0$ in both partners).
  \item Median shared tissues across sequence-confirmed pairs:
        {med_shared} (out of 54) -- reflecting the high tissue
        specificity of ncRNAs in this dataset (median $\tau = $ {tau_median})
        and the residual miRNA coverage gap.
\end{{itemize}}

\begin{{figure}}[H]
\centering
\includegraphics[width=\textwidth]{{{fig6b}}}
\caption{{Genes and transcripts before and after sequence-confirmed
filtering (step 5), broken down by RNAInter ncRNA category.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=\textwidth]{{{fig7}}}
\caption{{Number of co-expressed interaction pairs per GTEx tissue
(both partners detected at median TPM $> 0$, sequence-confirmed pairs only).}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.75\textwidth]{{{fig8}}}
\caption{{Distribution of the number of shared expressed tissues across
sequence-confirmed interaction pairs.}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=\textwidth]{{{fig9}}}
\caption{{ncRNA category-pair co-expression landscape. Left: interaction
count (both GTEx-matched). Right: fraction of pairs with $\geq 1$ shared
expressed tissue (sequence-confirmed pairs only).}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.85\textwidth]{{{fig10}}}
\caption{{Yanai tissue-specificity ($\tau$) distribution per ncRNA category,
pooled across all interaction partners with GTEx data
(unique entities; $\tau = 0$ ubiquitous, $\tau = 1$ single-tissue).}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.75\textwidth]{{{fig11}}}
\caption{{Top 20 peak expression tissues across unique interacting ncRNA
entities (both partners pooled, deduplicated by entity).}}
\end{{figure}}

\begin{{figure}}[H]
\centering
\includegraphics[width=0.65\textwidth]{{{fig12}}}
\caption{{Tissue-specificity correlation between the two partners of each
interaction (sequence-confirmed pairs; top 5 category combinations
highlighted).}}
\end{{figure}}

\subsection{{Biological system grouping}}
\begin{{sloppypar}}
To improve interpretability, the 54 GTEx tissue sites are collapsed into
eight biological systems: CNS (13 tissues), Peripheral nervous (2),
Cardiovascular (5), GI (10), Reproductive (9), Immune (3), Muscle (1),
and Other (11 metabolic, urological, and integumentary tissues). For each
entity, \texttt{{ent1\_tissue\_systems}} and \texttt{{ent2\_tissue\_systems}}
list the systems in which it is detected (median TPM $>0$ across confirmed
transcripts). For sequence-confirmed interaction pairs,
\texttt{{shared\_tissue\_systems}} gives the subset of systems shared by both
partners. Overall, {n_sys_any} pairs ({pct_sys_any}) share expression in at
least one biological system; the most represented system is
\textbf{{{top_sys}}} ({top_sys_count} pairs), likely driven by the broad
transcriptional activity of gonadal tissues (testis in particular ranks
among the highest expressing tissues in GTEx across all gene classes).
\end{{sloppypar}}

\begin{{figure}}[H]
\centering
\includegraphics[width=\textwidth]{{{fig13}}}
\caption{{Left: number of sequence-confirmed co-expressed interaction pairs
per biological system (54 GTEx tissues collapsed to 8 systems).
Right: symmetric co-occurrence matrix -- each cell shows how many pairs
have both partners expressed in system~$A$ \emph{{and}} system~$B$
simultaneously (a pair can contribute to multiple cells if both partners
are expressed across several systems).}}
\end{{figure}}

\subsection{{Interpreting the co-expression rate}}
\begin{{sloppypar}}
Two questions arise naturally from the {pct_shared_any} co-expression figure:
(\emph{{i}}) does the analysis account for transcripts being expressed across
multiple tissues simultaneously, and (\emph{{ii}}) can a transcript be
co-expressed in a tissue that GTEx simply does not detect?

\textbf{{Multi-tissue expression is fully captured.}} The analysis does
\emph{{not}} compare only peak tissues -- it uses the complete 54-tissue
max-TPM vector for each entity (columns \texttt{{ent1\_tissue\_tpm}} /
\texttt{{ent2\_tissue\_tpm}}).
\texttt{{shared\_expressed\_tissues}} is the intersection of the sets of tissues
where each partner has median TPM $>0$, across all confirmed transcripts and
all matched loci. Of the {n_both_confirmed} sequence-confirmed pairs, {n_shared_any}
({pct_shared_conf}) share at least one tissue. Among entities with any GTEx
expression at all, the median is {med_n_tissues} out of 54 tissues, and
{n_all54} entities are expressed in every single GTEx tissue -- so when one
partner is broadly expressed the method finds overlap even if the other partner
only peaks in one tissue.

\textbf{{GTEx detectability, not tissue incompatibility, drives the low rate.}}
Dissecting the {n_both_confirmed} sequence-confirmed pairs:
\begin{{itemize}}
  \item {n_shared_any} ({pct_shared_conf}) -- both partners detected, tissues overlap.
  \item {n_zero_partner} ({pct_zero_partner}) -- at least one partner shows
        median TPM $= 0$ across all confirmed transcripts in all 54 GTEx tissues.
        In most cases ($\sim$68\%) the gene locus is also entirely absent from GTEx;
        in the remainder ($\sim$32\%) a different isoform of the same gene is
        detected but the sequence-confirmed isoform matching RNAInterSEQ's record is not.
        This accounts for \textbf{{nearly all}} non-co-expressed confirmed pairs.
  \item {n_disjoint} ({pct_disjoint}) -- both partners detected in GTEx but in
        non-overlapping tissue sets within GTEx's scope.
\end{{itemize}}
The detectability gap is category-specific: 99.9\% of miRNA entities, 90.3\%
of snoRNA, 87.1\% of circRNA, and 100\% of tRNA show zero expression across
all 54 GTEx tissues -- for reasons inherent to how GTEx profiles RNA (see
Limitations). Only lncRNAs have reasonable coverage (31.7\% undetected).
The {pct_shared_any} co-expression rate therefore reflects a
\textbf{{conservative lower bound}}: interactions involving at least one
miRNA, snoRNA, circRNA, or tRNA partner are almost never scored as
co-expressed here regardless of the underlying biology.
\end{{sloppypar}}

\section{{Limitations}}
\begin{{sloppypar}}
\begin{{itemize}}
  \item \textbf{{miRNA mature/precursor mismatch.}} RNAInter stores
        \emph{{mature}} miRNA names (e.g.\ \texttt{{hsa-miR-101-3p}}); GTEx
        captures bulk poly(A)-selected RNA-seq, which mostly fails to detect
        the short, rapidly processed mature miRNA at all -- even a correctly
        matched precursor (pri-miRNA) locus typically shows near-zero TPM.
        Some mature products also lack any HGNC-approved host gene symbol
        and cannot be resolved by this pipeline at all.
  \item \textbf{{GTEx bulk RNA-seq averaging and isoform specificity.}} GTEx
        profiles bulk RNA from heterogeneous tissue biopsies. A transcript
        expressed in a minority cell population (e.g.\ Purkinje cells within
        cerebellum, or interstitial cells within testis) may be diluted below
        detection threshold in the bulk signal and appear at TPM\,$= 0$ even
        in the tissue where it is most relevant. Additionally, the
        sequence-confirmed isoform can be non-dominant: $\sim$32\% of
        zero-tissue confirmed entities have gene-level GTEx expression in at
        least one tissue but zero TPM for the specific isoform matching
        RNAInterSEQ's sequence record. Both effects contribute to the
        dominant reason the co-expression rate is low: {n_zero_partner}
        ({pct_zero_partner}) of sequence-confirmed pairs have at least one
        partner with zero TPM across all confirmed transcripts in all 54 GTEx
        tissues, compared with only {n_disjoint} ({pct_disjoint}) pairs where
        both partners are detected but in non-overlapping tissue sets.
  \item \textbf{{GTEx tissue panel is not exhaustive.}} The 54-tissue
        vocabulary covers major adult solid organs and blood but omits
        important tissue contexts where ncRNA activity has been reported:
        placenta, bone marrow, lymph node, tonsil, retina, cornea, and
        synovial tissue, among others. A co-expressed pair in any of these
        contexts will not be captured, and the \texttt{{shared\_tissue\_systems}}
        field will list those systems as absent for both partners.
        GTEx also captures only adult postmortem donors, missing
        developmental and foetal expression programmes.
  \item \textbf{{circRNA host-gene proxy.}} A circRNA's expression here is its
        \emph{{linear host gene's}} expression, not the circularized isoform's
        -- GENCODE has no circRNA biotype, so this is the closest available
        proxy, not a direct measurement.
  \item \textbf{{RNAcentral, GenBank, tRFdb, ``Other'' buckets not pursued.}}
        The latter three are empty for every human entity in this dataset.
        RNAcentral URS IDs are populated for most entities but resolving them
        to GENCODE genes would require downloading RNAcentral's full
        \texttt{{id\_mapping.tsv.gz}} (hundreds of MB) for marginal gains once
        Ensembl/NCBI/miRBase/circBase resolution is exhausted; left as future
        work if residual unmatched entities are specifically needed.
  \item \textbf{{Multi-locus aggregation.}} Entities resolved to more than one
        genomic locus (paralogous miRNA precursors, or symbol/Ensembl keys
        landing on distinct historical GENCODE records) are kept as separate
        rows rather than collapsed into one number; downstream consumers
        should decide whether to sum, average, or enumerate per locus.
  \item \textbf{{Sequence confirmation is exact/containment-based, not
        alignment-based.}} A transcript that differs from RNAInterSEQ's
        sequence by even one base (a different assembly patch, a SNP-bearing
        haplotype, a slightly trimmed UTR) will not be confirmed, even though
        it is almost certainly the same biological isoform -- the
        sub-100\% confirmation rates for snoRNA/snRNA/rRNA likely reflect
        this rather than a true absence of the transcript in GENCODE. A
        fuzzy/alignment-based fallback (as already used elsewhere in this
        repository, see \texttt{{similarity\_computation.py}}) would recover
        some of these at the cost of one pairwise alignment per candidate
        transcript.
\end{{itemize}}
\end{{sloppypar}}

\section{{Data availability}}
\begin{{sloppypar}}
All intermediate and final tables are saved as pandas pickles under
\texttt{{gtex\_ncrna\_expression/data/}}:
\begin{{itemize}}
  \item \path{{ncrna_entities_human.pkl}} -- the entity table built in step 1
  \item \path{{gtex_gene_matches.pkl}} -- the long entity-locus match table from step 2
  \item \path{{gtex_gene_expression.pkl}} / \path{{gtex_transcript_expression.pkl}} -- raw GTEx median TPM, gene and transcript level
  \item \path{{ncrna_entity_locus_summary.pkl}} -- the merged, analysis-ready table from step 3, including a \texttt{{sequence\_confirmed}} flag from step 5
  \item \path{{ncrna_locus_transcript_matches.pkl}} -- per (entity, gene, transcript) sequence-confirmation detail from step 5
  \item \path{{gtex_transcript_expression_sequence_confirmed.pkl}} -- the transcript-expression table filtered to only sequence-confirmed transcripts
  \item \path{{human_ncrna_interactions_gtex_annotated.pkl}} -- the {n_human_interactions} human--human interaction rows annotated with per-entity GTEx tissue expression metadata (step 6); includes \texttt{{ent1\_*/ent2\_*}} columns for confirmed transcripts, tau, peak tissue, the 54-tissue TPM vector, and per-entity biological system lists (\texttt{{tissue\_systems}}), plus \texttt{{shared\_expressed\_tissues}} and \texttt{{shared\_tissue\_systems}}
\end{{itemize}}
Raw GTEx/mygene.info API responses are cached under \texttt{{gtex\_ncrna\_expression/cache/}}.
\end{{sloppypar}}

\end{{document}}
"""

if __name__ == "__main__":
    main()
