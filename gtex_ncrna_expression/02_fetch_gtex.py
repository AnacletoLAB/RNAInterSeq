"""
Step 2: resolve each human ncRNA entity to one or more GENCODE genes known
to GTEx, then pull median tissue expression at both the gene level and the
transcript level for everything that matched.

RNAInterSEQ stores four identifier systems for these entities that are
useful here (the rest -- GenBank, tRFdb, "Other" -- are empty for every
human entity in the 19-05-2026 dataset, so there is nothing to resolve from
them): Ensembl Gene ID, NCBI Gene ID, miRBase mature accession, circBase
accession. Each gets its own resolution path down to a GTEx-queryable
symbol or Ensembl ID:

  Ensembl Gene ID   -> queried directly.
  NCBI Gene ID      -> NCBI's own gene2ensembl table (data/gene2ensembl,
                        human rows only) gives the Ensembl Gene ID.
  miRBase mature ID -> data/miRNA.dat gives the precursor locus name(s) that
                        produce that mature product (often >1: paralogous
                        loci such as hsa-mir-101-1 / hsa-mir-101-2 yield the
                        identical mature sequence); a single batched
                        mygene.info alias lookup turns each precursor name
                        into its official HGNC/GENCODE symbol (handles
                        irregular naming such as let-7a-1 -> MIRLET7A1).
  circBase ID       -> circBase's own FASTA headers
                        (data/circdata/circBase_9606.fa) already carry the
                        host gene symbol the circRNA is spliced from.

A symbol/Ensembl ID is additionally tried at two GENCODE versions: v39
(paired with the current gtex_v10 dataset) first, falling back to v26
(paired with gtex_v8) for genes GENCODE pruned out of v39 -- e.g. several
miRNA host genes such as MIRLET7A1 only have a v26 record. Whichever
version a gene is found at determines which GTEx dataset its expression is
pulled from.

Because several of these paths (miRBase precursors especially) can resolve
one RNAInter entity to more than one genomic locus, the matches table is
long: one row per (entity, locus) pair, not one row per entity.
"""
import pickle
from pathlib import Path

import pandas as pd

import gtex_api as gtex
import id_resolvers as idr

ROOT = Path(__file__).parent
ENTITIES_PKL = ROOT / "data" / "ncrna_entities_human.pkl"
MATCHES_PKL = ROOT / "data" / "gtex_gene_matches.pkl"
GENE_EXPR_PKL = ROOT / "data" / "gtex_gene_expression.pkl"
TX_EXPR_PKL = ROOT / "data" / "gtex_transcript_expression.pkl"


def add_search_keys(ents):
    """
    Build, per entity, an ordered list of (search_key, method) candidates to
    try against GTEx reference/gene. Keys are tried in priority order and
    every hit is kept (not just the first) for methods that can legitimately
    point at more than one genomic locus (miRBase precursors).
    """
    ncbi_to_ensembl = idr.load_ncbi_to_ensembl_human()
    circbase_host = idr.load_circbase_host_genes_human()
    mirbase_mature_to_precursor = idr.load_mirbase_mature_to_precursor_human()

    # only resolve the miRBase precursor names we actually need (entities
    # carrying a miRBase mature accession), not the whole flat file
    needed_precursors = set()
    for accs in ents["mirbase_ids"]:
        for acc in accs:
            needed_precursors.update(mirbase_mature_to_precursor.get(acc, []))
    print(f"Resolving {len(needed_precursors):,} miRBase precursor names via mygene.info ...")
    precursor_resolution = idr.resolve_precursor_symbols(needed_precursors)

    search_keys = []
    for ent in ents.itertuples(index=False):
        keys = []  # list of (value, method)
        if ent.ensembl_gene_id:
            keys.append((ent.ensembl_gene_id, "ensembl_id"))
        if ent.ncbi_gene_id and ent.ncbi_gene_id in ncbi_to_ensembl:
            eg = ncbi_to_ensembl[ent.ncbi_gene_id]
            if eg != ent.ensembl_gene_id:
                keys.append((eg, "ncbi_to_ensembl"))
        keys.append((ent.symbol.upper(), "symbol"))
        if ent.category == "miRNA":
            for acc in ent.mirbase_ids:
                for precursor in mirbase_mature_to_precursor.get(acc, []):
                    res = precursor_resolution.get(precursor)
                    if res and res.get("symbol"):
                        keys.append((res["symbol"].upper(), "mirbase_precursor"))
        if ent.category == "circRNA":
            for cb_id in ent.circbase_ids:
                host = circbase_host.get(cb_id)
                if host:
                    keys.append((host.upper(), "circbase_host"))
        search_keys.append(keys)
    ents = ents.copy()
    ents["search_keys"] = search_keys
    return ents


def match_entities(ents):
    ents = add_search_keys(ents)

    all_keys = sorted({k for keys in ents["search_keys"] for k, _ in keys})
    print(f"Looking up {len(all_keys):,} unique search keys against GTEx "
          f"reference/gene (v39/gtex_v10 first, v26/gtex_v8 fallback) ...")
    resolved = gtex.lookup_genes_with_fallback(all_keys)
    print(f"  resolved {len(resolved):,}/{len(all_keys):,} keys")

    rows = []
    for ent in ents.itertuples(index=False):
        seen_gencode_ids = set()
        any_hit = False
        for key, method in ent.search_keys:
            for rec, dataset_id in resolved.get(key, []):
                if rec["gencodeId"] in seen_gencode_ids:
                    continue  # same locus reached via more than one search key
                seen_gencode_ids.add(rec["gencodeId"])
                any_hit = True
                rows.append({
                    "entity_key": ent.entity_key,
                    "symbol": ent.symbol,
                    "category": ent.category,
                    "ncbi_gene_id": ent.ncbi_gene_id,
                    "matched": True,
                    "match_method": method,
                    "gencode_id": rec["gencodeId"],
                    "gencode_version": rec["gencodeVersion"],
                    "dataset_id": dataset_id,
                    "gtex_gene_symbol": rec["geneSymbol"],
                    "gtex_gene_type": rec["geneType"],
                    "gtex_entrez_id": rec["entrezGeneId"],
                })
        if not any_hit:
            rows.append({
                "entity_key": ent.entity_key, "symbol": ent.symbol,
                "category": ent.category, "ncbi_gene_id": ent.ncbi_gene_id,
                "matched": False, "match_method": None, "gencode_id": None,
                "gencode_version": None, "dataset_id": None,
                "gtex_gene_symbol": None, "gtex_gene_type": None, "gtex_entrez_id": None,
            })
    return pd.DataFrame(rows)


def main():
    with open(ENTITIES_PKL, "rb") as f:
        ents = pickle.load(f)

    matches = match_entities(ents)
    matches.to_pickle(MATCHES_PKL)

    n_entities = ents["entity_key"].nunique()
    matched_entities = matches.loc[matches["matched"], "entity_key"].nunique()
    n_loci = len(matches[matches["matched"]])
    print(f"\nMatched {matched_entities:,}/{n_entities:,} human ncRNA entities "
          f"to >=1 GTEx gene ({matched_entities/n_entities:.1%}), "
          f"{n_loci:,} entity-locus pairs total "
          f"({n_loci - matched_entities:,} from multi-locus resolution).")

    print("\nMatch rate by category (entity-level):")
    per_entity = matches.groupby("entity_key").agg(
        category=("category", "first"), matched=("matched", "any"))
    rate = per_entity.groupby("category")["matched"].agg(n_total="count", n_matched="sum")
    rate["match_rate"] = rate["n_matched"] / rate["n_total"]
    print(rate.sort_values("n_total", ascending=False).to_string())

    print("\nMatch method (entity-locus pairs):")
    print(matches.loc[matches["matched"], "match_method"].value_counts().to_string())
    print("\nGENCODE version used (entity-locus pairs):")
    print(matches.loc[matches["matched"], "gencode_version"].value_counts().to_string())

    gene_expr_frames, tx_expr_frames = [], []
    for dataset_id, group in matches[matches["matched"]].groupby("dataset_id"):
        gencode_ids = group["gencode_id"].unique().tolist()
        print(f"\nFetching expression for {len(gencode_ids):,} genes from {dataset_id} ...")
        gene_expr_frames.append(pd.DataFrame(
            gtex.median_gene_expression(gencode_ids, dataset_id=dataset_id)))
        tx_expr_frames.append(pd.DataFrame(
            gtex.median_transcript_expression(gencode_ids, dataset_id=dataset_id)))

    gene_expr_df = pd.concat(gene_expr_frames, ignore_index=True)
    tx_expr_df = pd.concat(tx_expr_frames, ignore_index=True)
    gene_expr_df.to_pickle(GENE_EXPR_PKL)
    tx_expr_df.to_pickle(TX_EXPR_PKL)
    print(f"\n{len(gene_expr_df):,} (gene, tissue) rows -> {GENE_EXPR_PKL}")
    print(f"{len(tx_expr_df):,} (transcript, tissue) rows, "
          f"{tx_expr_df['transcriptId'].nunique():,} unique transcripts -> {TX_EXPR_PKL}")


if __name__ == "__main__":
    main()
