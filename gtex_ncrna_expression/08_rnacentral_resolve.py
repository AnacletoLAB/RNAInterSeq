"""
Step 8: Extend GTEx gene matches using the RNAcentral id_mapping database.

The 2,297 entities that failed every resolution path in step 2 are revisited
here via three DuckDB-powered paths through the RNAcentral release-25 flat
files:

  Path A  URS (already in entity table) → id_map ENSEMBL_GENCODE → ENSG
          Bulk: one JOIN of all URS IDs against id_mapping, reads the parent
          ENSG directly from the `misc` column (e.g. "ENSG00000226803.10").
          Covers 95 % of unmatched entities that already carry URS IDs.

  Path B  miRBase MIMAT/MI → id_map MIRBASE → URS → ENSEMBL_GENCODE → ENSG
          Two-hop JOIN; fills gaps where the entity's stored URS is the
          mature-product sequence (no Ensembl entry) but the miRBase accession
          points at a precursor URS that does have one.

  Path C  Sequence MD5 → md5_map → URS → ENSEMBL_GENCODE → ENSG
          For any entity whose stored URS and miRBase IDs both failed: compute
          MD5 of the entity's own sequence from RNAInterSEQ, look it up in
          md5_map, then chain to ENSEMBL_GENCODE. Works for mature miRNAs
          whose exact sequence is catalogued as a separate URS in RNAcentral.

After resolution, newly found ENSG IDs are verified against the GTEx API
(v39/gtex_v10 primary, v26/gtex_v8 fallback -- same logic as step 2), and
gene-level + transcript-level median expression are fetched and appended to
the existing pickles so that downstream steps (05 → 03 → 06 → 07 → 04) can
simply be re-run without any other changes.

Prerequisites:
  rnacentral-etl/current_release/id_mapping/id_mapping.tsv.gz   (1.8 GB)
  rnacentral-etl/current_release/md5/md5.tsv.gz                  (1.0 GB)

  rsync -avz rsync.ebi.ac.uk::pub/databases/RNAcentral/releases/25.0/id_mapping/ current_release/id_mapping/
  rsync -avz rsync.ebi.ac.uk::pub/databases/RNAcentral/releases/25.0/md5/        current_release/md5/

Run from the repo root with:
  /home/mnicolini/miniconda3/envs/vGenRNA/bin/python gtex_ncrna_expression/08_rnacentral_resolve.py
"""
import hashlib
import sys
from pathlib import Path

import duckdb
import pandas as pd

ROOT          = Path(__file__).parent
DATA          = ROOT / "data"
RC_DIR        = ROOT.parent / "rnacentral-etl"
DB_PATH       = RC_DIR / "rnacentral.duckdb"
ID_MAP_GZ     = RC_DIR / "current_release" / "id_mapping" / "id_mapping.tsv.gz"
MD5_GZ        = RC_DIR / "current_release" / "md5" / "md5.tsv.gz"

ENTITIES_PKL  = DATA / "ncrna_entities_human.pkl"
MATCHES_PKL   = DATA / "gtex_gene_matches.pkl"
GENE_EXPR_PKL = DATA / "gtex_gene_expression.pkl"
TX_EXPR_PKL   = DATA / "gtex_transcript_expression.pkl"

sys.path.insert(0, str(ROOT))
import gtex_api as gtex


# ---------------------------------------------------------------------------
# DuckDB setup
# ---------------------------------------------------------------------------

def _build_db():
    """Create rnacentral.duckdb with id_map and md5_map from the release-25 flat files."""
    print(f"Building DuckDB at {DB_PATH} ...")
    print("  Loading id_mapping.tsv.gz (~1.8 GB compressed; a few minutes) …")
    con = duckdb.connect(str(DB_PATH))
    con.execute("""
        CREATE TABLE IF NOT EXISTS id_map AS
        SELECT
            column0 AS RNAcentral_ID,
            column1 AS external_db,
            column2 AS external_id,
            column3 AS taxid,
            column4 AS rna_type,
            column5 AS misc
        FROM read_csv_auto(?, delim='\t', header=False, ignore_errors=True)
    """, [str(ID_MAP_GZ)])
    print(f"  id_map: {con.execute('SELECT COUNT(*) FROM id_map').fetchone()[0]:,} rows")

    print("  Creating indexes on id_map …")
    con.execute("CREATE INDEX IF NOT EXISTS idx_urs        ON id_map (RNAcentral_ID)")
    con.execute("CREATE INDEX IF NOT EXISTS idx_db_ext_tax ON id_map (external_db, external_id, taxid)")
    con.execute("CREATE INDEX IF NOT EXISTS idx_urs_db     ON id_map (RNAcentral_ID, external_db, taxid)")

    print("  Loading md5.tsv.gz (~1.0 GB compressed) …")
    con.execute("""
        CREATE TABLE IF NOT EXISTS md5_map AS
        SELECT column0 AS RNAcentral_ID, column1 AS md5
        FROM read_csv_auto(?, delim='\t', header=False, ignore_errors=True)
    """, [str(MD5_GZ)])
    print(f"  md5_map: {con.execute('SELECT COUNT(*) FROM md5_map').fetchone()[0]:,} rows")
    con.execute("CREATE INDEX IF NOT EXISTS idx_md5     ON md5_map (md5)")
    con.execute("CREATE INDEX IF NOT EXISTS idx_md5_urs ON md5_map (RNAcentral_ID)")

    con.close()
    print(f"  DuckDB ready.")


def _open_db(read_only=True):
    return duckdb.connect(str(DB_PATH), read_only=read_only)


# ---------------------------------------------------------------------------
# Entity preparation
# ---------------------------------------------------------------------------

def _flatten_urs(rnacentral_ids):
    """Normalise the potentially nested list of URS IDs stored in the entity table."""
    out = []
    for item in rnacentral_ids:
        if isinstance(item, list):
            out.extend(item)
        else:
            out.append(item)
    return [u for u in out if isinstance(u, str) and u.startswith("URS")]


def _seq_md5(seq):
    return hashlib.md5(seq.upper().replace("U", "T").encode()).hexdigest()


def load_unmatched(ents, matches):
    """Return entity rows that have no matched=True row in gtex_gene_matches."""
    matched_keys = set(matches.loc[matches["matched"], "entity_key"])
    return ents[~ents["entity_key"].isin(matched_keys)].copy()


# ---------------------------------------------------------------------------
# Resolution paths (all return DataFrames with columns:
#   entity_key, symbol, category, ncbi_gene_id, gencode_id_raw, rna_type, match_method)
# ---------------------------------------------------------------------------

def _path_a(con, unmatched):
    """URS already in entity table → ENSEMBL_GENCODE misc field → ENSG (no version)."""
    # build URS → entity mapping (one URS can belong to multiple entities)
    urs_to_ents = {}
    for _, row in unmatched.iterrows():
        for urs in _flatten_urs(row["rnacentral_ids"]):
            urs_to_ents.setdefault(urs, []).append(row)

    if not urs_to_ents:
        return pd.DataFrame()

    urs_df = pd.DataFrame({"urs_id": list(urs_to_ents.keys())})
    con.register("_urs_a", urs_df)
    hits = con.execute("""
        SELECT DISTINCT
            u.urs_id,
            regexp_replace(i.misc, '\\..*$', '') AS gencode_id_raw,
            i.rna_type
        FROM _urs_a u
        JOIN id_map i ON i.RNAcentral_ID = u.urs_id
        WHERE i.taxid       = '9606'
          AND i.external_db = 'ENSEMBL_GENCODE'
          AND i.misc        IS NOT NULL
          AND i.misc        LIKE 'ENSG%'
    """).df()
    con.unregister("_urs_a")

    rows = []
    for _, h in hits.iterrows():
        for ent_row in urs_to_ents[h["urs_id"]]:
            rows.append({
                "entity_key":    ent_row["entity_key"],
                "symbol":        ent_row["symbol"],
                "category":      ent_row["category"],
                "ncbi_gene_id":  ent_row["ncbi_gene_id"],
                "gencode_id_raw": h["gencode_id_raw"],
                "rna_type":      h["rna_type"],
                "match_method":  "rnacentral_urs",
            })
    return pd.DataFrame(rows)


def _path_b(con, still_unmatched):
    """miRBase MIMAT/MI accession → URS → ENSEMBL_GENCODE → ENSG (two-hop JOIN)."""
    mirna = still_unmatched[still_unmatched["category"] == "miRNA"]
    mirbase_to_ents = {}
    for _, row in mirna.iterrows():
        for acc_list in row["mirbase_ids"]:
            accs = acc_list if isinstance(acc_list, list) else [acc_list]
            for acc in accs:
                if isinstance(acc, str) and acc.startswith(("MI", "MIMAT")):
                    mirbase_to_ents.setdefault(acc, []).append(row)

    if not mirbase_to_ents:
        return pd.DataFrame()

    mir_df = pd.DataFrame({"mirbase_id": list(mirbase_to_ents.keys())})
    con.register("_mir_b", mir_df)
    hits = con.execute("""
        SELECT DISTINCT
            m.mirbase_id,
            regexp_replace(e.misc, '\\..*$', '') AS gencode_id_raw,
            e.rna_type
        FROM _mir_b m
        JOIN id_map mb ON mb.external_db = 'MIRBASE'
                      AND mb.external_id = m.mirbase_id
                      AND mb.taxid       = '9606'
        JOIN id_map e  ON e.RNAcentral_ID = mb.RNAcentral_ID
                      AND e.external_db   = 'ENSEMBL_GENCODE'
                      AND e.taxid         = '9606'
                      AND e.misc          IS NOT NULL
                      AND e.misc          LIKE 'ENSG%'
    """).df()
    con.unregister("_mir_b")

    rows = []
    for _, h in hits.iterrows():
        for ent_row in mirbase_to_ents[h["mirbase_id"]]:
            rows.append({
                "entity_key":    ent_row["entity_key"],
                "symbol":        ent_row["symbol"],
                "category":      ent_row["category"],
                "ncbi_gene_id":  ent_row["ncbi_gene_id"],
                "gencode_id_raw": h["gencode_id_raw"],
                "rna_type":      h["rna_type"],
                "match_method":  "rnacentral_mirbase",
            })
    return pd.DataFrame(rows)


def _path_c(con, still_unmatched):
    """Sequence MD5 → md5_map → URS → ENSEMBL_GENCODE → ENSG."""
    seq_to_ents = {}
    for _, row in still_unmatched.iterrows():
        for seq in row["known_sequences"]:
            if isinstance(seq, str) and seq:
                md5 = _seq_md5(seq)
                seq_to_ents.setdefault(md5, []).append(row)

    if not seq_to_ents:
        return pd.DataFrame()

    md5_df = pd.DataFrame({"md5": list(seq_to_ents.keys())})
    con.register("_seq_c", md5_df)
    hits = con.execute("""
        SELECT DISTINCT
            s.md5,
            regexp_replace(e.misc, '\\..*$', '') AS gencode_id_raw,
            e.rna_type
        FROM _seq_c s
        JOIN md5_map md ON md.md5              = s.md5
        JOIN id_map e   ON e.RNAcentral_ID     = md.RNAcentral_ID
                       AND e.external_db       = 'ENSEMBL_GENCODE'
                       AND e.taxid             = '9606'
                       AND e.misc              IS NOT NULL
                       AND e.misc              LIKE 'ENSG%'
    """).df()
    con.unregister("_seq_c")

    rows = []
    for _, h in hits.iterrows():
        for ent_row in seq_to_ents[h["md5"]]:
            rows.append({
                "entity_key":    ent_row["entity_key"],
                "symbol":        ent_row["symbol"],
                "category":      ent_row["category"],
                "ncbi_gene_id":  ent_row["ncbi_gene_id"],
                "gencode_id_raw": h["gencode_id_raw"],
                "rna_type":      h["rna_type"],
                "match_method":  "rnacentral_seq",
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# GTEx verification + expression fetch
# ---------------------------------------------------------------------------

def verify_and_fetch(resolved_df, existing_gene_expr, existing_tx_expr):
    """
    For each candidate ENSG from the three paths:
      1. Verify it actually exists in GTEx (v39 primary / v26 fallback).
      2. Fetch gene-level and transcript-level median expression.
      3. Return new rows for gtex_gene_matches, gene_expr, tx_expr.

    Already-fetched genes (present in existing_gene_expr) are not re-fetched
    from the API — their expression rows are inherited as-is.
    """
    candidate_ensg = resolved_df["gencode_id_raw"].unique().tolist()
    already_have   = set(existing_gene_expr["gencodeId"].str.split(".").str[0])

    to_verify = [g for g in candidate_ensg if g not in already_have]
    print(f"  Verifying {len(to_verify):,} new ENSG IDs against GTEx API …")
    gtex_records = gtex.lookup_genes_with_fallback(to_verify)

    # base-ENSG → (versioned gencodeId, geneSymbol, geneType, entrezGeneId, gencode_version, dataset_id)
    ensg_to_gtex = {}
    for base_id, recs in gtex_records.items():
        for rec, dataset_id in recs:
            base = rec["gencodeId"].split(".")[0]
            ensg_to_gtex.setdefault(base, []).append(
                (rec["gencodeId"], rec["geneSymbol"], rec["geneType"],
                 rec.get("entrezGeneId"), rec["gencodeVersion"], dataset_id))

    # also look up already-cached genes via existing matches index
    print(f"  {len(gtex_records):,}/{len(to_verify):,} new ENSG IDs found in GTEx")

    # build match rows
    match_rows = []
    for _, row in resolved_df.iterrows():
        base = row["gencode_id_raw"]
        if base in ensg_to_gtex:
            for (gencode_id, gtex_sym, gtex_type, entrez, gv, ds) in ensg_to_gtex[base]:
                match_rows.append({
                    "entity_key":     row["entity_key"],
                    "symbol":         row["symbol"],
                    "category":       row["category"],
                    "ncbi_gene_id":   row["ncbi_gene_id"],
                    "matched":        True,
                    "match_method":   row["match_method"],
                    "gencode_id":     gencode_id,
                    "gencode_version": gv,
                    "dataset_id":     ds,
                    "gtex_gene_symbol": gtex_sym,
                    "gtex_gene_type":   gtex_type,
                    "gtex_entrez_id":   entrez,
                })
        # entities that resolve via RNAcentral but whose ENSG is unknown to GTEx
        # are simply omitted; their existing matched=False row remains in the file.

    new_matches = pd.DataFrame(match_rows)

    # fetch expression only for truly new genes not in any existing pickle
    new_gencode_ids = set(new_matches["gencode_id"]) - set(existing_gene_expr["gencodeId"])
    new_gene_expr_frames, new_tx_expr_frames = [], []
    if new_gencode_ids:
        for dataset_id, grp in new_matches[new_matches["gencode_id"].isin(new_gencode_ids)].groupby("dataset_id"):
            ids = grp["gencode_id"].unique().tolist()
            print(f"  Fetching gene expression for {len(ids):,} new genes from {dataset_id} …")
            new_gene_expr_frames.append(pd.DataFrame(gtex.median_gene_expression(ids, dataset_id=dataset_id)))
            print(f"  Fetching transcript expression for {len(ids):,} new genes from {dataset_id} …")
            new_tx_expr_frames.append(pd.DataFrame(gtex.median_transcript_expression(ids, dataset_id=dataset_id)))

    new_gene_expr = pd.concat(new_gene_expr_frames, ignore_index=True) if new_gene_expr_frames else pd.DataFrame()
    new_tx_expr   = pd.concat(new_tx_expr_frames,   ignore_index=True) if new_tx_expr_frames   else pd.DataFrame()
    return new_matches, new_gene_expr, new_tx_expr


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # ── 0. Build DuckDB if needed ─────────────────────────────────────────
    if not DB_PATH.exists():
        for f in (ID_MAP_GZ, MD5_GZ):
            if not f.exists():
                raise FileNotFoundError(
                    f"{f} not found.\n"
                    "Run:\n"
                    "  rsync -avz rsync.ebi.ac.uk::pub/databases/RNAcentral/releases/25.0/id_mapping/ "
                    "rnacentral-etl/current_release/id_mapping/\n"
                    "  rsync -avz rsync.ebi.ac.uk::pub/databases/RNAcentral/releases/25.0/md5/ "
                    "rnacentral-etl/current_release/md5/"
                )
        _build_db()
    else:
        print(f"Using existing DuckDB at {DB_PATH}")

    # ── 1. Load existing data ─────────────────────────────────────────────
    print("Loading entity table and existing GTEx matches …")
    ents    = pd.read_pickle(ENTITIES_PKL)
    matches = pd.read_pickle(MATCHES_PKL)
    gene_expr = pd.read_pickle(GENE_EXPR_PKL)
    tx_expr   = pd.read_pickle(TX_EXPR_PKL)

    unmatched = load_unmatched(ents, matches)
    print(f"  Unmatched entities: {len(unmatched):,}")
    print(unmatched["category"].value_counts().to_string())

    # ── 2. Three resolution paths ─────────────────────────────────────────
    con = _open_db(read_only=True)

    print("\nPath A: URS → ENSEMBL_GENCODE …")
    res_a = _path_a(con, unmatched)
    resolved_a_keys = set(res_a["entity_key"]) if not res_a.empty else set()
    still_b = unmatched[~unmatched["entity_key"].isin(resolved_a_keys)]
    print(f"  {res_a['entity_key'].nunique() if not res_a.empty else 0:,} entities resolved "
          f"({len(res_a) if not res_a.empty else 0:,} entity-locus pairs)")

    print("Path B: miRBase MIMAT → URS → ENSEMBL_GENCODE …")
    res_b = _path_b(con, still_b)
    resolved_b_keys = set(res_b["entity_key"]) if not res_b.empty else set()
    still_c = still_b[~still_b["entity_key"].isin(resolved_b_keys)]
    print(f"  {res_b['entity_key'].nunique() if not res_b.empty else 0:,} entities resolved "
          f"({len(res_b) if not res_b.empty else 0:,} entity-locus pairs)")

    print("Path C: sequence MD5 → md5_map → ENSEMBL_GENCODE …")
    res_c = _path_c(con, still_c)
    print(f"  {res_c['entity_key'].nunique() if not res_c.empty else 0:,} entities resolved "
          f"({len(res_c) if not res_c.empty else 0:,} entity-locus pairs)")

    con.close()

    # ── 3. Combine, deduplicate across paths ──────────────────────────────
    all_resolved = pd.concat(
        [df for df in (res_a, res_b, res_c) if not df.empty],
        ignore_index=True
    )
    if all_resolved.empty:
        print("\nNo new entities resolved. Nothing to update.")
        return

    # keep only the first match method per (entity_key, gencode_id_raw) pair
    all_resolved = all_resolved.drop_duplicates(subset=["entity_key", "gencode_id_raw"])

    n_ents_resolved = all_resolved["entity_key"].nunique()
    print(f"\nTotal: {n_ents_resolved:,} new entities have >=1 RNAcentral-based ENSG candidate")
    print("By category:")
    print(all_resolved.groupby("category")["entity_key"].nunique().sort_values(ascending=False).to_string())

    # ── 4. GTEx verification + expression fetch ───────────────────────────
    print("\nVerifying against GTEx API and fetching expression …")
    new_matches, new_gene_expr, new_tx_expr = verify_and_fetch(all_resolved, gene_expr, tx_expr)

    if new_matches.empty:
        print("No RNAcentral candidates matched any GTEx record.")
        return

    n_gtex = new_matches["entity_key"].nunique()
    print(f"\n{n_gtex:,}/{n_ents_resolved:,} entities confirmed in GTEx "
          f"({len(new_matches):,} entity-locus pairs)")

    # ── 5. Append to existing pickles ─────────────────────────────────────
    # Remove existing matched=False rows for entities that are now matched
    now_matched_keys = set(new_matches["entity_key"])
    matches_updated = matches[
        ~((matches["entity_key"].isin(now_matched_keys)) & (~matches["matched"]))
    ]
    matches_updated = pd.concat([matches_updated, new_matches], ignore_index=True)
    matches_updated.to_pickle(MATCHES_PKL)
    print(f"Updated {MATCHES_PKL.name}: {len(matches_updated):,} rows "
          f"(was {len(matches):,})")

    if not new_gene_expr.empty:
        gene_expr_updated = pd.concat([gene_expr, new_gene_expr], ignore_index=True)
        gene_expr_updated.to_pickle(GENE_EXPR_PKL)
        print(f"Updated {GENE_EXPR_PKL.name}: {len(gene_expr_updated):,} rows "
              f"(+{len(new_gene_expr):,})")

    if not new_tx_expr.empty:
        tx_expr_updated = pd.concat([tx_expr, new_tx_expr], ignore_index=True)
        tx_expr_updated.to_pickle(TX_EXPR_PKL)
        print(f"Updated {TX_EXPR_PKL.name}: {len(tx_expr_updated):,} rows "
              f"(+{len(new_tx_expr):,})")

    # ── 6. Summary ────────────────────────────────────────────────────────
    print("\n=== Match rate after step 8 ===")
    all_matches = pd.read_pickle(MATCHES_PKL)
    per_entity  = all_matches.groupby("entity_key").agg(
        category=("category", "first"), matched=("matched", "any"))
    rate = per_entity.groupby("category")["matched"].agg(n_total="count", n_matched="sum")
    rate["match_rate"] = rate["n_matched"] / rate["n_total"]
    print(rate.sort_values("n_total", ascending=False).to_string())
    overall_matched = per_entity["matched"].sum()
    overall_total   = len(per_entity)
    print(f"\nOverall: {overall_matched:,}/{overall_total:,} ({overall_matched/overall_total:.1%})")
    print("\nRe-run steps 05 → 03 → 06 → 07 → 04 to propagate the new matches.")


if __name__ == "__main__":
    main()
