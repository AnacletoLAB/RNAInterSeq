"""
Step 6: Annotate the human ncRNA-ncRNA interaction table with GTEx expression metadata.

For every human-human interaction row, adds:

Per-entity columns (ent1_* and ent2_*):
  gtex_matched       -- bool: entity resolved to >=1 GENCODE gene in GTEx
  seq_confirmed      -- bool: >=1 transcript sequence-confirmed against RNAInterSEQ
  confirmed_transcripts -- dict {gencode_id: [transcript_ids]} (sequence-confirmed only;
                           empty dict if no confirmed match)
  tau                -- Yanai tissue-specificity index of the best-expressed confirmed
                        locus (NaN if not detected in any tissue)
  peak_tissue        -- GTEx tissue with highest median TPM across confirmed transcripts
  peak_tpm           -- median TPM at peak tissue
  n_expressed_tissues -- number of 54 GTEx tissues where aggregated max-TPM > 0
  tissue_tpm         -- dict {tissue_id: max_median_tpm} aggregated across all confirmed
                        transcripts from all matched loci (54 entries; 0.0 where undetected)

Per-entity system columns (ent1_* and ent2_*):
  tissue_systems     -- sorted list of biological systems (CNS, GI, …) where
                        the entity is expressed (max-TPM > 0 in ≥1 tissue of
                        that system); empty list if no confirmed expression

Derived interaction columns:
  both_gtex_matched  -- both partners have a GTEx match
  both_seq_confirmed -- both partners have sequence-confirmed transcripts
  shared_expressed_tissues -- sorted list of tissues where both partners have TPM > 0
  n_shared_expressed_tissues -- len of the above
  shared_tissue_systems    -- sorted list of the systems to which the shared
                              expressed tissues belong (non-empty only when the
                              two partners are detected in the same tissue)
  n_shared_tissue_systems  -- len of the above
  shared_systems           -- sorted list of systems in which BOTH partners are
                              detected, i.e. the intersection of ent1_tissue_systems
                              and ent2_tissue_systems (the partners may be detected
                              in different tissues of the same system); a superset
                              of shared_tissue_systems, sequence-confirmed pairs only
  n_shared_systems         -- len of the above

Output: data/human_ncrna_interactions_gtex_annotated.pkl
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent
DATA = ROOT / "data"

# Mapping of the 54 GTEx tissue IDs to 8 biological systems.
TISSUE_SYSTEM = {
    # CNS (13)
    "Brain_Amygdala":                           "CNS",
    "Brain_Anterior_cingulate_cortex_BA24":     "CNS",
    "Brain_Caudate_basal_ganglia":              "CNS",
    "Brain_Cerebellar_Hemisphere":              "CNS",
    "Brain_Cerebellum":                         "CNS",
    "Brain_Cortex":                             "CNS",
    "Brain_Frontal_Cortex_BA9":                 "CNS",
    "Brain_Hippocampus":                        "CNS",
    "Brain_Hypothalamus":                       "CNS",
    "Brain_Nucleus_accumbens_basal_ganglia":    "CNS",
    "Brain_Putamen_basal_ganglia":              "CNS",
    "Brain_Spinal_cord_cervical_c-1":           "CNS",
    "Brain_Substantia_nigra":                   "CNS",
    # Peripheral nervous (2)
    "Nerve_Tibial":                             "Peripheral nervous",
    "Pituitary":                                "Peripheral nervous",
    # Cardiovascular (5)
    "Artery_Aorta":                             "Cardiovascular",
    "Artery_Coronary":                          "Cardiovascular",
    "Artery_Tibial":                            "Cardiovascular",
    "Heart_Atrial_Appendage":                   "Cardiovascular",
    "Heart_Left_Ventricle":                     "Cardiovascular",
    # GI (10)
    "Colon_Sigmoid":                            "GI",
    "Colon_Transverse":                         "GI",
    "Esophagus_Gastroesophageal_Junction":      "GI",
    "Esophagus_Mucosa":                         "GI",
    "Esophagus_Muscularis":                     "GI",
    "Liver":                                    "GI",
    "Minor_Salivary_Gland":                     "GI",
    "Pancreas":                                 "GI",
    "Small_Intestine_Terminal_Ileum":           "GI",
    "Stomach":                                  "GI",
    # Reproductive (9)
    "Breast_Mammary_Tissue":                    "Reproductive",
    "Cervix_Ectocervix":                        "Reproductive",
    "Cervix_Endocervix":                        "Reproductive",
    "Fallopian_Tube":                           "Reproductive",
    "Ovary":                                    "Reproductive",
    "Prostate":                                 "Reproductive",
    "Testis":                                   "Reproductive",
    "Uterus":                                   "Reproductive",
    "Vagina":                                   "Reproductive",
    # Immune (3)
    "Cells_EBV-transformed_lymphocytes":        "Immune",
    "Spleen":                                   "Immune",
    "Whole_Blood":                              "Immune",
    # Muscle (1)
    "Muscle_Skeletal":                          "Muscle",
    # Other (11)
    "Adipose_Subcutaneous":                     "Other",
    "Adipose_Visceral_Omentum":                 "Other",
    "Adrenal_Gland":                            "Other",
    "Bladder":                                  "Other",
    "Cells_Cultured_fibroblasts":               "Other",
    "Kidney_Cortex":                            "Other",
    "Kidney_Medulla":                           "Other",
    "Lung":                                     "Other",
    "Skin_Not_Sun_Exposed_Suprapubic":          "Other",
    "Skin_Sun_Exposed_Lower_leg":               "Other",
    "Thyroid":                                  "Other",
}

INTERACTIONS_PKL = ROOT.parent / "data" / "interactions_with_sequences_RNA_RNA_ncrna_ncrna_19-05-2026.pkl"
SUMMARY_PKL      = DATA / "ncrna_entity_locus_summary.pkl"
LOCUS_TX_PKL     = DATA / "ncrna_locus_transcript_matches.pkl"
TX_EXPR_CONF_PKL = DATA / "gtex_transcript_expression_sequence_confirmed.pkl"
OUT_PKL          = DATA / "human_ncrna_interactions_gtex_annotated.pkl"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _primary_locus_stats(group):
    """For a multi-locus entity, keep the row with the highest peak_median_tpm."""
    confirmed_rows = group[group["sequence_confirmed"]]
    best = confirmed_rows if not confirmed_rows.empty else group
    idx = best["peak_median_tpm"].fillna(0).idxmax()
    return best.loc[idx, ["tau", "peak_tissue", "peak_median_tpm"]]


def build_entity_lookup(summary, locus_tx, tx_expr_filtered):
    """Returns a DataFrame indexed by entity_key with one row per entity."""
    print("  Building confirmed-transcript map …")
    confirmed = locus_tx[locus_tx["sequence_confirmed"]].copy()

    # {entity_key: {gencode_id: [transcript_id, ...]}}
    entity_transcripts = (
        confirmed.groupby("entity_key")
        .apply(lambda g: g.groupby("gencode_id")["transcript_id"].apply(list).to_dict(),
               include_groups=False)
        .rename("confirmed_transcripts")
    )

    print("  Building per-entity tissue TPM vectors …")
    # join confirmed transcripts → expression rows
    tx_to_entity = (
        confirmed[["entity_key", "gencode_id", "transcript_id"]]
        .drop_duplicates()
        .rename(columns={"gencode_id": "gencodeId", "transcript_id": "transcriptId"})
    )
    tx_with_entity = tx_expr_filtered.merge(tx_to_entity, on=["gencodeId", "transcriptId"])

    # max TPM per entity per tissue (across all confirmed transcripts / loci)
    entity_tissue = (
        tx_with_entity
        .groupby(["entity_key", "tissueSiteDetailId"])["median"]
        .max()
        .unstack(fill_value=0.0)
    )
    all_tissues = sorted(tx_expr_filtered["tissueSiteDetailId"].unique())
    entity_tissue = entity_tissue.reindex(columns=all_tissues, fill_value=0.0)
    entity_tissue_tpm = entity_tissue.apply(lambda row: row.to_dict(), axis=1).rename("tissue_tpm")
    entity_n_expressed = (entity_tissue > 0).sum(axis=1).rename("n_expressed_tissues")

    print("  Computing primary-locus stats (tau, peak_tissue, peak_tpm) …")
    locus_stats = (
        summary[summary["matched"]]
        .groupby("entity_key")
        .apply(_primary_locus_stats, include_groups=False)
        .rename(columns={"peak_median_tpm": "peak_tpm"})
    )

    # GTEx match / sequence-confirmed flags
    matched_entities = set(summary[summary["matched"]]["entity_key"])
    confirmed_entities = set(confirmed["entity_key"])

    print("  Assembling entity lookup table …")
    lookup = pd.DataFrame({"entity_key": list(matched_entities | confirmed_entities)})
    lookup["gtex_matched"]  = lookup["entity_key"].isin(matched_entities)
    lookup["seq_confirmed"] = lookup["entity_key"].isin(confirmed_entities)
    lookup = lookup.join(entity_transcripts, on="entity_key")
    lookup = lookup.join(locus_stats,        on="entity_key")
    lookup = lookup.join(entity_tissue_tpm,  on="entity_key")
    lookup = lookup.join(entity_n_expressed, on="entity_key")

    lookup["confirmed_transcripts"] = lookup["confirmed_transcripts"].apply(
        lambda x: x if isinstance(x, dict) else {})
    lookup["tissue_tpm"] = lookup["tissue_tpm"].apply(
        lambda x: x if isinstance(x, dict) else {})
    lookup["n_expressed_tissues"] = lookup["n_expressed_tissues"].fillna(0).astype(int)

    return lookup.set_index("entity_key")


def annotate(interactions, lookup):
    """Left-join lookup columns onto interactions for both partners."""
    cols = ["gtex_matched", "seq_confirmed", "confirmed_transcripts",
            "tau", "peak_tissue", "peak_tpm", "n_expressed_tissues", "tissue_tpm"]

    result = interactions.copy()
    for prefix, key_col in [("ent1", "key1"), ("ent2", "key2")]:
        side = lookup.reindex(result[key_col])[cols]
        side.index = result.index
        side.columns = [f"{prefix}_{c}" for c in cols]
        result = pd.concat([result, side], axis=1)

    # fill bool cols for unmatched entities (where avoids FutureWarning from fillna on object dtype)
    for prefix in ("ent1", "ent2"):
        for col in (f"{prefix}_gtex_matched", f"{prefix}_seq_confirmed"):
            s = result[col]
            result[col] = s.where(s.notna(), other=False).astype(bool)
        result[f"{prefix}_confirmed_transcripts"] = result[f"{prefix}_confirmed_transcripts"].apply(
            lambda x: x if isinstance(x, dict) else {})
        result[f"{prefix}_tissue_tpm"] = result[f"{prefix}_tissue_tpm"].apply(
            lambda x: x if isinstance(x, dict) else {})
        result[f"{prefix}_n_expressed_tissues"] = result[f"{prefix}_n_expressed_tissues"].fillna(0).astype(int)

    result["both_gtex_matched"]  = result["ent1_gtex_matched"]  & result["ent2_gtex_matched"]
    result["both_seq_confirmed"] = result["ent1_seq_confirmed"] & result["ent2_seq_confirmed"]

    print("  Computing shared expressed tissues …")
    def _shared(row):
        t1 = {t for t, v in row["ent1_tissue_tpm"].items() if v > 0}
        t2 = {t for t, v in row["ent2_tissue_tpm"].items() if v > 0}
        return sorted(t1 & t2)

    both_conf = result["both_seq_confirmed"]
    shared = result[both_conf].apply(_shared, axis=1)
    result.loc[both_conf, "shared_expressed_tissues"] = shared
    result["shared_expressed_tissues"] = result["shared_expressed_tissues"].apply(
        lambda x: x if isinstance(x, list) else [])
    result["n_shared_expressed_tissues"] = result["shared_expressed_tissues"].str.len()

    # ── tissue system grouping ────────────────────────────────────────────────
    def _to_systems(tissues):
        return sorted(set(TISSUE_SYSTEM.get(t, "Other") for t in tissues))

    for prefix in ("ent1", "ent2"):
        expressed = result[f"{prefix}_tissue_tpm"].apply(
            lambda d: [t for t, v in d.items() if v > 0])
        result[f"{prefix}_tissue_systems"] = expressed.apply(_to_systems)

    result["shared_tissue_systems"]   = result["shared_expressed_tissues"].apply(_to_systems)
    result["n_shared_tissue_systems"] = result["shared_tissue_systems"].str.len()
    # system-level sharing proper: both partners detected somewhere in the same system
    result["shared_systems"] = [
        sorted(set(a) & set(b)) if conf else []
        for a, b, conf in zip(result["ent1_tissue_systems"], result["ent2_tissue_systems"], both_conf)]
    result["n_shared_systems"] = result["shared_systems"].str.len()

    return result


def main():
    print(f"Loading interaction table …")
    df = pd.read_pickle(INTERACTIONS_PKL)
    human = df[(df["Species1"] == "Homo sapiens") & (df["Species2"] == "Homo sapiens")].copy()
    print(f"  Human-human interactions: {len(human):,}")

    print("Loading GTEx artefacts …")
    summary      = pd.read_pickle(SUMMARY_PKL)
    locus_tx     = pd.read_pickle(LOCUS_TX_PKL)
    tx_expr_conf = pd.read_pickle(TX_EXPR_CONF_PKL)

    print("Building entity lookup …")
    lookup = build_entity_lookup(summary, locus_tx, tx_expr_conf)
    print(f"  Lookup covers {len(lookup):,} entities")

    print("Annotating interactions …")
    annotated = annotate(human, lookup)

    print("\n=== Annotation coverage ===")
    n = len(annotated)
    n_both = annotated["both_gtex_matched"].sum()
    n_both_conf = annotated["both_seq_confirmed"].sum()
    n_shared_any = (annotated["n_shared_expressed_tissues"] > 0).sum()
    print(f"  Total human-human interactions: {n:,}")
    print(f"  Both partners GTEx-matched:     {n_both:,} ({n_both/n:.1%})")
    print(f"  Both seq-confirmed:             {n_both_conf:,} ({n_both_conf/n:.1%})")
    print(f"  >=1 shared expressed tissue:    {n_shared_any:,} ({n_shared_any/n:.1%})")
    n_shared_sys = (annotated["n_shared_systems"] > 0).sum()
    print(f"  >=1 shared biological system:   {n_shared_sys:,} ({n_shared_sys/n:.1%})")
    med_shared = annotated.loc[annotated["both_seq_confirmed"], "n_shared_expressed_tissues"].median()
    print(f"  Median shared tissues (both-confirmed pairs): {med_shared:.0f}")

    annotated.to_pickle(OUT_PKL)
    print(f"\nSaved {OUT_PKL}")
    print(f"Final DataFrame shape: {annotated.shape}")
    print("New columns added:", [c for c in annotated.columns if c not in df.columns])


if __name__ == "__main__":
    main()
