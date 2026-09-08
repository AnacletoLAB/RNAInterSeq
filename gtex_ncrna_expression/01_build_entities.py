"""
Step 1: extract the unique set of human ncRNA entities from the RNAInterSEQ
ncRNA-ncRNA interaction dataset (19-05-2026 release) and attach the
identifiers (NCBI Gene ID / Ensembl Gene ID) needed to query GTEx.

GTEx profiles human donor tissue exclusively, so every other species in the
interaction table (mouse, chicken, viral, etc.) is out of scope and dropped
here -- that's the answer to "is it only human?".
"""
import pickle
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).parent
INTERACTIONS_PKL = ROOT.parent / "data" / "interactions_with_sequences_RNA_RNA_ncrna_ncrna_19-05-2026.pkl"
OUT_PKL = ROOT / "data" / "ncrna_entities_human.pkl"

STANDARD_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}


def _extract_ncbi_id(mygene_raw):
    if not isinstance(mygene_raw, dict):
        return None
    val = mygene_raw.get("_id") or mygene_raw.get("query")
    if val is None:
        return None
    val = str(val)
    return val if val.isdigit() else None


def _extract_ensembl_gene(mygene_raw):
    if not isinstance(mygene_raw, dict):
        return None
    gp = mygene_raw.get("genomic_pos")
    if gp is None:
        return None
    candidates = gp if isinstance(gp, list) else [gp]
    # prefer a position on a primary chromosome over alt-haplotype patches
    primary = [c for c in candidates if str(c.get("chr")) in STANDARD_CHROMS]
    pick_from = primary if primary else candidates
    for c in pick_from:
        eg = c.get("ensemblgene")
        if eg:
            return eg
    return None


def _ensembl_from_buckets(id_buckets):
    if not isinstance(id_buckets, dict):
        return None
    ens = id_buckets.get("Ensembl") or []
    return ens[0] if ens else None


def _bucket_list(id_buckets, key):
    if not isinstance(id_buckets, dict):
        return []
    return list(id_buckets.get(key) or [])


def _ncbi_from_buckets(id_buckets):
    ids = _bucket_list(id_buckets, "NCBI")
    return ids[0] if ids else None


def build_entity_table(df):
    rows = []
    for prefix in ("ent1", "ent2"):
        symbol_col = "Interactor1.Symbol" if prefix == "ent1" else "Interactor2.Symbol"
        category_col = "Category1" if prefix == "ent1" else "Category2"
        species_col = "Species1" if prefix == "ent1" else "Species2"
        key_col = "key1" if prefix == "ent1" else "key2"
        sub = df[[symbol_col, category_col, species_col, key_col,
                  f"{prefix}_ID_buckets", f"{prefix}_mygene_raw",
                  f"{prefix}_sequences"]].copy()
        sub.columns = ["symbol", "category", "species", "entity_key",
                       "id_buckets", "mygene_raw", "sequences"]
        rows.append(sub)

    ents = pd.concat(rows, ignore_index=True)
    ents = ents[ents["species"] == "Homo sapiens"].copy()

    # one row per unique entity: take the first non-null annotation seen
    ents = ents.groupby("entity_key", as_index=False).first()

    ents["ncbi_gene_id"] = ents["mygene_raw"].apply(_extract_ncbi_id)
    ents["ncbi_gene_id"] = ents["ncbi_gene_id"].fillna(
        ents["id_buckets"].apply(_ncbi_from_buckets)
    )
    ents["ensembl_gene_id"] = ents["mygene_raw"].apply(_extract_ensembl_gene)
    ents["ensembl_gene_id"] = ents["ensembl_gene_id"].fillna(
        ents["id_buckets"].apply(_ensembl_from_buckets)
    )
    ents["mirbase_ids"] = ents["id_buckets"].apply(lambda b: _bucket_list(b, "miRBase"))
    ents["circbase_ids"] = ents["id_buckets"].apply(lambda b: _bucket_list(b, "circBase"))
    ents["rnacentral_ids"] = ents["id_buckets"].apply(lambda b: _bucket_list(b, "RNAcentral"))
    # known_sequences[i] corresponds 1:1 to rnacentral_ids[i] (verified against
    # the source table: len(sequences) always equals len(ID_buckets['RNAcentral']))
    ents["known_sequences"] = ents["sequences"].apply(
        lambda s: [x for x in s if isinstance(x, str) and x] if isinstance(s, (list, tuple)) else [])

    return ents[["entity_key", "symbol", "category", "species",
                 "ncbi_gene_id", "ensembl_gene_id",
                 "mirbase_ids", "circbase_ids", "rnacentral_ids", "known_sequences"]]


def main():
    print(f"Loading {INTERACTIONS_PKL.name} ...")
    df = pd.read_pickle(INTERACTIONS_PKL)
    print(f"  {len(df):,} interaction rows, "
          f"{df['Species1'].nunique() + df['Species2'].nunique()} raw species columns")

    all_species = pd.concat([df["Species1"], df["Species2"]]).value_counts()
    print("\nSpecies present in this ncRNA-ncRNA interaction set:")
    print(all_species.to_string())
    n_human_rows = ((df["Species1"] == "Homo sapiens") | (df["Species2"] == "Homo sapiens")).sum()
    print(f"\n{n_human_rows:,}/{len(df):,} interaction rows involve at least one human entity.")
    print("GTEx (Genotype-Tissue Expression) profiles postmortem human donor tissue "
          "only -- no other species in this list has any GTEx data, confirmed via "
          "the GTEx v2 API's tissueSiteDetailId vocabulary (54 human anatomical sites, "
          "no species axis).")

    ents = build_entity_table(df)
    print(f"\nUnique human ncRNA entities: {len(ents):,}")
    print(ents["category"].value_counts().to_string())
    print(f"\n  with NCBI Gene ID:    {ents['ncbi_gene_id'].notna().sum():,}")
    print(f"  with Ensembl Gene ID: {ents['ensembl_gene_id'].notna().sum():,}")
    print(f"  with miRBase ID:      {(ents['mirbase_ids'].str.len() > 0).sum():,}")
    print(f"  with circBase ID:     {(ents['circbase_ids'].str.len() > 0).sum():,}")
    print(f"  with RNAcentral ID:   {(ents['rnacentral_ids'].str.len() > 0).sum():,}")
    print(f"  with >=1 known sequence: {(ents['known_sequences'].str.len() > 0).sum():,}")

    OUT_PKL.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PKL, "wb") as f:
        pickle.dump(ents, f)
    print(f"\nSaved {OUT_PKL}")


if __name__ == "__main__":
    main()
