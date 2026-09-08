"""
Step 5: narrow GTEx's per-gene transcript list down to the specific
transcript(s) RNAInterSEQ actually has on record for each ncRNA entity.

Steps 2-4 matched entities to GENCODE *genes* and reported expression for
every annotated transcript of that gene -- but RNAInterSEQ's own sequence
for an entity (sourced from miRBase, NCBI RefSeq, or RNAcentral, almost
never from Ensembl directly -- see id_resolvers.py) corresponds to one
specific isoform, not the gene's whole transcript set. A lncRNA host gene
can carry dozens of GENCODE transcripts; reporting all of them as "this
ncRNA's expression" overstates what we actually know.

There is no shared identifier to join on across these databases, so the
join key has to be the sequence itself: fetch the spliced cDNA sequence of
every GENCODE transcript belonging to a matched gene (ensembl_api.py), and
test it for containment against the entity's own known sequence(s) (exact
equality, or substring containment in either direction -- handles a mature
miRNA sequence sitting inside its precursor transcript, or an isoform
boundary differing slightly between database releases).
"""
import pickle
from pathlib import Path

import pandas as pd

import ensembl_api as ens

ROOT = Path(__file__).parent
DATA = ROOT / "data"

ENTITIES_PKL = DATA / "ncrna_entities_human.pkl"
MATCHES_PKL = DATA / "gtex_gene_matches.pkl"
TX_EXPR_PKL = DATA / "gtex_transcript_expression.pkl"
TX_SEQ_PKL = DATA / "ensembl_transcript_sequences.pkl"
LOCUS_TX_MATCH_PKL = DATA / "ncrna_locus_transcript_matches.pkl"
TX_EXPR_FILTERED_PKL = DATA / "gtex_transcript_expression_sequence_confirmed.pkl"


def to_dna(seq):
    """RNAInter stores ncRNA sequences as RNA (U); Ensembl cDNA sequences use T."""
    return seq.strip().upper().replace("U", "T")


def sequence_matches(known_dna_seqs, candidate_seq):
    if not candidate_seq:
        return False
    for known in known_dna_seqs:
        if known == candidate_seq or known in candidate_seq or candidate_seq in known:
            return True
    return False


def main():
    ents = pd.read_pickle(ENTITIES_PKL)
    matches = pd.read_pickle(MATCHES_PKL)
    tx_expr = pd.read_pickle(TX_EXPR_PKL)

    matched_loci = matches[matches["matched"]].copy()
    known_seqs = ents.set_index("entity_key")["known_sequences"].apply(
        lambda seqs: [to_dna(s) for s in seqs])

    gene_to_transcripts = tx_expr.groupby("gencodeId")["transcriptId"].unique()
    all_transcript_ids = tx_expr["transcriptId"].unique().tolist()
    print(f"Fetching Ensembl cDNA sequences for {len(all_transcript_ids):,} "
          f"unique GENCODE transcripts ...")
    tx_seqs = ens.fetch_transcript_sequences(all_transcript_ids)
    n_found = sum(1 for v in tx_seqs.values() if v)
    print(f"  retrieved {n_found:,}/{len(tx_seqs):,} sequences "
          f"({len(tx_seqs) - n_found:,} not found at Ensembl, e.g. retired IDs)")

    tx_seq_df = pd.Series(tx_seqs, name="sequence").rename_axis("transcript_id_unversioned").reset_index()
    tx_seq_df.to_pickle(TX_SEQ_PKL)

    rows = []
    for locus in matched_loci.itertuples(index=False):
        entity_seqs = known_seqs.get(locus.entity_key, [])
        candidates = gene_to_transcripts.get(locus.gencode_id, [])
        for tid in candidates:
            seq = tx_seqs.get(tid.split(".")[0])
            confirmed = bool(entity_seqs) and sequence_matches(entity_seqs, seq)
            rows.append({
                "entity_key": locus.entity_key, "symbol": locus.symbol,
                "category": locus.category, "gencode_id": locus.gencode_id,
                "transcript_id": tid, "sequence_confirmed": confirmed,
            })
    locus_tx = pd.DataFrame(rows)
    locus_tx.to_pickle(LOCUS_TX_MATCH_PKL)

    n_loci = len(matched_loci)
    confirmed_per_locus = locus_tx.groupby(["entity_key", "gencode_id"])["sequence_confirmed"].any()
    n_loci_confirmed = int(confirmed_per_locus.sum())
    print(f"\n{n_loci_confirmed:,}/{n_loci:,} matched entity-locus pairs ({n_loci_confirmed/n_loci:.1%}) "
          f"have >=1 sequence-confirmed GENCODE transcript.")

    print("\nSequence-confirmation rate by category (entity-locus level):")
    rate = locus_tx.groupby(["entity_key", "gencode_id", "category"])["sequence_confirmed"].any() \
        .reset_index().groupby("category")["sequence_confirmed"].agg(n="count", confirmed="sum")
    rate["rate"] = rate["confirmed"] / rate["n"]
    print(rate.sort_values("n", ascending=False).to_string())

    tx_per_gene_total = locus_tx.groupby("gencode_id")["transcript_id"].nunique()
    tx_per_gene_confirmed = (locus_tx[locus_tx["sequence_confirmed"]]
                              .groupby("gencode_id")["transcript_id"].nunique())
    narrowed = (tx_per_gene_confirmed / tx_per_gene_total.reindex(tx_per_gene_confirmed.index)).dropna()
    print(f"\nFor genes with >=1 confirmed transcript, the confirmed set is on average "
          f"{narrowed.mean():.1%} of that gene's total annotated transcripts "
          f"(median {narrowed.median():.1%}) -- i.e. sequence confirmation meaningfully "
          f"narrows multi-isoform genes down to the specific isoform RNAInterSEQ has.")

    confirmed_index = pd.MultiIndex.from_frame(
        locus_tx.loc[locus_tx["sequence_confirmed"], ["gencode_id", "transcript_id"]]
        .drop_duplicates()
        .rename(columns={"gencode_id": "gencodeId", "transcript_id": "transcriptId"})
    )
    tx_pairs = pd.MultiIndex.from_frame(tx_expr[["gencodeId", "transcriptId"]])
    tx_expr_filtered = tx_expr[tx_pairs.isin(confirmed_index)]
    tx_expr_filtered.to_pickle(TX_EXPR_FILTERED_PKL)
    print(f"\n{len(tx_expr_filtered):,}/{len(tx_expr):,} (transcript, tissue) expression rows "
          f"survive sequence filtering -> {TX_EXPR_FILTERED_PKL}")


if __name__ == "__main__":
    main()
