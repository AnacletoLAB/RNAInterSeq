"""RNA_entities_plot: resolution of RNA-entities to transcript sequences.

Replaces the ASCII-art rendering with the same visual language as
03_flow_diagram.py (rounded boxes, cream input, blue ncRNA, green mRNA, orange
result). The mRNA count is read from the newest processed mRNA table, so the
figure follows the pipeline: 197,490 (fallback symbol mapping) -> 197,491
(symmetric symbol mapping, Sep 2026) -> the value after the re-fetch of the
RefSeq accessions lost to failed efetch batches (2c_refetch_refseq_failures.py).

Fixed counts, verified against the pipeline tables:
  data/nc_rna_entities_to_be_processed.pkl   76,961 rows
  data/m_rna_entities_to_be_processed.pkl   213,723 rows
  data/nc_rna_entities_processed.pkl         40,911 rows, all with >=1 sequence
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

CREAM="#FFF2CC"; CREAM_E="#D6B656"
ORANGE="#FFE6CC"; ORANGE_E="#D79B00"
GREEN="#D5E8D4"; GREEN_E="#82B366"
BLUE="#DAE8FC"; BLUE_E="#6C8EBF"

import glob, json, pandas as pd
NC_ALL, MR_ALL = 76_961, 213_723            # data/*_rna_entities_to_be_processed.pkl
NC_RES = 40_911                             # data/nc_rna_entities_processed.pkl
# mRNA: the newest processed table (2c_refetch_refseq_failures.py rebuilt it with the re-fetched accessions)
MR_TABLE = sorted(glob.glob("data/m_rna_entities_processed_union*.pkl"))[-1]
MR_RES = len(pd.read_pickle(MR_TABLE).index)
ALL = NC_ALL + MR_ALL
RES = NC_RES + MR_RES
print(f"{MR_TABLE}: {MR_RES:,} sequence-resolved mRNA-entities -> {RES:,} overall")

fig, ax = plt.subplots(figsize=(9.2, 8.2))
ax.set_xlim(0, 10); ax.set_ylim(0, 9.6); ax.axis("off")

def box(x, y, w, h, title, value, fc, ec, fs=13, vfs=15):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        facecolor=fc, edgecolor=ec, linewidth=1.6))
    ax.text(x, y+h*0.16, title, ha="center", va="center", fontsize=fs)
    ax.text(x, y-h*0.22, value, ha="center", va="center",
            fontsize=vfs, fontweight="bold")

def arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
        arrowstyle="-|>", mutation_scale=17, linewidth=1.5, color="#444444"))

def note(x, y, s, ha="center", fs=9.5):
    ax.text(x, y, s, ha=ha, va="center", fontsize=fs, color="#555555")

H = 1.15
LX, RX = 2.55, 7.45          # ncRNA branch / mRNA branch

# --- non-redundant RNA-entity table -------------------------------------
box(5, 8.85, 5.6, H, "Non-redundant RNA-entities", f"{ALL:,}", CREAM, CREAM_E)

# split into the two biotype branches
ax.plot([5, 5], [8.27, 7.85], color="#444444", linewidth=1.5)
ax.plot([LX, RX], [7.85, 7.85], color="#444444", linewidth=1.5)
arrow(LX, 7.85, LX, 7.32)
arrow(RX, 7.85, RX, 7.32)

box(LX, 6.75, 4.3, H, "ncRNA-entities", f"{NC_ALL:,}", BLUE, BLUE_E)
box(RX, 6.75, 4.3, H, "mRNA-entities", f"{MR_ALL:,}", GREEN, GREEN_E)

# --- the two sequence-mapping pipelines ---------------------------------
arrow(LX, 6.17, LX, 4.53)
note(LX-0.12, 5.35, "five-step identifier\nresolution via\nRNAcentral", ha="right")
arrow(RX, 6.17, RX, 4.53)
note(RX+0.12, 5.35, "three-step mapping via\nEnsembl, NCBI Gene\nand gene symbol", ha="left")

box(LX, 3.95, 4.3, H, "sequence-resolved", f"{NC_RES:,}", BLUE, BLUE_E)
note(LX, 3.08, f"{100*NC_RES/NC_ALL:.1f} % of ncRNA-entities\n"
               f"({NC_ALL-NC_RES:,} unresolved)")
box(RX, 3.95, 4.3, H, "sequence-resolved", f"{MR_RES:,}", GREEN, GREEN_E)
note(RX, 3.08, f"{100*MR_RES/MR_ALL:.1f} % of mRNA-entities\n"
               f"({MR_ALL-MR_RES:,} unresolved)")

# --- merge back into the sequence-resolved entity table -----------------
ax.plot([LX, LX], [2.60, 2.28], color="#444444", linewidth=1.5)
ax.plot([RX, RX], [2.60, 2.28], color="#444444", linewidth=1.5)
ax.plot([LX, RX], [2.28, 2.28], color="#444444", linewidth=1.5)
arrow(5, 2.28, 5, 1.83)

box(5, 1.25, 5.6, H, "Sequence-resolved RNA-entities", f"{RES:,}", ORANGE, ORANGE_E)
note(5, 0.42, f"{100*RES/ALL:.1f} % of the non-redundant RNA-entities resolved "
              f"to at least one transcript sequence", fs=10.5)

plt.tight_layout()
import os, sys
OUT = sys.argv[1] if len(sys.argv) > 1 else "RNA_entities_plot.png"
os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
print("written:", OUT)
