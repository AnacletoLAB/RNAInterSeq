"""
Local-resource resolvers that turn the non-Ensembl identifiers stored in the
RNAInterSEQ ID buckets (NCBI Gene ID, miRBase mature accession, circBase
accession) into a GTEx/GENCODE-queryable gene symbol or Ensembl Gene ID,
without round-tripping through a live mygene.info gene-by-gene call for
every entity. Two of the three sources are entirely local files already
shipped in this repo (`data/gene2ensembl`, `data/circdata/circBase_9606.fa`,
`data/miRNA.dat`); only the miRBase precursor-name -> official gene symbol
step needs a (batched, cached) mygene.info lookup, because miRBase precursor
names don't always match the HGNC/GENCODE symbol by a fixed string rule
(e.g. "let-7a-1" -> "MIRLET7A1", not "MIRLET-7A-1").

GenBank/tRFdb/"Other" ID buckets are not handled here: they are empty for
every human entity in the 19-05-2026 ncRNA-ncRNA dataset (checked directly
against the source table), so there is nothing to resolve.
"""
import json
import re
import subprocess
from pathlib import Path

import requests

ROOT = Path(__file__).parent
CACHE_DIR = ROOT / "cache"
DATA_DIR = ROOT.parent / "data"

GENE2ENSEMBL = DATA_DIR / "gene2ensembl"
CIRCBASE_FASTA = DATA_DIR / "circdata" / "circBase_9606.fa"
MIRBASE_DAT = DATA_DIR / "miRNA.dat"
HUMAN_TAXID = "9606"

_NCBI_TO_ENSEMBL_CACHE = CACHE_DIR / "ncbi_to_ensembl_human.json"
_MIRBASE_MAP_CACHE = CACHE_DIR / "mirbase_mature_to_precursor_human.json"
_PRECURSOR_RESOLUTION_CACHE = CACHE_DIR / "mirbase_precursor_to_gene.json"


# ---------------------------------------------------------------------------
# NCBI Gene ID -> Ensembl Gene ID (NCBI's own gene2ensembl mapping, human only)
# ---------------------------------------------------------------------------
def load_ncbi_to_ensembl_human():
    if _NCBI_TO_ENSEMBL_CACHE.exists():
        return json.loads(_NCBI_TO_ENSEMBL_CACHE.read_text())

    # gene2ensembl is ~1.6GB across all species; awk-filter to human rows
    # only instead of paying a full pandas read.
    out = {}
    proc = subprocess.run(
        ["awk", "-F\t", f"$1==\"{HUMAN_TAXID}\" {{print $2\"\\t\"$3}}", str(GENE2ENSEMBL)],
        capture_output=True, text=True, check=True,
    )
    for line in proc.stdout.splitlines():
        ncbi_id, ensembl_id = line.split("\t")
        out.setdefault(ncbi_id, ensembl_id)  # keep first mapping; ambiguous cases are rare (<0.1%)

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    _NCBI_TO_ENSEMBL_CACHE.write_text(json.dumps(out))
    return out


# ---------------------------------------------------------------------------
# circBase accession -> host gene symbol (parsed straight from FASTA headers:
# ">hsa_circ_0000650|NM_020210|SEMA4B")
# ---------------------------------------------------------------------------
def load_circbase_host_genes_human():
    out = {}
    with open(CIRCBASE_FASTA) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            parts = line[1:].strip().split("|")
            if len(parts) >= 3 and parts[2]:
                out[parts[0]] = parts[2]
    return out


# ---------------------------------------------------------------------------
# miRBase mature accession (MIMATxxxxxxx) -> candidate precursor name(s)
# A mature product can be produced by more than one genomic precursor locus
# (e.g. hsa-miR-101-3p comes from both hsa-mir-101-1 and hsa-mir-101-2).
# ---------------------------------------------------------------------------
_ID_RE = re.compile(r"^ID   (\S+)")
_ACC_RE = re.compile(r'/accession="(MIMAT\d+)"')


def load_mirbase_mature_to_precursor_human():
    if _MIRBASE_MAP_CACHE.exists():
        return json.loads(_MIRBASE_MAP_CACHE.read_text())

    out = {}
    precursor_name = None
    in_hsa_entry = False
    with open(MIRBASE_DAT, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ID   "):
                m = _ID_RE.match(line)
                full_id = m.group(1) if m else None
                in_hsa_entry = bool(full_id and full_id.startswith("hsa-"))
                precursor_name = full_id[4:] if in_hsa_entry else None  # strip "hsa-"
            elif in_hsa_entry and line.startswith("FT") and "/accession=" in line:
                m = _ACC_RE.search(line)
                if m:
                    out.setdefault(m.group(1), []).append(precursor_name)
            elif line.startswith("//"):
                in_hsa_entry = False
                precursor_name = None

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    _MIRBASE_MAP_CACHE.write_text(json.dumps(out))
    return out


# ---------------------------------------------------------------------------
# precursor name (e.g. "mir-101-1", "let-7a-1") -> official symbol / IDs,
# via a single batched mygene.info alias lookup (deterministic: aliases
# store the bare miRBase precursor name, so this is an exact match, not a
# fuzzy text search).
# ---------------------------------------------------------------------------
def resolve_precursor_symbols(precursor_names, batch_size=900):
    cache = {}
    if _PRECURSOR_RESOLUTION_CACHE.exists():
        cache = json.loads(_PRECURSOR_RESOLUTION_CACHE.read_text())

    todo = sorted(set(precursor_names) - cache.keys())
    for i in range(0, len(todo), batch_size):
        chunk = todo[i:i + batch_size]
        resp = requests.post(
            "https://mygene.info/v3/query",
            data={
                "q": ",".join(chunk),
                "scopes": "alias",
                "species": "human",
                "fields": "symbol,entrezgene,ensembl.gene",
            },
            timeout=60,
        )
        resp.raise_for_status()
        for rec in resp.json():
            ensembl = rec.get("ensembl")
            if isinstance(ensembl, list):
                ensembl = ensembl[0] if ensembl else {}
            cache[rec["query"]] = None if rec.get("notfound") else {
                "symbol": rec.get("symbol"),
                "entrezgene": rec.get("entrezgene"),
                "ensembl_gene": (ensembl or {}).get("gene"),
            }

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    _PRECURSOR_RESOLUTION_CACHE.write_text(json.dumps(cache))
    return {k: cache[k] for k in precursor_names if k in cache}
