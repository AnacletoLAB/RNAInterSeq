"""
Thin client for the GTEx Portal API v2 (https://gtexportal.org/api/v2).

GTEx (Genotype-Tissue Expression) profiles postmortem tissue samples from
human donors only -- there is no other-species data behind this API. This
client targets the gtex_v10 dataset (GENCODE v39 / GRCh38), the most
current GTEx release, and exposes the three endpoints needed to go from an
RNAInter ncRNA entity to a per-tissue expression profile:

  - reference/gene                  symbol/Ensembl ID  -> GENCODE gene record
  - expression/medianGeneExpression gencodeId          -> median TPM per tissue
  - expression/medianTranscriptExpression gencodeId    -> median TPM per
                                                            transcript per tissue

All calls are paginated, batched, retried on transient errors, and cached
to disk (one JSON file per request) so repeated runs don't re-hit the API.
"""
import hashlib
import json
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter, Retry

BASE_URL = "https://gtexportal.org/api/v2"
DEFAULT_DATASET = "gtex_v10"
DEFAULT_GENCODE_VERSION = "v39"
CACHE_DIR = Path(__file__).parent / "cache"

GTEX_TISSUES = [
    "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland",
    "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Bladder",
    "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus", "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra",
    "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts",
    "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts",
    "Cervix_Ectocervix", "Cervix_Endocervix", "Colon_Sigmoid",
    "Colon_Transverse", "Esophagus_Gastroesophageal_Junction",
    "Esophagus_Mucosa", "Esophagus_Muscularis", "Fallopian_Tube",
    "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex",
    "Kidney_Medulla", "Liver", "Lung", "Minor_Salivary_Gland",
    "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary",
    "Prostate", "Skin_Not_Sun_Exposed_Suprapubic",
    "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum",
    "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina",
    "Whole_Blood",
]


def _session():
    s = requests.Session()
    retries = Retry(total=5, backoff_factor=0.5,
                     status_forcelist=[429, 500, 502, 503, 504])
    s.mount("https://", HTTPAdapter(max_retries=retries))
    return s


_SESSION = _session()


def _cache_path(endpoint, params):
    key = json.dumps(params, sort_keys=True, default=str)
    digest = hashlib.sha1(f"{endpoint}|{key}".encode()).hexdigest()
    return CACHE_DIR / f"{digest}.json"


def _get(endpoint, params, use_cache=True, sleep=0.1):
    cache_file = _cache_path(endpoint, params)
    if use_cache and cache_file.exists():
        return json.loads(cache_file.read_text())
    resp = _SESSION.get(f"{BASE_URL}/{endpoint}", params=params, timeout=60)
    resp.raise_for_status()
    data = resp.json()
    if use_cache:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(data))
    if sleep:
        time.sleep(sleep)
    return data


def _paginate(endpoint, params, items_per_page=10000, use_cache=True):
    """Fetch every page for a request, returning the concatenated `data` list."""
    page = 0
    out = []
    while True:
        p = dict(params, page=page, itemsPerPage=items_per_page)
        result = _get(endpoint, p, use_cache=use_cache)
        out.extend(result["data"])
        n_pages = result["paging_info"]["numberOfPages"]
        page += 1
        if page >= n_pages:
            break
    return out


def _chunks(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def lookup_genes(gene_ids, gencode_version=DEFAULT_GENCODE_VERSION, batch_size=300,
                  use_cache=True):
    """
    Resolve a list of gene symbols and/or (versioned or unversioned) GENCODE
    IDs to GTEx reference/gene records. A given symbol may resolve to more
    than one record (e.g. a coding/non-coding pair sharing a name).
    Returns a flat list of dicts as returned by the API.
    """
    out = []
    for chunk in _chunks(sorted(set(gene_ids)), batch_size):
        out.extend(_paginate(
            "reference/gene",
            {"geneId": chunk, "gencodeVersion": gencode_version},
            use_cache=use_cache,
        ))
    return out


# GENCODE annotation pruning between releases means a handful of (mostly
# small RNA) genes that exist at v26 are absent from v39: e.g. MIRLET7A1 has
# a v26 record but no v39 one. Each version is paired with the GTEx dataset
# built on that annotation, so a gene only findable at v26 can still only be
# queried for expression against gtex_v8, not gtex_v10.
VERSION_DATASET_PAIRS = [
    (DEFAULT_GENCODE_VERSION, DEFAULT_DATASET),  # ("v39", "gtex_v10") -- preferred: latest, largest cohort
    ("v26", "gtex_v8"),                          # fallback for genes pruned out of v39
]


def lookup_genes_with_fallback(gene_ids, batch_size=300, use_cache=True):
    """
    Like lookup_genes, but tries each (gencode_version, dataset_id) pair in
    VERSION_DATASET_PAIRS in order and only falls back to the next version
    for IDs that found nothing in the previous one.

    Returns dict: gene_id (the literal symbol/Ensembl ID queried) -> list of
    (record, dataset_id) tuples, since one query value can still resolve to
    more than one GENCODE record within a single version.
    """
    remaining = set(gene_ids)
    out = {}
    for gencode_version, dataset_id in VERSION_DATASET_PAIRS:
        if not remaining:
            break
        hits = lookup_genes(remaining, gencode_version=gencode_version,
                             batch_size=batch_size, use_cache=use_cache)
        found_for = set()
        for rec in hits:
            # the literal query value may have been a symbol or an Ensembl ID;
            # match back on both so the caller's lookup key always works
            for key in (rec.get("geneSymbol"), rec.get("geneSymbolUpper"),
                        rec.get("gencodeId", "").split(".")[0]):
                if key in remaining:
                    out.setdefault(key, []).append((rec, dataset_id))
                    found_for.add(key)
        remaining -= found_for
    return out


def median_gene_expression(gencode_ids, dataset_id=DEFAULT_DATASET, batch_size=500,
                            use_cache=True):
    """Median TPM per tissue for each (versioned) gencodeId. Returns a flat list."""
    out = []
    for chunk in _chunks(sorted(set(gencode_ids)), batch_size):
        out.extend(_paginate(
            "expression/medianGeneExpression",
            {"gencodeId": chunk, "datasetId": dataset_id},
            use_cache=use_cache,
        ))
    return out


def median_transcript_expression(gencode_ids, dataset_id=DEFAULT_DATASET, batch_size=50,
                                  use_cache=True):
    """Median TPM per tissue per transcript for each (versioned) gencodeId.

    Batch size is kept small relative to median_gene_expression because the
    response fans out by transcript count as well as tissue count (54), and
    some ncRNA host genes carry hundreds of annotated transcripts.
    """
    out = []
    for chunk in _chunks(sorted(set(gencode_ids)), batch_size):
        out.extend(_paginate(
            "expression/medianTranscriptExpression",
            {"gencodeId": chunk, "datasetId": dataset_id},
            use_cache=use_cache,
        ))
    return out
