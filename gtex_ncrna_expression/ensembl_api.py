"""
Thin client for the Ensembl REST API (https://rest.ensembl.org), used here
for exactly one thing: fetching the spliced cDNA sequence of GENCODE/Ensembl
transcripts (ENST IDs) so they can be compared, by sequence content rather
than by ID, against the sequences RNAInterSEQ already has on file for each
ncRNA entity.

This cross-database comparison is necessary because RNAInterSEQ's own
sequence for an entity may have been sourced from miRBase, NCBI RefSeq, or
RNAcentral rather than Ensembl -- an ID match against GTEx's GENCODE
transcript IDs is not possible in those cases, but the underlying sequence
is the same biological molecule regardless of which database first
deposited it, so sequence identity is the only join key that works across
all of them.

Single dict-based disk cache (transcript_id -> sequence or None) rather than
the gtex_api.py per-request-hash cache, since callers ask for overlapping
transcript sets across runs and an incremental dict lets every previously
seen ID short-circuit immediately.
"""
import json
import time
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter, Retry

BASE_URL = "https://rest.ensembl.org"
CACHE_FILE = Path(__file__).parent / "cache" / "ensembl_transcript_sequences.json"
BATCH_SIZE = 50  # hard server-side limit on POST /sequence/id


def _session():
    s = requests.Session()
    retries = Retry(total=5, backoff_factor=0.5,
                     status_forcelist=[429, 500, 502, 503, 504])
    s.mount("https://", HTTPAdapter(max_retries=retries))
    return s


_SESSION = _session()


def fetch_transcript_sequences(transcript_ids, seq_type="cdna", sleep=0.1):
    """
    Return {unversioned_transcript_id: sequence_or_None} for every ID in
    transcript_ids (cached across calls/runs). `seq_type="cdna"` returns the
    spliced exonic sequence, matching how RNAInter's own ncRNA sequences are
    represented.
    """
    transcript_ids = sorted({t.split(".")[0] for t in transcript_ids})

    cache = {}
    if CACHE_FILE.exists():
        cache = json.loads(CACHE_FILE.read_text())

    todo = [t for t in transcript_ids if t not in cache]
    n_failed_batches = 0
    for i in range(0, len(todo), BATCH_SIZE):
        chunk = todo[i:i + BATCH_SIZE]
        for attempt in range(3):
            try:
                resp = _SESSION.post(
                    f"{BASE_URL}/sequence/id",
                    headers={"Content-Type": "application/json", "Accept": "application/json"},
                    json={"ids": chunk, "type": seq_type},
                    timeout=120,
                )
                resp.raise_for_status()
                results = resp.json()
                found = {r["query"]: r["seq"] for r in results if "error" not in r}
                for tid in chunk:
                    cache[tid] = found.get(tid)
                break
            except (requests.exceptions.RequestException, ValueError):
                if attempt == 2:
                    n_failed_batches += 1  # leave this chunk uncached; next run retries it
                else:
                    time.sleep(2 * (attempt + 1))
        if (i // BATCH_SIZE) % 20 == 0:
            CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
            CACHE_FILE.write_text(json.dumps(cache))
        if sleep:
            time.sleep(sleep)

    CACHE_FILE.parent.mkdir(parents=True, exist_ok=True)
    CACHE_FILE.write_text(json.dumps(cache))
    if n_failed_batches:
        print(f"  warning: {n_failed_batches} batch(es) failed after retries and were skipped "
              f"(their IDs remain unresolved this run; rerun to retry just those).")
    return {t: cache.get(t) for t in transcript_ids}
