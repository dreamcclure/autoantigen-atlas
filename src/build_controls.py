"""
Phase 6 (support) — a principled PROTEOME-BACKGROUND control set.
================================================================
Our original 5 negative controls (GAPDH, actin, ...) are fine for a quick specificity
check, but they're hand-picked, abundant, well-folded proteins — so comparing
autoantigens against them risks measuring "housekeeping proteins are ordered" rather
than "autoantigens are disordered" (selection bias). And 5 is far too few for real
statistical power.

This script builds the honest comparison group: a RANDOM sample of reviewed human
proteins — a "typical human protein" background. We:
  1. Pull the pool of all reviewed (Swiss-Prot) human accessions in one stream request.
  2. Exclude every protein already in our disease registry (so the groups stay disjoint),
     and entries whose text mentions "autoantigen" (a light guard, done in the query).
  3. Randomly sample from the pool with a FIXED SEED (reproducible), fetching each one's
     UniProt + AlphaFold files the same way Phase 0 did, until we have TARGET_N proteins
     that actually have an AlphaFold model.

Residual contamination (a random protein that is secretly an unrecognized autoantigen) is
rare and would only make our later test more CONSERVATIVE, never inflate it — we note this.

Run:
    python3 src/build_controls.py            # ~40 proteins, seed 42
    python3 src/build_controls.py 60 7       # target 60 proteins, seed 7
"""

import json
import os
import random
import sys
import time

import requests

# Reuse Phase 0's exact fetch + maths — one source of truth (see learn/fetch_proteins.py).
from fetch_proteins import (
    HEADERS, RAW_DIR, fetch_uniprot, summarize_uniprot, fetch_alphafold_pdb, mean_plddt,
)
from diseases import DISEASES

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

# Query for the sampling pool: reviewed human proteins, minus obvious autoantigen mentions.
POOL_QUERY = "reviewed:true AND organism_id:9606 NOT autoantigen"
STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
POOL_CACHE = os.path.join(RAW_DIR, "_human_proteome_accessions.txt")  # cache the pool list

TARGET_N = 40   # how many background proteins (with models) we want
DEFAULT_SEED = 42


def known_accessions():
    """Every accession already in our registry — excluded so the groups never overlap."""
    accs = set()
    for meta in DISEASES.values():
        for p in meta["seed"]:
            accs.add(p["accession"])
    return accs


def load_pool():
    """Fetch (and cache) the full list of candidate human accessions."""
    if os.path.exists(POOL_CACHE):
        with open(POOL_CACHE) as f:
            return [line.strip() for line in f if line.strip()]
    print("Downloading reviewed-human accession pool from UniProt (one-time)...")
    resp = requests.get(STREAM_URL, params={"query": POOL_QUERY, "format": "list"},
                        headers=HEADERS, timeout=120)
    resp.raise_for_status()
    accessions = [line.strip() for line in resp.text.splitlines() if line.strip()]
    os.makedirs(RAW_DIR, exist_ok=True)
    with open(POOL_CACHE, "w") as f:
        f.write("\n".join(accessions))
    return accessions


def fetch_one(acc):
    """Fetch UniProt facts + AlphaFold mean pLDDT for one accession.

    Idempotent: if the raw files already exist we read them instead of re-downloading.
    Returns a record dict, or None if the protein has no usable AlphaFold model.
    """
    uni_path = os.path.join(RAW_DIR, "{}_uniprot.json".format(acc))
    pdb_path = os.path.join(RAW_DIR, "{}_alphafold.pdb".format(acc))

    # --- UniProt (reuse cached raw file if present) ---
    if os.path.exists(uni_path):
        with open(uni_path) as f:
            record = json.load(f)
    else:
        record = fetch_uniprot(acc)
        with open(uni_path, "w") as f:
            json.dump(record, f, indent=2)
    facts = summarize_uniprot(record)

    # --- AlphaFold (reuse cached raw file if present) ---
    if os.path.exists(pdb_path):
        with open(pdb_path) as f:
            pdb_text = f.read()
    else:
        pdb_text = fetch_alphafold_pdb(acc)   # raises if no model exists
        with open(pdb_path, "w") as f:
            f.write(pdb_text)

    plddt = mean_plddt(pdb_text)
    if plddt is None:
        return None
    return {
        "short_name": facts["gene"] if facts["gene"] != "unknown" else acc,
        "accession": acc,
        "gene": facts["gene"],
        "protein_name": facts["name"],
        "length_aa": facts["length"],
        "role": "background",
        "structure_mean_plddt": plddt,
        "source_uniprot": "https://www.uniprot.org/uniprotkb/{}/entry".format(acc),
    }


def main():
    target = int(sys.argv[1]) if len(sys.argv) > 1 else TARGET_N
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_SEED

    pool = load_pool()
    exclude = known_accessions()
    candidates = [a for a in pool if a not in exclude]

    # Fixed-seed shuffle = a reproducible "random" order. Anyone re-running gets the
    # same sample, so the result is verifiable rather than a different draw each time.
    random.seed(seed)
    random.shuffle(candidates)

    print("Pool: {} reviewed human proteins ({} excluded as known) — target {} with models, seed {}\n"
          .format(len(pool), len(exclude), target, seed))

    records = []
    tried = 0
    for acc in candidates:
        if len(records) >= target:
            break
        tried += 1
        try:
            rec = fetch_one(acc)
        except Exception as error:
            rec = None
            print("  [skip] {} — {}".format(acc, str(error)[:60]))
        if rec is None:
            continue
        records.append(rec)
        print("  [{:>2}/{}] {:<12} {:<10} {} aa  pLDDT {}".format(
            len(records), target, rec["accession"], rec["gene"],
            rec["length_aa"], rec["structure_mean_plddt"]))
        time.sleep(0.5)   # polite pause only on real (non-cached) fetches is ideal; keep simple

    payload = {
        "_metadata": {
            "project": "Autoantigen Discovery Platform",
            "set": "proteome_background",
            "phase": "6 - structural analysis (control set)",
            "description": "Random sample of reviewed human proteins as a 'typical protein' "
                           "background for structural comparison.",
            "pool_query": POOL_QUERY,
            "pool_size": len(pool),
            "seed": seed,
            "requested": target,
            "obtained": len(records),
            "candidates_tried": tried,
            "excluded_known": len(exclude),
            "caveat": "A random draw may rarely include an unrecognized autoantigen; this is "
                      "infrequent and would only make downstream tests more conservative.",
        },
        "background": records,
    }
    out = os.path.join(PROCESSED_DIR, "proteome_background.json")
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)

    print("\nDone. {} background proteins -> data/processed/proteome_background.json".format(len(records)))


if __name__ == "__main__":
    main()
