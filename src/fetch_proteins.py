"""
Phase 0 — Foundation
====================
Fetch our 6 seed lupus autoantigens and their AlphaFold 3D structures,
using only FREE public APIs (no API keys needed!).

For each protein this script:
  1. Downloads the UniProt record (name, gene, sequence, function)  -> data/raw/<ACC>_uniprot.json
  2. Downloads the AlphaFold 3D structure                            -> data/raw/<ACC>_alphafold.pdb
  3. Calculates the average confidence (mean pLDDT) of that structure

Run it from the project folder with:
    python3 src/fetch_proteins.py

Everything it downloads lands in  data/raw/  — the "untouched originals" folder.
We never hand-edit those files; later phases will read them and build a clean dataset.
"""

import json
import os
import sys
import time

import requests

# The seed list now lives in the disease registry (src/diseases.py), so this one
# script can fetch ANY disease's proteins. Pass a disease key on the command line
# (defaults to lupus).
from diseases import seed_for

# Where to save everything. We build the path relative to THIS file so the
# script works no matter which folder you run it from.
HERE = os.path.dirname(os.path.abspath(__file__))          # .../src
PROJECT_ROOT = os.path.dirname(HERE)                        # .../autoantigen-atlas
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")

# A short, polite label so the servers know who's calling (good manners for free APIs).
HEADERS = {"User-Agent": "autoantigen-atlas (educational project)"}


def fetch_uniprot(accession):
    """Download one protein's full record from UniProt as JSON."""
    url = "https://rest.uniprot.org/uniprotkb/{}.json".format(accession)
    response = requests.get(url, headers=HEADERS, timeout=30)
    response.raise_for_status()   # stop loudly if the download failed
    return response.json()


def summarize_uniprot(record):
    """Pull a few friendly, human-readable facts out of the big UniProt record."""
    # UniProt records are deeply nested, so we dig carefully and fall back to
    # "unknown" if a field happens to be missing.
    name = (record.get("proteinDescription", {})
                  .get("recommendedName", {})
                  .get("fullName", {})
                  .get("value", "unknown"))

    genes = record.get("genes", [])
    gene = genes[0].get("geneName", {}).get("value", "unknown") if genes else "unknown"

    length = record.get("sequence", {}).get("length", "unknown")

    # The "function" is stored as a FUNCTION comment; grab the first one if present.
    function = "unknown"
    for comment in record.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                function = texts[0].get("value", "unknown")
            break

    return {"name": name, "gene": gene, "length": length, "function": function}


def fetch_alphafold_pdb(accession):
    """Ask AlphaFold where this protein's structure lives, then download the PDB file."""
    # Step 1: the metadata endpoint tells us the download URL for the structure.
    meta_url = "https://alphafold.ebi.ac.uk/api/prediction/{}".format(accession)
    meta_resp = requests.get(meta_url, headers=HEADERS, timeout=30)
    meta_resp.raise_for_status()
    predictions = meta_resp.json()
    if not predictions:
        raise RuntimeError("AlphaFold has no structure for {}".format(accession))

    pdb_url = predictions[0]["pdbUrl"]   # the ".pdb" is the classic 3D-structure file

    # Step 2: actually download the structure file (it's plain text).
    pdb_resp = requests.get(pdb_url, headers=HEADERS, timeout=60)
    pdb_resp.raise_for_status()
    return pdb_resp.text


def mean_plddt(pdb_text):
    """
    Work out the average confidence of the structure.

    AlphaFold stores its per-residue confidence (pLDDT, 0-100) in the
    "B-factor" column of the PDB file. We read that column for each
    main backbone atom (the "CA" / alpha-carbon) and average it.
    Higher = AlphaFold is more sure about that region's shape.
    """
    scores = []
    for line in pdb_text.splitlines():
        # Only ATOM lines describe atoms; we pick the CA atom (one per residue).
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            scores.append(float(line[60:66]))   # the B-factor column
    if not scores:
        return None
    return round(sum(scores) / len(scores), 1)


def main():
    # Which disease? First command-line word, defaulting to lupus (so the old
    # `python3 src/fetch_proteins.py` command still works exactly as before).
    disease = sys.argv[1] if len(sys.argv) > 1 else "lupus"
    seed = seed_for(disease)

    # Make sure the destination folder exists (harmless if it already does).
    os.makedirs(RAW_DIR, exist_ok=True)

    print("Fetching {} seed proteins for '{}' into {}\n".format(len(seed), disease, RAW_DIR))

    for protein in seed:
        short = protein["short"]
        acc = protein["accession"]
        print("--- {} ({}) ---".format(short, acc))

        try:
            # 1) UniProt record -------------------------------------------------
            record = fetch_uniprot(acc)
            uniprot_path = os.path.join(RAW_DIR, "{}_uniprot.json".format(acc))
            with open(uniprot_path, "w") as f:
                json.dump(record, f, indent=2)
            facts = summarize_uniprot(record)
            print("  UniProt : {}  |  gene {}  |  {} aa".format(
                facts["name"], facts["gene"], facts["length"]))

            # 2) AlphaFold structure -------------------------------------------
            pdb_text = fetch_alphafold_pdb(acc)
            pdb_path = os.path.join(RAW_DIR, "{}_alphafold.pdb".format(acc))
            with open(pdb_path, "w") as f:
                f.write(pdb_text)
            confidence = mean_plddt(pdb_text)
            print("  AlphaFold: structure saved  |  mean pLDDT {}".format(confidence))

        except Exception as error:
            # One protein failing shouldn't stop the whole run — report and move on.
            print("  !! Problem with {}: {}".format(acc, error))

        # Be a polite API citizen: a short pause between proteins.
        time.sleep(1)
        print()

    print("Done. Raw files are in data/raw/ — ready for Phase 1 (building the clean dataset).")


if __name__ == "__main__":
    main()
