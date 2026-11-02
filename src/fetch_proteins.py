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
import time

import requests


# ---------------------------------------------------------------------------
# Our starter seed list — six well-known lupus autoantigens.
# (UniProt "accession" is the stable ID that every database agrees on.)
# ---------------------------------------------------------------------------
#   role: "positive" = an established SLE autoantigen (from the 2025 review,
#         PMC12386345); "negative_control" = an abundant protein NOT established
#         as a lupus autoantigen, used in Phase 3 to test the agent's specificity.
#   Every accession verified live against UniProt on 2026-11-01.
SEED_PROTEINS = [
    # --- Ro / La ---
    {"short": "Ro60",        "accession": "P10155", "gene": "RO60",    "antibody_test": "anti-Ro/SSA",          "role": "positive"},
    {"short": "Ro52",        "accession": "P19474", "gene": "TRIM21",  "antibody_test": "anti-Ro52",            "role": "positive"},
    {"short": "La/SSB",      "accession": "P05455", "gene": "SSB",     "antibody_test": "anti-La/SSB",          "role": "positive"},
    # --- Sm ---
    {"short": "Sm-B/B'",     "accession": "P14678", "gene": "SNRPB",   "antibody_test": "anti-Sm",              "role": "positive"},
    {"short": "Sm-D1",       "accession": "P62314", "gene": "SNRPD1",  "antibody_test": "anti-Sm",              "role": "positive"},
    {"short": "Sm-D2",       "accession": "P62316", "gene": "SNRPD2",  "antibody_test": "anti-Sm",              "role": "positive"},
    {"short": "Sm-D3",       "accession": "P62318", "gene": "SNRPD3",  "antibody_test": "anti-Sm",              "role": "positive"},
    # --- U1-RNP ---
    {"short": "U1-70K",      "accession": "P08621", "gene": "SNRNP70", "antibody_test": "anti-U1-RNP",          "role": "positive"},
    {"short": "U1-A",        "accession": "P09012", "gene": "SNRPA",   "antibody_test": "anti-U1-RNP",          "role": "positive"},
    {"short": "U1-C",        "accession": "P09234", "gene": "SNRPC",   "antibody_test": "anti-U1-RNP",          "role": "positive"},
    # --- nucleosome / histones ---
    {"short": "Histone H4",  "accession": "P62805", "gene": "H4C1",    "antibody_test": "anti-histone/nucleosome", "role": "positive"},
    {"short": "Histone H3.1","accession": "P68431", "gene": "H3C1",    "antibody_test": "anti-histone/nucleosome", "role": "positive"},
    {"short": "Histone H2B", "accession": "P62807", "gene": "H2BC4",   "antibody_test": "anti-histone/nucleosome", "role": "positive"},
    {"short": "Histone H2A", "accession": "P0C0S8", "gene": "H2AC11",  "antibody_test": "anti-histone/nucleosome", "role": "positive"},
    # --- ribosomal P ---
    {"short": "Ribosomal P0","accession": "P05388", "gene": "RPLP0",   "antibody_test": "anti-ribosomal P",     "role": "positive"},
    {"short": "Ribosomal P1","accession": "P05386", "gene": "RPLP1",   "antibody_test": "anti-ribosomal P",     "role": "positive"},
    {"short": "Ribosomal P2","accession": "P05387", "gene": "RPLP2",   "antibody_test": "anti-ribosomal P",     "role": "positive"},
    # --- complement C1q ---
    {"short": "C1q A",       "accession": "P02745", "gene": "C1QA",    "antibody_test": "anti-C1q",             "role": "positive"},
    {"short": "C1q B",       "accession": "P02746", "gene": "C1QB",    "antibody_test": "anti-C1q",             "role": "positive"},
    {"short": "C1q C",       "accession": "P02747", "gene": "C1QC",    "antibody_test": "anti-C1q",             "role": "positive"},
    # --- DNA-repair (overlap) ---
    {"short": "Ku70",        "accession": "P12956", "gene": "XRCC6",   "antibody_test": "anti-Ku",              "role": "positive"},
    {"short": "Ku80",        "accession": "P13010", "gene": "XRCC5",   "antibody_test": "anti-Ku",              "role": "positive"},
    # --- proliferation ---
    {"short": "PCNA",        "accession": "P12004", "gene": "PCNA",    "antibody_test": "anti-PCNA",            "role": "positive"},
    # --- antiphospholipid ---
    {"short": "Beta-2-GPI",  "accession": "P02749", "gene": "APOH",    "antibody_test": "anti-beta2GPI",        "role": "positive"},
    # --- negative controls (NOT established lupus autoantigens) ---
    {"short": "GAPDH",           "accession": "P04406", "gene": "GAPDH", "antibody_test": "(none)", "role": "negative_control"},
    {"short": "Beta-actin",      "accession": "P60709", "gene": "ACTB",  "antibody_test": "(none)", "role": "negative_control"},
    {"short": "Serum albumin",   "accession": "P02768", "gene": "ALB",   "antibody_test": "(none)", "role": "negative_control"},
    {"short": "Hemoglobin beta", "accession": "P68871", "gene": "HBB",   "antibody_test": "(none)", "role": "negative_control"},
    {"short": "Tubulin beta",    "accession": "P07437", "gene": "TUBB",  "antibody_test": "(none)", "role": "negative_control"},
]

# Where to save everything. We build the path relative to THIS file so the
# script works no matter which folder you run it from.
HERE = os.path.dirname(os.path.abspath(__file__))          # .../src
PROJECT_ROOT = os.path.dirname(HERE)                        # .../lupus-autoantigen-assistant
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")

# A short, polite label so the servers know who's calling (good manners for free APIs).
HEADERS = {"User-Agent": "lupus-autoantigen-assistant (educational project)"}


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
    # Make sure the destination folder exists (harmless if it already does).
    os.makedirs(RAW_DIR, exist_ok=True)

    print("Fetching {} seed proteins into {}\n".format(len(SEED_PROTEINS), RAW_DIR))

    for protein in SEED_PROTEINS:
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
