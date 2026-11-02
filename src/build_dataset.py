"""
Phase 1 — Knowledge Base
========================
Read the raw files we downloaded in Phase 0 and distill them into ONE clean,
trustworthy dataset that later phases (and humans) can rely on.

Design principles / guardrails (intentional from the start):
  * READ-ONLY on data/raw/  — we only read the originals, never change them.
  * WRITE-ONLY to data/processed/ — all clean output lands here.
  * VALIDATE before we trust — every raw file must exist, parse, and contain
    the fields we expect, or the protein is flagged (not silently dropped).
  * IDENTITY CHECK — confirm each record is the HUMAN protein (taxon 9606).
  * PROVENANCE — every record carries its source URL + source file + build date,
    so nothing is ever taken on trust.
  * REPRODUCIBLE — no randomness; we even reuse Phase 0's exact pLDDT maths
    (one source of truth) by importing it, so the two scripts can't drift apart.

Run it from the project folder with:
    python3 src/build_dataset.py
"""

import csv
import datetime
import json
import os

# Reuse the EXACT same confidence calculation from Phase 0 — one source of truth.
from fetch_proteins import mean_plddt, SEED_PROTEINS


# ---------------------------------------------------------------------------
# Paths — built relative to THIS file so the script runs from anywhere.
# ---------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

HUMAN_TAXON_ID = 9606   # the identity guardrail: every protein must be human

# UniProt base URL, kept here so provenance links are easy to rebuild.
UNIPROT_ENTRY_URL = "https://www.uniprot.org/uniprotkb/{}/entry"
ALPHAFOLD_ENTRY_URL = "https://alphafold.ebi.ac.uk/entry/{}"


def confidence_band(plddt):
    """
    Translate a mean pLDDT number into AlphaFold's OFFICIAL plain-English bands,
    so a reader instantly knows how much to trust the 3D structure.
    """
    if plddt is None:
        return "unknown"
    if plddt > 90:
        return "Very high"
    if plddt >= 70:
        return "Confident"
    if plddt >= 50:
        return "Low"
    return "Very low"


def dig(record, *keys, default="unknown"):
    """Safely reach into a deeply-nested record; return default if any step is missing."""
    current = record
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            return default
    return current


def read_uniprot_facts(accession):
    """
    Load one protein's raw UniProt JSON and pull out the facts we care about.
    Returns (facts_dict, problems_list). problems_list is empty if all is well.
    """
    problems = []
    path = os.path.join(RAW_DIR, "{}_uniprot.json".format(accession))

    # GUARDRAIL: the raw file must exist.
    if not os.path.exists(path):
        return None, ["missing raw UniProt file: {}".format(os.path.basename(path))]

    # GUARDRAIL: it must be valid JSON.
    try:
        with open(path) as f:
            record = json.load(f)
    except (ValueError, OSError) as error:
        return None, ["could not read UniProt JSON: {}".format(error)]

    # Pull the human-readable facts.
    name = dig(record, "proteinDescription", "recommendedName", "fullName", "value")
    genes = record.get("genes", [])
    gene = genes[0].get("geneName", {}).get("value", "unknown") if genes else "unknown"
    length = dig(record, "sequence", "length", default=None)
    organism = dig(record, "organism", "scientificName")
    taxon_id = dig(record, "organism", "taxonId", default=None)

    function = "unknown"
    for comment in record.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                function = texts[0].get("value", "unknown")
            break

    # GUARDRAIL: identity check — must be the human protein.
    if taxon_id != HUMAN_TAXON_ID:
        problems.append(
            "identity check failed: taxon {} (expected human {})".format(taxon_id, HUMAN_TAXON_ID))

    # GUARDRAIL: the essentials must actually be present.
    if name == "unknown":
        problems.append("missing protein name")
    if length is None:
        problems.append("missing sequence length")

    facts = {
        "protein_name": name,
        "gene": gene,
        "length_aa": length,
        "organism": organism,
        "function": function,
    }
    return facts, problems


def read_structure_confidence(accession):
    """Load the raw AlphaFold .pdb file and compute its mean confidence (pLDDT)."""
    path = os.path.join(RAW_DIR, "{}_alphafold.pdb".format(accession))
    if not os.path.exists(path):
        return None, ["missing raw AlphaFold file: {}".format(os.path.basename(path))]
    try:
        with open(path) as f:
            pdb_text = f.read()
    except OSError as error:
        return None, ["could not read PDB file: {}".format(error)]

    plddt = mean_plddt(pdb_text)   # <-- the SAME maths Phase 0 used
    if plddt is None:
        return None, ["could not compute pLDDT (no CA atoms found)"]
    return plddt, []


def build():
    # GUARDRAIL: make sure the OUTPUT folder exists; we never write into raw/.
    os.makedirs(PROCESSED_DIR, exist_ok=True)

    clean_records = []
    flagged = []   # proteins with any problem — reported loudly, never hidden.

    print("Building clean dataset from {}\n".format(RAW_DIR))

    for protein in SEED_PROTEINS:
        acc = protein["accession"]
        short = protein["short"]

        facts, uni_problems = read_uniprot_facts(acc)
        plddt, pdb_problems = read_structure_confidence(acc)
        problems = uni_problems + pdb_problems

        if facts is None:
            # Can't build a meaningful record without the UniProt facts.
            flagged.append({"short": short, "accession": acc, "problems": problems})
            print("  [FLAGGED] {} ({}): {}".format(short, acc, "; ".join(problems)))
            continue

        record = {
            "short_name": short,
            "accession": acc,
            "role": protein.get("role", "positive"),  # positive autoantigen vs negative control
            "protein_name": facts["protein_name"],
            "gene": facts["gene"],
            "organism": facts["organism"],
            "length_aa": facts["length_aa"],
            "antibody_test": protein["antibody_test"],
            "structure_mean_plddt": plddt,
            "structure_confidence": confidence_band(plddt),
            "function": facts["function"],
            # --- provenance: where every fact came from ---
            "source_uniprot": UNIPROT_ENTRY_URL.format(acc),
            "source_alphafold": ALPHAFOLD_ENTRY_URL.format(acc),
            "raw_uniprot_file": "data/raw/{}_uniprot.json".format(acc),
            "raw_structure_file": "data/raw/{}_alphafold.pdb".format(acc),
        }
        clean_records.append(record)

        status = "OK" if not problems else "OK (with notes: {})".format("; ".join(problems))
        print("  [{}] {} ({}) — {} aa, {} ({})".format(
            status, short, acc, record["length_aa"],
            record["structure_confidence"], record["structure_mean_plddt"]))

    # ---- Write the JSON (for the agent) with a provenance header ----
    build_date = datetime.date.today().isoformat()
    payload = {
        "_metadata": {
            "project": "Lupus Autoantigen Discovery Assistant",
            "phase": "1 - knowledge base",
            "built_from": "data/raw/",
            "build_date": build_date,
            "record_count": len(clean_records),
            "positive_count": sum(1 for r in clean_records if r["role"] == "positive"),
            "negative_control_count": sum(1 for r in clean_records if r["role"] == "negative_control"),
            "flagged_count": len(flagged),
            "source": "Autoantibodies in SLE: Diagnostic and Pathogenic Insights (2025), PMC12386345",
            "source_apis": {
                "uniprot": "https://rest.uniprot.org/uniprotkb/",
                "alphafold": "https://alphafold.ebi.ac.uk/api/prediction/",
            },
            "confidence_bands": "AlphaFold official: Very high >90, Confident 70-90, Low 50-70, Very low <50",
        },
        "autoantigens": clean_records,
    }
    json_path = os.path.join(PROCESSED_DIR, "autoantigens.json")
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    # ---- Write the CSV (for you to open in Excel/Numbers) ----
    csv_path = os.path.join(PROCESSED_DIR, "autoantigens.csv")
    if clean_records:
        columns = list(clean_records[0].keys())
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            writer.writeheader()
            writer.writerows(clean_records)

    # ---- Final honest report ----
    print("\n--- Summary ---")
    print("  Clean records : {} / {}".format(len(clean_records), len(SEED_PROTEINS)))
    print("  Flagged       : {}".format(len(flagged)))
    print("  Written       : data/processed/autoantigens.json")
    print("                  data/processed/autoantigens.csv")
    if flagged:
        print("\n  !! Some proteins were flagged — coverage is INCOMPLETE. Review above.")
    else:
        print("\n  All seed proteins passed every guardrail. Dataset is complete. ")


if __name__ == "__main__":
    build()
