"""Tool 2: fetch a protein's AlphaFold 3D-structure confidence (free, no API key).

Reuses the EXACT pLDDT maths and confidence bands from Phases 0-1 (one source of
truth), so the agent's structural claims match our dataset.
"""

import requests
from langchain_core.tools import tool

# Reuse the shared helpers from the earlier phases (src/ is on the import path
# when the agent runs, so these top-level modules are importable).
from fetch_proteins import mean_plddt
from build_dataset import confidence_band

HEADERS = {"User-Agent": "autoantigen-atlas (educational project)"}
META_URL = "https://alphafold.ebi.ac.uk/api/prediction/{}"
ENTRY_URL = "https://alphafold.ebi.ac.uk/entry/{}"


@tool
def get_structure(accession: str) -> str:
    """Get the AlphaFold predicted 3D-structure confidence for a protein, given its
    UniProt accession (example: "P10155"). Returns the mean pLDDT confidence score
    (0-100), a plain-English confidence band, and a citation URL. Call this after
    lookup_protein (which gives you the accession) to report how reliable the
    predicted structure is.
    """
    try:
        meta = requests.get(META_URL.format(accession), headers=HEADERS, timeout=30)
        meta.raise_for_status()
        predictions = meta.json()
        if not predictions:
            return "AlphaFold has no predicted structure for accession '{}'.".format(accession)
        pdb_url = predictions[0]["pdbUrl"]

        pdb = requests.get(pdb_url, headers=HEADERS, timeout=60)
        pdb.raise_for_status()
        plddt = mean_plddt(pdb.text)
    except Exception as error:
        return "TOOL ERROR (AlphaFold): {}. Could not fetch structure for '{}'.".format(
            error, accession)

    if plddt is None:
        return "Could not compute a confidence score for '{}'.".format(accession)

    band = confidence_band(plddt)
    citation = ENTRY_URL.format(accession)
    return (
        "AlphaFold structure confidence for {acc}:\n"
        "- Mean pLDDT: {plddt} / 100\n"
        "- Confidence band: {band}  (AlphaFold: Very high >90, Confident 70-90, "
        "Low 50-70, Very low <50)\n"
        "- Source (cite this): {citation}"
    ).format(acc=accession, plddt=plddt, band=band, citation=citation)
