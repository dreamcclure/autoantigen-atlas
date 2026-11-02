"""Tool 1: look up a human protein's facts in UniProt (free, no API key).

Returns the facts AND a citation URL, so the agent can cite every claim.

Lookup strategy (improved in the Phase 3 recall pass):
  1. If the query is a UniProt accession (e.g. "P09012"), fetch it directly.
  2. Otherwise search by GENE symbol first (exact, reliable), then fall back to
     a general name search.
We also surface UniProt's curated SUBUNIT / complex annotation, because many lupus
autoantigens are recognized at the *complex* level (the Sm core, U1 snRNP, etc.).
"""

import re

import requests
from langchain_core.tools import tool

HEADERS = {"User-Agent": "lupus-autoantigen-assistant (educational project)"}
SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{}.json"
ACCESSION_RE = re.compile(r"[A-Z][0-9][A-Z0-9]{3}[0-9]([A-Z][A-Z0-9]{2}[0-9]){0,2}$")


def _dig(record, *keys, default="unknown"):
    """Safely reach into a nested record; return default if any step is missing."""
    current = record
    for key in keys:
        if isinstance(current, dict) and key in current:
            current = current[key]
        else:
            return default
    return current


def _find_entry(query):
    """Resolve a query to one UniProt record, most-reliable method first."""
    q = query.strip()
    # 1. Direct accession fetch (most reliable).
    if ACCESSION_RE.fullmatch(q.upper()):
        try:
            r = requests.get(ENTRY_URL.format(q.upper()), headers=HEADERS, timeout=30)
            if r.ok:
                return r.json()
        except Exception:
            pass
    # 2. Gene symbol first, then a general search.
    for expr in ("gene:{}".format(q), "({})".format(q)):
        try:
            params = {"query": "{} AND organism_id:9606 AND reviewed:true".format(expr),
                      "format": "json", "size": 1}
            r = requests.get(SEARCH_URL, params=params, headers=HEADERS, timeout=30)
            r.raise_for_status()
            results = r.json().get("results", [])
            if results:
                return results[0]
        except Exception:
            continue
    return None


@tool
def lookup_protein(query: str) -> str:
    """Look up a human protein in the UniProt database by name, gene symbol, or
    accession (examples: "Ro60", "TRIM21", "SNRPA", "P10155"). Returns the official
    protein name, gene, UniProt accession, length, biological function, its
    subunit/complex membership, and a citation URL. Prefer searching by GENE symbol
    for reliability. Use this to get factual information about any protein before
    making claims about it.
    """
    try:
        rec = _find_entry(query)
    except Exception as error:
        return "TOOL ERROR (UniProt): {}. Could not look up '{}'.".format(error, query)

    if rec is None:
        return ("No reviewed human protein found in UniProt for '{}'. "
                "Try the gene symbol or the UniProt accession.").format(query)

    acc = rec.get("primaryAccession", "unknown")
    name = _dig(rec, "proteinDescription", "recommendedName", "fullName", "value")
    genes = rec.get("genes", [])
    gene = genes[0].get("geneName", {}).get("value", "unknown") if genes else "unknown"
    length = _dig(rec, "sequence", "length", default="unknown")

    function = "No function annotation available."
    subunit = "No subunit/complex annotation available."
    for comment in rec.get("comments", []):
        ctype = comment.get("commentType")
        texts = comment.get("texts", [])
        if ctype == "FUNCTION" and texts:
            function = texts[0].get("value", function)
        elif ctype == "SUBUNIT" and texts:
            subunit = texts[0].get("value", subunit)

    citation = "https://www.uniprot.org/uniprotkb/{}/entry".format(acc)
    return (
        "UniProt record for '{query}':\n"
        "- Protein name: {name}\n"
        "- Gene: {gene}\n"
        "- UniProt accession: {acc}\n"
        "- Length: {length} amino acids\n"
        "- Function: {function}\n"
        "- Subunit / complex: {subunit}\n"
        "- Source (cite this): {citation}"
    ).format(query=query, name=name, gene=gene, acc=acc, length=length,
             function=function, subunit=subunit, citation=citation)
