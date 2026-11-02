"""Tool 3: search the biomedical literature via Europe PMC (free, no API key).

Returns real papers with citation links, so the agent can ground claims in
published evidence instead of its own memory.
"""

import requests
from langchain_core.tools import tool

HEADERS = {"User-Agent": "lupus-autoantigen-assistant (educational project)"}
SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


@tool
def search_literature(query: str) -> str:
    """Search the biomedical research literature (Europe PMC) for published papers
    about a topic (examples: "Ro60 lupus autoantigen", "anti-C1q lupus nephritis").
    Returns up to 5 papers with title, authors, year, and a citation link. Use this
    to find published evidence supporting any claim about lupus or its autoantigens.
    """
    try:
        params = {"query": query, "format": "json", "pageSize": 6, "resultType": "lite"}
        resp = requests.get(SEARCH_URL, params=params, headers=HEADERS, timeout=30)
        resp.raise_for_status()
        results = resp.json().get("resultList", {}).get("result", [])
    except Exception as error:
        return "TOOL ERROR (Europe PMC): {}. Could not search for '{}'.".format(error, query)

    if not results:
        return "No papers found in Europe PMC for '{}'. Try different search terms.".format(query)

    lines = ["Europe PMC results for '{}':".format(query)]
    for i, paper in enumerate(results, 1):
        title = paper.get("title", "(no title)")
        authors = paper.get("authorString", "(authors unavailable)")
        year = paper.get("pubYear", "n.d.")
        source = paper.get("source", "")
        pid = paper.get("id", "")
        # Build a stable citation link (cite this)
        if source and pid:
            link = "https://europepmc.org/article/{}/{}".format(source, pid)
        elif paper.get("doi"):
            link = "https://doi.org/{}".format(paper["doi"])
        else:
            link = "(no stable link)"
        lines.append(
            "{i}. {title} ({year})\n   {authors}\n   Cite: {link}".format(
                i=i, title=title, year=year, authors=authors, link=link))
    return "\n".join(lines)
