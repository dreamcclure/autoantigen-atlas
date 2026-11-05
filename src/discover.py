"""
Phase 7 — Discovery mode (partner-arm core + hold-out validation)
=================================================================
Flip the question. Instead of "is protein X an established autoantigen?", ask "which
proteins — NOT already in our established set — might be autoantigens?" and return a
ranked, evidence-backed shortlist of CANDIDATE HYPOTHESES for expert review (cited;
never "confirmed discoveries").

THIS FILE is the foundation: the PARTNER ARM. Its biological rationale is EPITOPE
SPREADING — a real phenomenon in lupus where the immune response, having targeted one
member of a molecular complex, spreads to its physical partners. So a protein that
physically associates with many known autoantigens is a natural candidate. We read those
associations from UniProt's curated interaction data (local files, free, no agent calls).

HONESTY FIRST — we prove the method works before trusting it:
  * HOLD-OUT REDISCOVERY (leave-one-out): for each known autoantigen, hide it, then ask
    whether the partner arm — using only the OTHER known antigens — proposes it back.
    The fraction recovered is a real, measurable "discovery recall" (no answer-key cheat).
  * Everything a candidate is proposed on is spelled out (which known antigens it partners
    with), so a human can judge the evidence.

Run:
    python3 src/discover.py                       # lupus
    python3 src/discover.py rheumatoid_arthritis
"""

import json
import os
import re
import sys
import time
import datetime

import requests

from diseases import DISEASES

# Discovery-framed question for the (optional, paid) agent-assessment step. Unlike the blind
# evaluation (established YES/NO), here we expect most candidates to be NOVEL — so we ask the
# agent what evidence EXISTS, tiering it: already established / emerging hints / none found.
ASSESS_QUESTION = (
    "The protein {gene} ({acc}) {rationale}. Using your tools, investigate whether there is any "
    "PUBLISHED evidence that {gene} is itself targeted by autoantibodies or implicated as an "
    "autoantigen in {disease_name} (or a related autoimmune disease). Give a short, cited "
    "assessment of what you find.\n\n"
    "End with a single final line in exactly one of these formats:\n"
    "EVIDENCE: ESTABLISHED  -- tools found it is already a recognized autoantigen in the disease\n"
    "EVIDENCE: EMERGING     -- tools found some published hints/reports, but not established\n"
    "EVIDENCE: NONE-FOUND   -- tools found no published autoantigen evidence (a pure structural candidate)"
)
EVIDENCE_RE = re.compile(r"EVIDENCE:\s*(ESTABLISHED|EMERGING|NONE-FOUND)", re.IGNORECASE)

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

HEADERS = {"User-Agent": "autoantigen-atlas (educational project)"}
ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{}.json"
MIN_SUPPORT = 2   # a candidate must link to >=2 known antigens to be ranked (one link = fluke)

# Phase 6 biological prior: lupus autoantigens are enriched for nuclear / nucleic-acid biology.
# We boost candidates with that property — but ONLY for diseases where we actually showed the
# signal (lupus). RA antigens are NOT nuclear (Phase 6), so applying it there would be wrong.
NA_KEYWORDS = {"Nucleus", "DNA-binding", "RNA-binding"}
NA_PRIOR_DISEASES = {"lupus"}

# Two documented filters that control the confounds in candidate ranking:
MIN_DEGREE = 4        # thin-data floor: need >=4 recorded interactions for a trustworthy ratio
MIN_SPECIFICITY = 0.05  # hub floor: drop promiscuous hubs (e.g. HTT, 6/860 = 0.007)

# --- Literature arm (Europe PMC text-mining) ---
EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EPMC_ANNOT = "https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds"
UNIPROT_URI_RE = re.compile(r"uniprotkb/([A-Z0-9]+)/entry")
# Keep only mentions sitting near autoantibody/reactivity language — this is what separates
# actual autoantigen targets from drug targets / risk genes / cytokines that also co-occur
# with the disease (tested: it drops CD19/CD20/TLR7/IRF5, keeps MDA5/DNASE1L3/PRTN3/RNPC3).
AUTOAB_CUE = re.compile(r"autoantib|autoantig|autoimmun|immunoreact|seroreact|seroposit|reactiv", re.I)
LIT_MIN_ARTICLES = 2  # a candidate must recur in >=2 autoantibody-context papers (washes out mislinks)


def positives(disease):
    """The disease's established autoantigens (our 'known' set)."""
    return [p for p in DISEASES[disease]["seed"] if p.get("role", "positive") == "positive"]


def all_seed_accessions(disease):
    """Every accession already in the registry (positives + controls) — the 'known' universe."""
    return {p["accession"] for p in DISEASES[disease]["seed"]}


def partners_of(accession):
    """
    Curated binary interaction partners of one protein, from its local UniProt JSON.
    Returns a dict {partner_accession: partner_gene}. Empty if the file/annotation is absent.
    """
    path = os.path.join(RAW_DIR, "{}_uniprot.json".format(accession))
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        record = json.load(f)
    found = {}
    for comment in record.get("comments", []):
        if comment.get("commentType") != "INTERACTION":
            continue
        for inter in comment.get("interactions", []):
            two = inter.get("interactantTwo", {})
            pacc = two.get("uniProtKBAccession")
            # keep real, cross-protein, same-protein-excluded partners
            if pacc and pacc != accession and not inter.get("organismDiffer", False):
                found[pacc] = two.get("geneName", "?")
    return found


def generate_candidates(disease, known_accessions=None):
    """
    PARTNER ARM: pool the interaction partners of the 'known' autoantigens.

    `known_accessions` lets the caller restrict which antigens seed the search (used by the
    hold-out test to withhold some). Defaults to all of the disease's positives.

    Returns {candidate_acc: {"gene":..., "partner_support": N, "connects_to": [genes]}}
    where partner_support = how many distinct known antigens this candidate associates with
    (the core ranking signal — more links to known antigens = stronger candidate).
    """
    pos = positives(disease)
    if known_accessions is None:
        known_accessions = {p["accession"] for p in pos}
    gene_of = {p["accession"]: p["short"] for p in pos}

    candidates = {}
    for anchor in pos:
        if anchor["accession"] not in known_accessions:
            continue  # this known antigen is being withheld (hold-out) — can't use it to search
        for pacc, pgene in partners_of(anchor["accession"]).items():
            entry = candidates.setdefault(pacc, {"gene": pgene, "partner_support": 0, "connects_to": []})
            entry["partner_support"] += 1
            entry["connects_to"].append(gene_of.get(anchor["accession"], anchor["short"]))
    return candidates


def _ensure_cached(acc):
    """Download a candidate's UniProt record into data/raw/ if we don't already have it."""
    path = os.path.join(RAW_DIR, "{}_uniprot.json".format(acc))
    if os.path.exists(path):
        return True
    try:
        r = requests.get(ENTRY_URL.format(acc), headers=HEADERS, timeout=30)
        r.raise_for_status()
        with open(path, "w") as f:
            json.dump(r.json(), f)
        time.sleep(0.3)   # polite pause for the free API
        return True
    except Exception:
        return False


def interaction_degree(acc):
    """Total distinct interaction partners of a protein (its 'connectivity' / hub-ness)."""
    if not _ensure_cached(acc):
        return None
    return len(partners_of(acc))


def is_na_associated(acc):
    """True if the protein carries a nuclear / DNA- / RNA-binding UniProt keyword (Phase 6 prior)."""
    path = os.path.join(RAW_DIR, "{}_uniprot.json".format(acc))
    if not os.path.exists(path):
        return False
    with open(path) as f:
        record = json.load(f)
    return bool({k.get("name") for k in record.get("keywords", [])} & NA_KEYWORDS)


def holdout_recall(disease):
    """
    LEAVE-ONE-OUT rediscovery: hide each known antigen in turn and ask whether the partner
    arm (using only the others) proposes it back. Honest, answer-key-free proof of signal.
    """
    pos = positives(disease)
    all_acc = {p["accession"] for p in pos}
    recovered, missed = [], []
    for held in pos:
        others = all_acc - {held["accession"]}
        cands = generate_candidates(disease, known_accessions=others)
        if held["accession"] in cands:
            recovered.append((held["short"], cands[held["accession"]]["partner_support"]))
        else:
            missed.append(held["short"])
    recall = round(len(recovered) / len(pos), 3) if pos else None
    return recall, recovered, missed


# ---------------------------------------------------------------------------
# Shared agent-assessment step (used by both arms). `--assess N` turns it on.
# ---------------------------------------------------------------------------
def assess_count(default=5):
    """Return N if `--assess [N]` was passed (default 5), else None (assessment off)."""
    if "--assess" not in sys.argv:
        return None
    i = sys.argv.index("--assess")
    if i + 1 < len(sys.argv) and sys.argv[i + 1].isdigit():
        return int(sys.argv[i + 1])
    return default


def assess_top(disease, disease_name, items):
    """Run the trusted cited agent on each item; items = [{accession, gene, rationale}]."""
    print("\n--- Agent assessment of top {} candidates (live, cited) ---".format(len(items)))
    from agent import build_agent, ask
    agent = build_agent(disease)
    out = []
    for it in items:
        q = ASSESS_QUESTION.format(gene=it["gene"], acc=it["accession"],
                                   disease_name=disease_name, rationale=it["rationale"])
        try:
            answer = ask(agent, q)
        except Exception as error:
            answer = "AGENT ERROR: {}".format(error)
        matches = EVIDENCE_RE.findall(answer or "")
        tier = matches[-1].upper() if matches else "UNPARSED"
        print("  {:<11} {:<10} -> EVIDENCE: {}".format(it["gene"], it["accession"], tier))
        out.append({"accession": it["accession"], "gene": it["gene"],
                    "evidence_tier": tier, "transcript": answer})
    return out


# ---------------------------------------------------------------------------
# LITERATURE ARM — proteins recurring in the disease's autoantibody literature.
# Complements the partner arm: it finds candidates that already have a literature
# footing (and from adjacent autoimmune diseases), which structure-based partners miss.
# ---------------------------------------------------------------------------
def _epmc_articles(disease_name):
    """Collect Europe PMC article ids about the disease's autoantibodies/autoantigens."""
    queries = ["{} autoantibody target protein".format(disease_name),
               "{} autoantigen".format(disease_name),
               "{} anti antibody target novel".format(disease_name)]
    ids = []
    for q in queries:
        try:
            r = requests.get(EPMC_SEARCH, params={"query": q, "format": "json",
                             "pageSize": 40, "resultType": "lite"}, headers=HEADERS, timeout=30)
            for p in r.json().get("resultList", {}).get("result", []):
                if p.get("source") and p.get("id"):
                    ids.append("{}:{}".format(p["source"], p["id"]))
        except Exception:
            continue
    return list(dict.fromkeys(ids))   # dedup, preserve order


def generate_literature_candidates(disease):
    """
    Pull Gene_Proteins text-mining annotations from the disease's autoantibody papers,
    keep only mentions near autoantibody-context cues, map to reviewed HUMAN UniProt,
    and score by the number of distinct papers. Returns {acc: {gene, lit_support}}.
    """
    disease_name = DISEASES[disease]["name"]
    cache_path = os.path.join(PROCESSED_DIR, "_literature_cache_{}.json".format(disease))

    def from_cache(reason):
        """Fall back to a prior good run when Europe PMC is unavailable — fail honestly, not silently."""
        if os.path.exists(cache_path):
            print("  [note] {} — using cached literature candidates (delete {} to force a refresh)."
                  .format(reason, os.path.basename(cache_path)))
            return json.load(open(cache_path))["candidates"]
        print("  [warn] {} — and no cache to fall back on. Try again when the service is up.".format(reason))
        return {}

    articles = _epmc_articles(disease_name)
    if not articles:
        return from_cache("Europe PMC returned no articles (service may be unavailable)")

    acc_articles = {}
    for i in range(0, len(articles), 8):   # annotations API takes a small batch per call
        try:
            resp = requests.get(EPMC_ANNOT, params={"articleIds": articles[i:i + 8],
                                "type": "Gene_Proteins", "format": "JSON"}, headers=HEADERS, timeout=40)
            for art in resp.json():
                aid = art.get("extId") or art.get("id") or str(i)
                for ann in art.get("annotations", []):
                    ctx = "{} {} {}".format(ann.get("prefix", ""), ann.get("exact", ""), ann.get("postfix", ""))
                    if not AUTOAB_CUE.search(ctx):    # the context filter — the key quality step
                        continue
                    for tag in ann.get("tags", []):
                        m = UNIPROT_URI_RE.search(tag.get("uri", ""))
                        if m:
                            acc_articles.setdefault(m.group(1), set()).add(aid)
        except Exception:
            continue
        time.sleep(0.3)

    candidates = {}
    for acc, arts in acc_articles.items():
        if len(arts) < LIT_MIN_ARTICLES or not _ensure_cached(acc):
            continue
        rec = json.load(open(os.path.join(RAW_DIR, "{}_uniprot.json".format(acc))))
        if rec.get("organism", {}).get("taxonId") != 9606:   # drop the mouse/other-species mislinks
            continue
        genes = rec.get("genes", [])
        gene = genes[0].get("geneName", {}).get("value") if genes else "?"
        candidates[acc] = {"gene": gene, "lit_support": len(arts)}

    # Cache a good run (reproducible + resilient to later outages); if the live pass came back
    # empty but we have a cache, prefer the cache rather than reporting a misleading zero.
    if candidates:
        with open(cache_path, "w") as f:
            json.dump({"cached_date": datetime.date.today().isoformat(), "candidates": candidates}, f, indent=2)
    elif os.path.exists(cache_path):
        return from_cache("live extraction returned nothing")
    return candidates


def run_literature(disease, disease_name):
    """The literature-arm flow: generate, validate, rank, optionally assess, write."""
    print("Discovery mode (literature arm) — {}\n".format(disease_name))
    seed_acc = all_seed_accessions(disease)
    pos_acc = {p["accession"] for p in positives(disease)}
    cands = generate_literature_candidates(disease)
    known = {a: c for a, c in cands.items() if a in pos_acc}
    novel = {a: c for a, c in cands.items() if a not in seed_acc}

    # Soft validation: the literature arm should re-surface some KNOWN antigens.
    print("--- Validation: known antigens re-surfaced in the literature ---")
    if known:
        for a, c in sorted(known.items(), key=lambda kv: -kv[1]["lit_support"]):
            print("  {:>2} papers  {:<10} {}".format(c["lit_support"], a, c["gene"]))
    else:
        print("  (none of the seed antigens cleared the {}-paper threshold)".format(LIT_MIN_ARTICLES))

    ranked = sorted(novel.items(), key=lambda kv: -kv[1]["lit_support"])
    print("\n--- Novel literature candidates (ranked by # autoantibody-context papers) ---")
    print("  {:<6} {:<10} {}".format("papers", "accession", "gene"))
    for a, c in ranked[:20]:
        print("  {:<6} {:<10} {}".format(c["lit_support"], a, c["gene"]))

    assessments = []
    k = assess_count()
    if k is not None:
        items = [{"accession": a, "gene": c["gene"],
                  "rationale": "recurs in the {} autoantibody literature ({} papers)".format(
                      disease_name, c["lit_support"])}
                 for a, c in ranked[:k]]
        assessments = assess_top(disease, disease_name, items)

    payload = {
        "_metadata": {
            "project": "Autoantigen Discovery Platform",
            "disease_key": disease, "disease": disease_name,
            "phase": "7 - discovery mode (literature arm)",
            "run_date": datetime.date.today().isoformat(),
            "method": "literature arm: Europe PMC text-mined Gene_Proteins annotations from the "
                      "disease's autoantibody papers, kept only when near autoantibody-context cues, "
                      "mapped to reviewed human UniProt, scored by # distinct papers.",
            "context_cue": AUTOAB_CUE.pattern,
            "min_articles": LIT_MIN_ARTICLES,
            "framing": "Candidates are CITED HYPOTHESES for expert review, not confirmed discoveries.",
            "n_novel_candidates": len(novel),
            "n_known_resurfaced": len(known),
            "caveat": "Broad net: the agent-assessment step is the precision filter. Coverage of the "
                      "classic antigens is partial (recent-paper-biased); several candidates are real "
                      "autoantigens from ADJACENT autoimmune diseases, not necessarily this one.",
        },
        "known_resurfaced": [{"accession": a, "gene": c["gene"], "lit_support": c["lit_support"]}
                             for a, c in sorted(known.items(), key=lambda kv: -kv[1]["lit_support"])],
        "candidates": [{"accession": a, "gene": c["gene"], "lit_support": c["lit_support"]}
                       for a, c in ranked],
        "agent_assessments": assessments,
    }
    out = os.path.join(PROCESSED_DIR, "discovery_{}_literature.json".format(disease))
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)
    print("\n  {} novel candidates — full list in data/processed/discovery_{}_literature.json"
          .format(len(novel), disease))


def main():
    source = (sys.argv[sys.argv.index("--source") + 1]
              if "--source" in sys.argv and sys.argv.index("--source") + 1 < len(sys.argv)
              else "partner")
    # disease = first positional that isn't a flag or a flag's value
    flag_vals = set()
    for flag in ("--source", "--assess"):
        if flag in sys.argv:
            j = sys.argv.index(flag) + 1
            if j < len(sys.argv):
                flag_vals.add(sys.argv[j])
    positionals = [a for a in sys.argv[1:] if not a.startswith("--") and a not in flag_vals]
    disease = positionals[0] if positionals else "lupus"
    disease_name = DISEASES[disease]["name"]

    if source == "literature":
        run_literature(disease, disease_name)
        return

    print("Discovery mode (partner arm) — {}\n".format(disease_name))

    seed_acc = all_seed_accessions(disease)
    cands = generate_candidates(disease)

    # Split: candidates that ARE already known antigens (the validation signal) vs genuinely
    # NOVEL candidates (not in our seed at all) — the shortlist for expert review.
    rediscovered_known = {a: c for a, c in cands.items()
                          if a in {p["accession"] for p in positives(disease)}}
    novel = {a: c for a, c in cands.items() if a not in seed_acc}

    # ---- hold-out validation (does the method re-find known antigens?) ----
    recall, recovered, missed = holdout_recall(disease)
    print("--- Hold-out rediscovery (leave-one-out) ---")
    print("  Discovery recall: {}  ({} / {} known antigens re-proposed from their partners)"
          .format(recall, len(recovered), len(positives(disease))))
    if missed:
        print("  Not recoverable by partners (standalone / no interaction data): {}"
              .format(", ".join(missed)))

    # ---- CONTROL FOR HUB BIAS: rank by SPECIFICITY, not raw count ----
    # A candidate linking to many antigens is only meaningful if it isn't linked to
    # *everything*. specificity = (# known antigens it touches) / (its total partners).
    # So we fetch each shortlist candidate's total connectivity and normalize.
    use_prior = disease in NA_PRIOR_DISEASES
    scored = {a: c for a, c in novel.items() if c["partner_support"] >= MIN_SUPPORT}
    print("\n  Scoring {} candidates (support>={}) — fetching connectivity{}..."
          .format(len(scored), MIN_SUPPORT, "; flagging nuclear/NA prior" if use_prior else ""))
    for a, c in scored.items():
        deg = interaction_degree(a)
        c["total_partners"] = deg
        c["specificity"] = round(c["partner_support"] / deg, 3) if deg else None
        c["na_associated"] = is_na_associated(a) if use_prior else None

    # apply the two confound filters, then rank by support (epitope-spreading strength),
    # specificity (hub control), and the nuclear/NA prior as a final tiebreak.
    shortlist = {a: c for a, c in scored.items()
                 if c["total_partners"] and c["total_partners"] >= MIN_DEGREE
                 and c["specificity"] is not None and c["specificity"] >= MIN_SPECIFICITY}
    ranked = sorted(shortlist.items(),
                    key=lambda kv: (kv[1]["partner_support"], kv[1]["specificity"],
                                    1 if kv[1]["na_associated"] else 0), reverse=True)

    print("  {} pass the filters (degree>={}, specificity>={})\n"
          .format(len(shortlist), MIN_DEGREE, MIN_SPECIFICITY))
    print("--- Top candidate autoantigens (ranked by support, then specificity) ---")
    print("  {:<10} {:<11} {:<4} {:<4} {:<6} {:<4} {}".format(
        "accession", "gene", "sup", "deg", "spec", "NA", "connects to known antigens"))
    for acc, c in ranked[:15]:
        conns = ", ".join(sorted(set(c["connects_to"])))
        na = ("yes" if c["na_associated"] else "no") if use_prior else "-"
        print("  {:<10} {:<11} {:<4} {:<4} {:<6} {:<4} {}".format(
            acc, c["gene"], c["partner_support"], c["total_partners"], c["specificity"], na, conns))

    # ---- OPTIONAL agent-assessment pilot (costs money) ----
    # `--assess N` runs the trusted, cited agent on the top N candidates to see what published
    # evidence actually exists for each. Off by default so the ranking stays free.
    assessments = []
    k = assess_count()
    if k is not None:
        items = [{"accession": acc, "gene": c["gene"],
                  "rationale": "physically associates with known {} autoantigens ({})".format(
                      disease_name, ", ".join(sorted(set(c["connects_to"]))))}
                 for acc, c in ranked[:k]]
        assessments = assess_top(disease, disease_name, items)

    payload = {
        "_metadata": {
            "project": "Autoantigen Discovery Platform",
            "disease_key": disease, "disease": disease_name,
            "phase": "7 - discovery mode (partner arm)",
            "run_date": datetime.date.today().isoformat(),
            "method": "partner arm: pool UniProt curated interaction partners of known autoantigens "
                      "(epitope-spreading prior); filter hubs (specificity>=min) and thin annotation "
                      "(degree>=min); rank by support (# known antigens linked), then specificity, then "
                      "the nuclear/NA prior. specificity = links-to-known-antigens / total partners.",
            "filters": {"min_support": MIN_SUPPORT, "min_degree": MIN_DEGREE,
                        "min_specificity": MIN_SPECIFICITY},
            "framing": "Candidates are CITED HYPOTHESES for expert review, not confirmed discoveries.",
            "validation": {"method": "leave-one-out hold-out rediscovery",
                           "discovery_recall": recall,
                           "recovered": [r[0] for r in recovered],
                           "not_recoverable_by_partners": missed},
            "n_known_antigens": len(positives(disease)),
            "n_novel_candidates_total": len(novel),
            "n_scored_shortlist": len(shortlist),
            "n_rediscovered_known": len(rediscovered_known),
            "caveat": "Partner arm finds complex-members well; standalone antigens (e.g. Ro60, "
                      "no curated interactions) are invisible to it — a known coverage gap. "
                      "Specificity controls for hub bias (e.g. HTT, 860 partners, correctly demoted).",
        },
        "ranking": {"order": "support, then specificity, then nuclear/NA prior",
                    "nuclear_na_prior_applied": use_prior,
                    "note": "Transparent multi-signal ranking; filters/order are a documented first pass."},
        "candidates": [
            {"accession": a, "gene": c["gene"], "partner_support": c["partner_support"],
             "total_partners": c["total_partners"], "specificity": c["specificity"],
             "na_associated": c["na_associated"],
             "connects_to": sorted(set(c["connects_to"]))}
            for a, c in ranked
        ],
        "agent_assessments": assessments,
    }
    out = os.path.join(PROCESSED_DIR, "discovery_{}.json".format(disease))
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)
    print("\n  {} scored candidates — full ranked list in data/processed/discovery_{}.json"
          .format(len(shortlist), disease))


if __name__ == "__main__":
    main()
