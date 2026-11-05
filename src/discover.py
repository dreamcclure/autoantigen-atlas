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
    "The protein {gene} ({acc}) physically associates with known {disease_name} autoantigens "
    "({connects}). Using your tools, investigate whether there is any PUBLISHED evidence that "
    "{gene} is itself targeted by autoantibodies or implicated as an autoantigen in {disease_name} "
    "(or a related autoimmune disease). Give a short, cited assessment of what you find.\n\n"
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

HEADERS = {"User-Agent": "lupus-autoantigen-assistant (educational project)"}
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


def main():
    disease = sys.argv[1] if len(sys.argv) > 1 else "lupus"
    disease_name = DISEASES[disease]["name"]
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
    if "--assess" in sys.argv:
        i = sys.argv.index("--assess")
        k = int(sys.argv[i + 1]) if i + 1 < len(sys.argv) and sys.argv[i + 1].isdigit() else 5
        top = ranked[:k]
        print("\n--- Agent assessment of top {} candidates (live, cited) ---".format(len(top)))
        from agent import build_agent, ask
        agent = build_agent(disease)
        for acc, c in top:
            conns = ", ".join(sorted(set(c["connects_to"])))
            q = ASSESS_QUESTION.format(gene=c["gene"], acc=acc, disease_name=disease_name, connects=conns)
            try:
                answer = ask(agent, q)
            except Exception as error:
                answer = "AGENT ERROR: {}".format(error)
            m = EVIDENCE_RE.findall(answer or "")
            tier = m[-1].upper() if m else "UNPARSED"
            print("  {:<11} {:<10} -> EVIDENCE: {}".format(c["gene"], acc, tier))
            assessments.append({"accession": acc, "gene": c["gene"], "evidence_tier": tier,
                                "transcript": answer})

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
