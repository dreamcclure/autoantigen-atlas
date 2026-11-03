"""
Phase 3 — Evaluation
====================
Score the agent against our golden set (24 autoantigens + 5 negative controls).

For each protein we ask the agent — BLIND (it is never told the right answer or
that a seed list exists) — whether it's an established lupus autoantigen. We parse
its VERDICT, compare to the known label, and compute recall, specificity, precision,
and the full confusion matrix. Every transcript is saved so each verdict is auditable.

Run a quick pilot first:
    python3 src/evaluate.py --pilot        # 5 proteins, well under $1
Then the full run:
    python3 src/evaluate.py                # all 29 (~a few dollars)
"""

import json
import os
import re
import sys
import datetime

from agent import build_agent, ask
from diseases import DISEASES

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")

# The blind question template. It names the protein but NEVER its label or role.
# The decision rule is explicit so the agent can classify (not prove a universal
# negative): "no established evidence found" maps to NO, and UNCERTAIN is reserved
# for genuinely conflicting evidence. This lets specificity actually be measured
# while keeping the agent fully honest about what it did and didn't find.
QUESTION = (
    "Is the protein {short} (gene {gene}) an established, recognized autoantigen in "
    "{disease_name}? Investigate using your tools and give a short, cited assessment "
    "of what you find.\n\n"
    "Then classify with a single final line in exactly one of these formats:\n"
    "VERDICT: YES        -- your tools returned evidence it IS an established autoantigen in {disease_name}\n"
    "VERDICT: NO         -- your tools returned no evidence establishing it as one "
    "(treat absence of established evidence as NO; you are classifying, not proving a "
    "universal negative)\n"
    "VERDICT: UNCERTAIN  -- ONLY if the tools returned genuinely conflicting evidence "
    "you cannot resolve"
)


def pilot_subset(records):
    """A cheap harness-test subset: the first 3 positives + first 2 controls."""
    positives = [r for r in records if r["role"] == "positive"][:3]
    controls = [r for r in records if r["role"] == "negative_control"][:2]
    return positives + controls

VERDICT_RE = re.compile(r"VERDICT:\s*(YES|NO|UNCERTAIN)", re.IGNORECASE)


def parse_verdict(text):
    """Pull the agent's final YES/NO/UNCERTAIN out of its reply (last match wins)."""
    matches = VERDICT_RE.findall(text or "")
    return matches[-1].upper() if matches else "UNPARSED"


def outcome_for(role, verdict):
    """Compare a verdict to the known label and name the outcome."""
    is_positive = (role == "positive")
    if verdict == "YES":
        return "true_positive" if is_positive else "false_positive"
    if verdict == "NO":
        return "false_negative" if is_positive else "true_negative"
    # UNCERTAIN or UNPARSED — an honest abstention, not a positive/negative call
    return "abstain_positive" if is_positive else "abstain_negative"


def main():
    # Disease = first non-flag argument (default lupus); --pilot = cheap subset.
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    disease = args[0] if args else "lupus"
    pilot = "--pilot" in sys.argv

    dataset = os.path.join(PROCESSED_DIR, "{}.json".format(disease))
    out = os.path.join(PROCESSED_DIR, "evaluation_{}.json".format(disease))
    disease_name = DISEASES[disease]["name"]

    with open(dataset) as f:
        records = json.load(f)["autoantigens"]
    if pilot:
        records = pilot_subset(records)

    scope = "PILOT ({} proteins)".format(len(records)) if pilot else "FULL ({} proteins)".format(len(records))
    print("Running evaluation for '{}' — {}\n".format(disease, scope))

    agent = build_agent(disease)
    results = []
    counts = {k: 0 for k in
              ["true_positive", "false_negative", "true_negative", "false_positive",
               "abstain_positive", "abstain_negative"]}

    for rec in records:
        short, gene, acc, role = rec["short_name"], rec["gene"], rec["accession"], rec["role"]
        question = QUESTION.format(short=short, gene=gene, disease_name=disease_name)
        try:
            answer = ask(agent, question)
        except Exception as error:
            answer = "AGENT ERROR: {}".format(error)
        verdict = parse_verdict(answer)
        outcome = outcome_for(role, verdict)
        counts[outcome] += 1

        mark = {"true_positive": "OK", "true_negative": "OK"}.get(outcome, "!!")
        print("  [{}] {:<16} role={:<16} verdict={:<10} -> {}".format(
            mark, short, role, verdict, outcome))

        results.append({
            "short_name": short, "gene": gene, "accession": acc,
            "true_role": role, "verdict": verdict, "outcome": outcome,
            "transcript": answer,
        })

    # --- metrics (UNCERTAIN/UNPARSED count as "not identified", reported separately) ---
    tp, fn = counts["true_positive"], counts["false_negative"]
    tn, fp = counts["true_negative"], counts["false_positive"]
    pos_total = tp + fn + counts["abstain_positive"]
    neg_total = tn + fp + counts["abstain_negative"]

    def rate(num, den):
        return round(num / den, 3) if den else None

    metrics = {
        "recall_sensitivity": rate(tp, pos_total),       # of real autoantigens, how many caught
        "specificity": rate(tn, neg_total),              # of controls, how many correctly rejected
        "precision": rate(tp, tp + fp),                  # of "YES" calls, how many were right
        "positives_total": pos_total,
        "negatives_total": neg_total,
        "counts": counts,
    }

    payload = {
        "_metadata": {
            "project": "Autoantigen Discovery Platform",
            "disease_key": disease,
            "disease": disease_name,
            "phase": "3 - evaluation",
            "scope": "pilot" if pilot else "full",
            "run_date": datetime.date.today().isoformat(),
            "model": "claude-opus-4-8",
            "blind": True,
            "note": "Agent was never told the label or that a seed list exists.",
        },
        "metrics": metrics,
        "results": results,
    }
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)

    # --- report ---
    print("\n--- Confusion matrix ---")
    print("                 agent: YES   agent: NO   abstain")
    print("  real autoantigen   {:>3}         {:>3}        {:>3}".format(
        tp, fn, counts["abstain_positive"]))
    print("  control (not)      {:>3}         {:>3}        {:>3}".format(
        fp, tn, counts["abstain_negative"]))
    print("\n--- Metrics ---")
    print("  Recall / sensitivity : {}  ({}/{} real autoantigens identified)".format(
        metrics["recall_sensitivity"], tp, pos_total))
    print("  Specificity          : {}  ({}/{} controls correctly rejected)".format(
        metrics["specificity"], tn, neg_total))
    print("  Precision            : {}".format(metrics["precision"]))
    print("\n  Full transcripts saved to data/processed/evaluation_{}.json".format(disease))


if __name__ == "__main__":
    main()
