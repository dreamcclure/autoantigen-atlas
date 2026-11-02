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

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
DATASET = os.path.join(PROJECT_ROOT, "data", "processed", "autoantigens.json")
OUT = os.path.join(PROJECT_ROOT, "data", "processed", "evaluation.json")

# A fixed pilot subset (3 real autoantigens + 2 controls) for cheap harness testing.
PILOT = ["P10155", "P62314", "P02745", "P04406", "P60709"]

# The blind question template. It names the protein but NEVER its label or role.
# The decision rule is explicit so the agent can classify (not prove a universal
# negative): "no established evidence found" maps to NO, and UNCERTAIN is reserved
# for genuinely conflicting evidence. This lets specificity actually be measured
# while keeping the agent fully honest about what it did and didn't find.
QUESTION = (
    "Is the protein {short} (gene {gene}) an established, recognized autoantigen in "
    "systemic lupus erythematosus (SLE)? Investigate using your tools and give a "
    "short, cited assessment of what you find.\n\n"
    "Then classify with a single final line in exactly one of these formats:\n"
    "VERDICT: YES        -- your tools returned evidence it IS an established SLE autoantigen\n"
    "VERDICT: NO         -- your tools returned no evidence establishing it as an SLE "
    "autoantigen (treat absence of established evidence as NO; you are classifying, "
    "not proving a universal negative)\n"
    "VERDICT: UNCERTAIN  -- ONLY if the tools returned genuinely conflicting evidence "
    "you cannot resolve"
)

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
    pilot = "--pilot" in sys.argv
    with open(DATASET) as f:
        records = json.load(f)["autoantigens"]
    if pilot:
        records = [r for r in records if r["accession"] in PILOT]

    scope = "PILOT ({} proteins)".format(len(records)) if pilot else "FULL ({} proteins)".format(len(records))
    print("Running evaluation — {}\n".format(scope))

    agent = build_agent()
    results = []
    counts = {k: 0 for k in
              ["true_positive", "false_negative", "true_negative", "false_positive",
               "abstain_positive", "abstain_negative"]}

    for rec in records:
        short, gene, acc, role = rec["short_name"], rec["gene"], rec["accession"], rec["role"]
        question = QUESTION.format(short=short, gene=gene)
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
            "project": "Lupus Autoantigen Discovery Assistant",
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
    with open(OUT, "w") as f:
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
    print("\n  Full transcripts saved to data/processed/evaluation.json")


if __name__ == "__main__":
    main()
