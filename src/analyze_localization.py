"""
Phase 6 — Structural/biological signature, take 2: LOCALIZATION & NUCLEIC-ACID BINDING
=====================================================================================
The bulk-disorder hypothesis didn't survive a proper control (see analyze_structure.py).
This tests a hypothesis with a real mechanism behind it.

THE BIOLOGY (why we expect a signal): the classic Rosen / Casciola-Rosen model of lupus
holds that NUCLEAR antigens are cleaved during apoptosis and cluster in surface blebs,
where the immune system finally encounters them. The canonical lupus autoantigens are
components of nucleic-acid–protein complexes — nucleosomes (histones), the spliceosomal
snRNPs (Sm, U1), the Ro/La RNP. So we predict autoantigens are strongly enriched for
NUCLEAR localization and/or NUCLEIC-ACID BINDING vs a typical human protein.

WHY THIS ISN'T CIRCULAR: our autoantigen list was chosen by clinical autoantibody tests
(anti-Ro, anti-Sm, ...), NOT by where the proteins live. So enrichment for nuclear /
NA-binding is a genuine discovered property, not a consequence of how we picked the set.

HOW WE TEST IT: every protein already has a local UniProt record (data/raw/). UniProt
assigns controlled-vocabulary KEYWORDS; we classify each protein by three of them —
"Nucleus", "DNA-binding", "RNA-binding" — plus a combined primary endpoint (any of the
three). We compare autoantigens vs the random proteome background with FISHER'S EXACT
TEST (the exact test for a 2x2 proportion table; computed from the hypergeometric
distribution with math.comb — no scipy).

Run:
    python3 src/analyze_localization.py                      # lupus
    python3 src/analyze_localization.py rheumatoid_arthritis
"""

import json
import os
import sys
import datetime
from math import comb

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diseases import DISEASES

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
FIG_DIR = os.path.join(PROJECT_ROOT, "figures")

# UniProt keyword name -> our short label. These are controlled-vocabulary terms.
KEYWORD_PROPS = ["Nucleus", "DNA-binding", "RNA-binding"]

C_AUTO, C_CTRL = "#65CBF3", "#FF7D45"


def keywords_for(accession):
    """Return the set of UniProt keyword names for one protein (from its raw JSON)."""
    path = os.path.join(RAW_DIR, "{}_uniprot.json".format(accession))
    if not os.path.exists(path):
        return None
    with open(path) as f:
        record = json.load(f)
    return {k.get("name") for k in record.get("keywords", [])}


def classify(accession):
    """Boolean flags for each property + the combined 'nucleic-acid-associated'."""
    kws = keywords_for(accession)
    if kws is None:
        return None
    flags = {prop: (prop in kws) for prop in KEYWORD_PROPS}
    flags["nucleic_acid_associated"] = any(flags.values())   # primary endpoint
    return flags


def fisher_exact_two_sided(a, b, c, d):
    """
    Two-sided Fisher's exact test for the 2x2 table
        [[a, b],   a = autoantigen & has property,  b = autoantigen & lacks
         [c, d]]   c = background  & has property,  d = background  & lacks

    p = sum of hypergeometric probabilities of all tables (same margins) whose
    probability is <= the observed table's. Odds ratio = a*d / b*c.
    """
    n1, n2 = a + b, c + d          # row totals (group sizes)
    k, N = a + c, a + b + c + d    # col total (# with property), grand total

    def hyper(x):
        # P(top-left cell == x), hypergeometric with these margins
        lo = max(0, n1 - (N - k))
        hi = min(n1, k)
        if x < lo or x > hi:
            return 0.0
        return comb(k, x) * comb(N - k, n1 - x) / comb(N, n1)

    p_obs = hyper(a)
    lo = max(0, n1 - (N - k))
    hi = min(n1, k)
    p = sum(hyper(x) for x in range(lo, hi + 1) if hyper(x) <= p_obs + 1e-12)

    odds = (a * d) / (b * c) if (b * c) else float("inf")
    return round(min(p, 1.0), 4), (round(odds, 2) if odds != float("inf") else None)


def load_group(accessions):
    flags, skipped = [], []
    for acc in accessions:
        f = classify(acc)
        (skipped if f is None else flags).append(acc if f is None else f)
    return flags, skipped


def build_figure(disease_name, props, auto_n, ctrl_n, results, out_path):
    labels = [p["label"] for p in props]
    auto_fracs = [p["auto_has"] / auto_n for p in props]
    ctrl_fracs = [p["ctrl_has"] / ctrl_n for p in props]

    x = np.arange(len(labels))
    w = 0.38
    fig, ax = plt.subplots(figsize=(1.9 * len(labels) + 3, 6))
    fig.patch.set_facecolor("#0A0B10")
    ax.set_facecolor("#0A0B10")

    ax.bar(x - w / 2, auto_fracs, w, color=C_AUTO, label="autoantigen (n={})".format(auto_n))
    ax.bar(x + w / 2, ctrl_fracs, w, color=C_CTRL, label="proteome background (n={})".format(ctrl_n))

    for xi, p in zip(x, props):
        top = max(p["auto_has"] / auto_n, p["ctrl_has"] / ctrl_n)
        star = "p={}".format(p["p"]) if p["p"] >= 0.001 else "p<0.001"
        ax.text(xi, top + 0.03, star, ha="center", color="#E6E9F2", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, color="#E6E9F2", fontsize=9)
    ax.set_ylabel("fraction of proteins with the property", color="#E6E9F2", fontsize=10)
    ax.set_ylim(0, 1.05)
    ax.set_title("{} · are autoantigens enriched for nuclear / nucleic-acid biology?\n"
                 "(UniProt keywords; Fisher's exact test vs a random proteome background)"
                 .format(disease_name), color="#E6E9F2", fontsize=12, pad=12)
    ax.tick_params(colors="#8A91A8")
    for sp in ax.spines.values():
        sp.set_color("#232838")
    leg = ax.legend(facecolor="#141824", edgecolor="#232838", labelcolor="#E6E9F2", fontsize=10)
    leg.get_frame().set_alpha(0.9)

    os.makedirs(FIG_DIR, exist_ok=True)
    fig.savefig(out_path, dpi=150, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)


def main():
    disease = sys.argv[1] if len(sys.argv) > 1 else "lupus"
    disease_name = DISEASES[disease]["name"]

    with open(os.path.join(PROCESSED_DIR, "{}.json".format(disease))) as f:
        pos_acc = [r["accession"] for r in json.load(f)["autoantigens"] if r["role"] == "positive"]
    with open(os.path.join(PROCESSED_DIR, "proteome_background.json")) as f:
        bg_acc = [r["accession"] for r in json.load(f)["background"]]

    auto, sk1 = load_group(pos_acc)
    ctrl, sk2 = load_group(bg_acc)
    auto_n, ctrl_n = len(auto), len(ctrl)

    print("Localization / NA-binding — {}  (autoantigens n={} vs proteome background n={})\n"
          .format(disease_name, auto_n, ctrl_n))

    # primary endpoint first, then the three components
    endpoints = [
        ("nucleic_acid_associated", "Nuclear and/or\nNA-binding (primary)"),
        ("Nucleus", "Nucleus"),
        ("DNA-binding", "DNA-binding"),
        ("RNA-binding", "RNA-binding"),
    ]

    props = []
    print("  {:<34} {:<16} {:<16} {:<8} {}".format("endpoint", "autoantigen", "background", "OR", "Fisher p"))
    for key, label in endpoints:
        a = sum(1 for f in auto if f[key]);  b = auto_n - a
        c = sum(1 for f in ctrl if f[key]);  d = ctrl_n - c
        p, odds = fisher_exact_two_sided(a, b, c, d)
        props.append({"key": key, "label": label.replace("\n", " "),
                      "auto_has": a, "ctrl_has": c, "odds_ratio": odds, "p": p})
        or_s = "inf" if odds is None else "{:.1f}".format(odds)
        print("  {:<34} {:<16} {:<16} {:<8} {}{}".format(
            label.replace("\n", " "),
            "{}/{} ({:.0%})".format(a, auto_n, a / auto_n),
            "{}/{} ({:.0%})".format(c, ctrl_n, c / ctrl_n),
            or_s, p, "  *primary" if key == "nucleic_acid_associated" else ""))

    # figure uses the display labels (with line breaks)
    fig_props = []
    for (key, label), p in zip(endpoints, props):
        fig_props.append({**p, "label": label})
    fig_path = os.path.join(FIG_DIR, "{}_localization_signature.png".format(disease))
    build_figure(disease_name, fig_props, auto_n, ctrl_n, props, fig_path)

    payload = {
        "_metadata": {
            "project": "Autoantigen Discovery Platform",
            "disease_key": disease, "disease": disease_name,
            "phase": "6 - structural analysis (localization)",
            "run_date": datetime.date.today().isoformat(),
            "hypothesis": "Autoantigens are enriched for nuclear localization and/or "
                          "nucleic-acid binding vs a typical human protein (apoptotic-bleb model).",
            "classification": "UniProt controlled-vocabulary keywords: Nucleus, DNA-binding, RNA-binding",
            "primary_endpoint": "nucleic_acid_associated (any of the three keywords)",
            "stat_test": "Fisher's exact test, two-sided (hypergeometric; no scipy)",
            "non_circular": "Autoantigens were selected by clinical autoantibody tests, not by "
                            "localization, so this enrichment is a discovered property, not a selection loop.",
            "n_autoantigen": auto_n, "n_background": ctrl_n,
            "skipped_no_uniprot": sk1 + sk2,
            "caveat": "Keyword-based classification is coarse (e.g. Ro60 is keyworded cytoplasmic "
                      "but is an RNA-binding RNP) — the combined primary endpoint is the robust read.",
        },
        "results": props,
    }
    out = os.path.join(PROCESSED_DIR, "localization_analysis_{}.json".format(disease))
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)

    print("\n  Figure : figures/{}_localization_signature.png".format(disease))
    print("  Report : data/processed/localization_analysis_{}.json".format(disease))


if __name__ == "__main__":
    main()
