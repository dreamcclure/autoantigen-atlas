"""
Phase 6 — Structural Analysis
=============================
Ask a real research question of our own data: do the established autoantigens have a
STRUCTURAL SIGNATURE that non-autoantigens lack? This is the "beyond retrieval ->
analysis" step — we stop just fetching facts and start interrogating the structures.

The lead comes from the literature (and our own Phase 1 note that several controls were
MORE ordered than real autoantigens): autoantigens may be enriched in intrinsic DISORDER
/ flexibility. AlphaFold's per-residue pLDDT is a well-established disorder proxy (low
pLDDT ~ disordered), so we can test the hypothesis on files we already downloaded.

Per protein we compute, from its AlphaFold model:
  * mean pLDDT                — overall confidence  (PRIMARY metric; higher = more ordered)
  * fraction disordered       — % residues pLDDT < 50 (AlphaFold "very low")
  * fraction low-confidence   — % residues pLDDT < 70
  * longest disordered run    — longest consecutive stretch of pLDDT < 50 (residues)

Two comparison groups are available:
  * default  — the disease's 5 curated negative controls (a quick specificity check)
  * --background — a RANDOM sample of ~40 reviewed human proteins (build_controls.py):
       the honest "typical protein" comparison, WITH length-confound controls, because
       disorder scales with length and we must show the signal is not just a size effect.

Stats are non-parametric and assumption-free (appropriate at this sample size):
  * difference in means tested by permutation (exact when small, else Monte-Carlo)
  * the length confound is (1) measured, (2) adjusted for with a length-STRATIFIED
    permutation test, and (3) confirmed with a length-MATCHED paired subset.

Run:
    python3 src/analyze_structure.py                        # lupus vs 5 curated controls
    python3 src/analyze_structure.py lupus --background     # lupus vs proteome background (rigorous)
    python3 src/analyze_structure.py rheumatoid_arthritis --background
"""

import itertools
import json
import os
import sys
import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec

from diseases import DISEASES

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
PROCESSED_DIR = os.path.join(PROJECT_ROOT, "data", "processed")
FIG_DIR = os.path.join(PROJECT_ROOT, "figures")

DISORDER_CUTOFF = 50
LOW_CONF_CUTOFF = 70
PRIMARY_METRIC = "mean_plddt"   # the one metric the confound-adjustment focuses on (least arbitrary)

EXACT_COMBO_CAP = 200_000       # enumerate exactly below this many splits, else Monte-Carlo
MC_SHUFFLES = 20_000            # Monte-Carlo permutation count
RNG_SEED = 42                   # reproducible permutations / matching

METRICS = [
    ("mean_plddt", "Mean pLDDT", "higher = more ordered"),
    ("frac_disordered", "Fraction disordered (pLDDT<50)", "higher = more disordered"),
    ("frac_low_conf", "Fraction low-conf (pLDDT<70)", "higher = more disordered"),
    ("longest_disordered_run", "Longest disordered run (aa)", "higher = more disordered"),
]

# group colors: autoantigen (blue), control/background (orange)
C_AUTO, C_CTRL = "#65CBF3", "#FF7D45"


# --------------------------------------------------------------------------- metrics
def residue_plddt(accession):
    """Per-residue pLDDT list (CA atoms) from a raw AlphaFold .pdb, or None if absent."""
    path = os.path.join(RAW_DIR, "{}_alphafold.pdb".format(accession))
    if not os.path.exists(path):
        return None
    scores = []
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                scores.append(float(line[60:66]))
    return scores or None


def longest_run_below(scores, cutoff):
    best = run = 0
    for s in scores:
        run = run + 1 if s < cutoff else 0
        best = max(best, run)
    return best


def metrics_for(scores):
    arr = np.array(scores, dtype=float)
    return {
        "length_aa": int(arr.size),
        "mean_plddt": round(float(arr.mean()), 1),
        "frac_disordered": round(float((arr < DISORDER_CUTOFF).mean()), 3),
        "frac_low_conf": round(float((arr < LOW_CONF_CUTOFF).mean()), 3),
        "longest_disordered_run": longest_run_below(scores, DISORDER_CUTOFF),
    }


# ----------------------------------------------------------------------------- stats
def permutation_pvalue(pos_vals, neg_vals, rng):
    """
    Two-sided permutation p-value for the difference in group MEANS.

    Exact (enumerate every split) when the number of splits is small; otherwise a
    Monte-Carlo permutation with MC_SHUFFLES random relabelings. No distributional
    assumption — honest at any sample size.
    """
    pooled = np.array(pos_vals + neg_vals, dtype=float)
    n_total, n_neg = len(pooled), len(neg_vals)
    observed = abs(float(np.mean(pos_vals)) - float(np.mean(neg_vals)))

    from math import comb
    n_combos = comb(n_total, n_neg)
    if n_combos <= EXACT_COMBO_CAP:
        extreme = count = 0
        for neg_combo in itertools.combinations(range(n_total), n_neg):
            mask = np.zeros(n_total, dtype=bool)
            mask[list(neg_combo)] = True
            diff = abs(pooled[~mask].mean() - pooled[mask].mean())
            extreme += diff >= observed - 1e-9
            count += 1
        return round(extreme / count, 4), observed, "exact ({} splits)".format(count)

    extreme = 0
    for _ in range(MC_SHUFFLES):
        perm = rng.permutation(n_total)
        neg = pooled[perm[:n_neg]]
        pos = pooled[perm[n_neg:]]
        if abs(pos.mean() - neg.mean()) >= observed - 1e-9:
            extreme += 1
    # +1/+1 (add-one) smoothing so a p-value is never reported as exactly 0.
    return round((extreme + 1) / (MC_SHUFFLES + 1), 4), observed, "monte-carlo ({})".format(MC_SHUFFLES)


def length_strata(lengths, n_bins=3):
    """Assign each protein to a length bin by pooled quantiles (so bins are balanced)."""
    lengths = np.asarray(lengths, dtype=float)
    edges = np.quantile(lengths, np.linspace(0, 1, n_bins + 1))
    edges[0], edges[-1] = -np.inf, np.inf   # open the ends so every point lands in a bin
    return np.digitize(lengths, edges[1:-1]), edges


def stratified_permutation(metric_vals, is_auto, lengths, rng, n_bins=3):
    """
    Length-ADJUSTED test: compare autoantigens vs controls WITHIN length bins, so a raw
    size difference can't masquerade as a disorder difference.

    Statistic = size-weighted average of (auto_mean - ctrl_mean) within each bin that
    has both groups. We permute labels WITHIN each bin to build the null — holding the
    length distribution fixed by construction.
    """
    metric_vals = np.asarray(metric_vals, dtype=float)
    is_auto = np.asarray(is_auto, dtype=bool)
    strata, _ = length_strata(lengths, n_bins)

    def stat(auto_mask):
        num = den = 0.0
        for b in np.unique(strata):
            sel = strata == b
            a = metric_vals[sel & auto_mask]
            c = metric_vals[sel & ~auto_mask]
            if len(a) and len(c):
                w = sel.sum()
                num += w * (a.mean() - c.mean())
                den += w
        return num / den if den else 0.0

    observed = stat(is_auto)
    extreme = 0
    for _ in range(MC_SHUFFLES):
        shuffled = is_auto.copy()
        for b in np.unique(strata):
            idx = np.where(strata == b)[0]
            shuffled[idx] = rng.permutation(shuffled[idx])
        if abs(stat(shuffled)) >= abs(observed) - 1e-9:
            extreme += 1
    p = round((extreme + 1) / (MC_SHUFFLES + 1), 4)
    return round(float(observed), 3), p, int(n_bins)


def matched_subset(auto, ctrl, metric_key, rng):
    """
    Length-MATCHED robustness check: pair each autoantigen with the nearest-length
    (unused) control, then test the paired differences by sign-flipping (a non-parametric
    paired test). Confirms the effect survives when sizes are matched one-for-one.
    """
    ctrl_pool = sorted(ctrl, key=lambda r: r["length_aa"])
    used = [False] * len(ctrl_pool)
    pairs = []
    for a in sorted(auto, key=lambda r: r["length_aa"]):
        best_i, best_d = None, None
        for i, c in enumerate(ctrl_pool):
            if used[i]:
                continue
            d = abs(c["length_aa"] - a["length_aa"])
            if best_d is None or d < best_d:
                best_i, best_d = i, d
        if best_i is None:
            break
        used[best_i] = True
        pairs.append((a[metric_key], ctrl_pool[best_i][metric_key],
                      abs(a["length_aa"] - ctrl_pool[best_i]["length_aa"])))

    diffs = np.array([p[0] - p[1] for p in pairs], dtype=float)
    observed = float(diffs.mean())
    extreme = 0
    for _ in range(MC_SHUFFLES):
        signs = rng.choice([-1, 1], size=len(diffs))
        if abs(float((diffs * signs).mean())) >= abs(observed) - 1e-9:
            extreme += 1
    p = round((extreme + 1) / (MC_SHUFFLES + 1), 4)
    return {
        "n_pairs": len(pairs),
        "median_length_gap_aa": float(np.median([p[2] for p in pairs])),
        "mean_paired_diff": round(observed, 3),
        "permutation_p": p,
    }


# ---------------------------------------------------------------------------- loading
def load_proteins(disease, use_background):
    """Return (autoantigens, controls, control_label) each as lists of metric dicts."""
    with open(os.path.join(PROCESSED_DIR, "{}.json".format(disease))) as f:
        records = json.load(f)["autoantigens"]
    auto_seed = [r for r in records if r["role"] == "positive"]

    if use_background:
        with open(os.path.join(PROCESSED_DIR, "proteome_background.json")) as f:
            ctrl_seed = json.load(f)["background"]
        label = "proteome background"
    else:
        ctrl_seed = [r for r in records if r["role"] == "negative_control"]
        label = "curated controls"

    def enrich(seed):
        out, skipped = [], []
        for r in seed:
            scores = residue_plddt(r["accession"])
            if scores is None:
                skipped.append(r["short_name"])
                continue
            m = metrics_for(scores)
            m.update({"short_name": r["short_name"], "accession": r["accession"],
                      "gene": r.get("gene", "")})
            out.append(m)
        return out, skipped

    auto, sk1 = enrich(auto_seed)
    ctrl, sk2 = enrich(ctrl_seed)
    return auto, ctrl, label, sk1 + sk2


# ----------------------------------------------------------------------------- figure
def build_figure(title, auto, ctrl, ctrl_label, out_path, with_scatter):
    n = len(METRICS)
    if with_scatter:
        fig = plt.figure(figsize=(4 * n, 9))
        gs = gridspec.GridSpec(2, n, height_ratios=[1, 1.15], hspace=0.32, wspace=0.28)
        strip_axes = [fig.add_subplot(gs[0, i]) for i in range(n)]
        scatter_ax = fig.add_subplot(gs[1, :])
    else:
        fig, strip_axes = plt.subplots(1, n, figsize=(4 * n, 5))
        scatter_ax = None
    fig.patch.set_facecolor("#0A0B10")

    for ax, (key, label, _) in zip(strip_axes, METRICS):
        ax.set_facecolor("#0A0B10")
        groups = [[r[key] for r in auto], [r[key] for r in ctrl]]
        for x, (vals, c) in enumerate(zip(groups, [C_AUTO, C_CTRL]), 1):
            jitter = np.linspace(-0.09, 0.09, len(vals)) if len(vals) > 1 else [0]
            ax.scatter([x + j for j in jitter], vals, c=c, s=34, alpha=0.8,
                       edgecolors="none", zorder=3)
            ax.hlines(np.median(vals), x - 0.22, x + 0.22, color=c, linewidth=2.4, zorder=2)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["autoantigen\n(n={})".format(len(auto)),
                            "{}\n(n={})".format(ctrl_label.split()[-1], len(ctrl))],
                           color="#E6E9F2", fontsize=9)
        ax.set_title(label, color="#E6E9F2", fontsize=10, pad=8)
        ax.tick_params(colors="#8A91A8")
        for sp in ax.spines.values():
            sp.set_color("#232838")

    if scatter_ax is not None:
        # The money plot: disorder (mean pLDDT) vs length. If autoantigens sit LOWER across
        # the whole length range, the signal is disorder — not just a size difference.
        scatter_ax.set_facecolor("#0A0B10")
        for grp, c, lbl in [(auto, C_AUTO, "autoantigen"), (ctrl, C_CTRL, ctrl_label)]:
            xs = [r["length_aa"] for r in grp]
            ys = [r["mean_plddt"] for r in grp]
            scatter_ax.scatter(xs, ys, c=c, s=46, alpha=0.82, edgecolors="none", label=lbl)
        scatter_ax.set_xscale("log")
        scatter_ax.set_xlabel("protein length (aa, log scale)", color="#E6E9F2", fontsize=10)
        scatter_ax.set_ylabel("mean pLDDT  (lower = more disordered)", color="#E6E9F2", fontsize=10)
        scatter_ax.set_title("Mean pLDDT vs protein length, by group "
                             "(do the clouds separate, or just overlap?)",
                             color="#E6E9F2", fontsize=11, pad=8)
        scatter_ax.tick_params(colors="#8A91A8")
        for sp in scatter_ax.spines.values():
            sp.set_color("#232838")
        leg = scatter_ax.legend(facecolor="#141824", edgecolor="#232838",
                                labelcolor="#E6E9F2", fontsize=10, loc="lower right")
        leg.get_frame().set_alpha(0.9)

    fig.suptitle(title, color="#E6E9F2", fontsize=13, y=0.995)
    os.makedirs(FIG_DIR, exist_ok=True)
    fig.savefig(out_path, dpi=150, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------- main
def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    disease = args[0] if args else "lupus"
    use_background = "--background" in sys.argv
    meta_name = DISEASES[disease]["name"]
    rng = np.random.default_rng(RNG_SEED)

    auto, ctrl, ctrl_label, skipped = load_proteins(disease, use_background)
    mode = "proteome background" if use_background else "curated controls"
    print("Structural analysis — {} vs {}  ({})\n".format(meta_name, ctrl_label, mode))
    if skipped:
        print("  skipped (no AlphaFold model): {}\n".format(", ".join(skipped)))

    # ---- per-metric comparison (all four), permutation test on the mean difference ----
    print("--- Autoantigen vs {} (median [min-max]) ---".format(ctrl_label))
    print("  {:<34} {:<20} {:<20} {:<9} {}".format("metric", "autoantigen", "control", "perm p", "method"))
    comparison = {}
    for key, label, direction in METRICS:
        av = [r[key] for r in auto]
        cv = [r[key] for r in ctrl]
        p, observed, method = permutation_pvalue(av, cv, rng)
        comparison[key] = {
            "label": label, "direction": direction,
            "autoantigen_median": round(float(np.median(av)), 3),
            "control_median": round(float(np.median(cv)), 3),
            "mean_difference": round(float(np.mean(av) - np.mean(cv)), 3),
            "permutation_p": p, "method": method,
            "primary": key == PRIMARY_METRIC,
        }
        f = lambda v: "{:.3f}".format(v) if isinstance(v, float) else str(v)
        star = " *primary" if key == PRIMARY_METRIC else ""
        print("  {:<34} {:<20} {:<20} {:<9} {}{}".format(
            label,
            "{} [{}-{}]".format(f(np.median(av)), f(min(av)), f(max(av))),
            "{} [{}-{}]".format(f(np.median(cv)), f(min(cv)), f(max(cv))),
            p, method, star))

    # ---- the length-confound analysis (only meaningful against the random background) ----
    confound = None
    if use_background:
        print("\n--- Length confound (does size differ, and does disorder survive adjusting for it?) ---")
        all_len = [r["length_aa"] for r in auto] + [r["length_aa"] for r in ctrl]
        lp, lobs, lmethod = permutation_pvalue(
            [r["length_aa"] for r in auto], [r["length_aa"] for r in ctrl], rng)
        print("  Length (median aa)    autoantigen {}  vs  control {}   perm p={} ({})".format(
            int(np.median([r["length_aa"] for r in auto])),
            int(np.median([r["length_aa"] for r in ctrl])), lp, lmethod))

        metric_vals = [r[PRIMARY_METRIC] for r in auto] + [r[PRIMARY_METRIC] for r in ctrl]
        is_auto = [True] * len(auto) + [False] * len(ctrl)
        strat_diff, strat_p, n_bins = stratified_permutation(metric_vals, is_auto, all_len, rng)
        print("  Length-STRATIFIED test on {} ({} length bins): adj. diff {}  perm p={}".format(
            PRIMARY_METRIC, n_bins, strat_diff, strat_p))

        matched = matched_subset(auto, ctrl, PRIMARY_METRIC, rng)
        print("  Length-MATCHED subset ({} pairs, median gap {:.0f} aa): paired diff {}  perm p={}".format(
            matched["n_pairs"], matched["median_length_gap_aa"],
            matched["mean_paired_diff"], matched["permutation_p"]))

        confound = {
            "length_median_autoantigen": int(np.median([r["length_aa"] for r in auto])),
            "length_median_control": int(np.median([r["length_aa"] for r in ctrl])),
            "length_permutation_p": lp,
            "stratified": {"primary_metric": PRIMARY_METRIC, "n_bins": n_bins,
                           "adjusted_mean_difference": strat_diff, "permutation_p": strat_p},
            "length_matched": matched,
        }

    # ---- figure ----
    suffix = "_background" if use_background else ""
    fig_path = os.path.join(FIG_DIR, "{}_structural_signature{}.png".format(disease, suffix))
    title = "{} · structural signature: autoantigens vs {}\n(horizontal bar = group median)".format(
        meta_name, ctrl_label)
    build_figure(title, auto, ctrl, ctrl_label, fig_path, with_scatter=use_background)

    # ---- report ----
    payload = {
        "_metadata": {
            "project": "Autoantigen Discovery Platform",
            "disease_key": disease, "disease": meta_name,
            "phase": "6 - structural analysis",
            "comparison_group": ctrl_label,
            "run_date": datetime.date.today().isoformat(),
            "hypothesis": "Autoantigens are enriched in intrinsic disorder (low pLDDT) "
                          "relative to non-autoantigens.",
            "disorder_proxy": "AlphaFold pLDDT (<50 very low ~ disordered; <70 low confidence)",
            "primary_metric": PRIMARY_METRIC,
            "stat_test": "permutation (exact when small, else Monte-Carlo); "
                         "length confound handled by stratified + matched tests",
            "rng_seed": RNG_SEED,
            "n_autoantigen": len(auto), "n_control": len(ctrl),
            "skipped_no_model": skipped,
            "caveat": ("Curated-control mode uses only 5 controls (exploratory). "
                       "Background mode uses a random proteome sample for real power; a "
                       "random draw may rarely include an unknown autoantigen, which is "
                       "infrequent and only makes the test more conservative."),
        },
        "comparison": comparison,
        "length_confound": confound,
        "per_protein": {
            "autoantigen": sorted(auto, key=lambda r: -r["frac_disordered"]),
            "control": sorted(ctrl, key=lambda r: -r["frac_disordered"]),
        },
    }
    out = os.path.join(PROCESSED_DIR, "structural_analysis_{}{}.json".format(disease, suffix))
    with open(out, "w") as f:
        json.dump(payload, f, indent=2)

    print("\n  Figure : figures/{}_structural_signature{}.png".format(disease, suffix))
    print("  Report : data/processed/structural_analysis_{}{}.json".format(disease, suffix))
    if not use_background:
        print("\n  NOTE: {} curated controls is a small sample — exploratory. Re-run with "
              "--background for the powered, length-controlled analysis.".format(len(ctrl)))


if __name__ == "__main__":
    main()
