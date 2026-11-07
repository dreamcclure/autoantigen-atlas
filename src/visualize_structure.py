"""
Visualize an AlphaFold structure, colored by confidence (pLDDT).
=================================================================
Reads a protein's AlphaFold .pdb file (already in data/raw/ from Phase 0) and draws
its backbone as a 3D trace, colored with AlphaFold's official confidence palette.
Saves a PNG into figures/.

Run:
    python3 src/visualize_structure.py                 # defaults to Ro60 (P10155)
    python3 src/visualize_structure.py P19474 Ro52     # any accession we've fetched
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")  # render to a file, no screen needed
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (enables 3D projection)

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(HERE)
RAW_DIR = os.path.join(PROJECT_ROOT, "data", "raw")
FIG_DIR = os.path.join(PROJECT_ROOT, "figures")

# AlphaFold's official pLDDT color bands.
BANDS = [
    (90, "#0053D6", "Very high (>90)"),
    (70, "#65CBF3", "Confident (70-90)"),
    (50, "#FFDB13", "Low (50-70)"),
    (0,  "#FF7D45", "Very low (<50)"),
]


def color_for(plddt):
    for threshold, color, _ in BANDS:
        if plddt > threshold:
            return color
    return BANDS[-1][1]


def read_ca_atoms(pdb_path):
    """Pull the backbone alpha-carbons: (x, y, z, pLDDT) per residue."""
    xs, ys, zs, bs = [], [], [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
                bs.append(float(line[60:66]))  # B-factor column = pLDDT
    return xs, ys, zs, bs


def render(acc, label):
    """Render one accession's AlphaFold structure to figures/{acc}_alphafold.png.
    Returns the output path, or None if we have no structure file for it. Factored out of
    main() so a batch job can render many proteins in one process (one matplotlib import)."""
    pdb_path = os.path.join(RAW_DIR, "{}_alphafold.pdb".format(acc))
    if not os.path.exists(pdb_path):
        return None

    xs, ys, zs, bs = read_ca_atoms(pdb_path)
    colors = [color_for(b) for b in bs]
    mean_plddt = round(sum(bs) / len(bs), 1)

    fig = plt.figure(figsize=(9, 8))
    fig.patch.set_facecolor("#0A0B10")            # dark "instrument" ground
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("#0A0B10")

    # Backbone drawn segment-by-segment so each bit carries its confidence color,
    # then the residue nodes on top.
    for i in range(len(xs) - 1):
        ax.plot(xs[i:i + 2], ys[i:i + 2], zs[i:i + 2],
                color=colors[i], linewidth=1.6, alpha=0.75, zorder=1)
    ax.scatter(xs, ys, zs, c=colors, s=16, depthshade=True, edgecolors="none", zorder=2)

    ax.set_title("{} ({}) · AlphaFold structure coloured by confidence (pLDDT)\n"
                 "mean pLDDT {} / 100".format(label, acc, mean_plddt),
                 color="#E6E9F2", fontsize=12, pad=14)
    ax.set_axis_off()  # clean hero — no bounding box or axis ticks
    leg = ax.legend(handles=[Line2D([0], [0], marker="o", linestyle="", markerfacecolor=c,
                                    markeredgecolor="none", markersize=8, label=lbl)
                             for _, c, lbl in BANDS],
                    title="AlphaFold confidence", loc="upper left", fontsize=9,
                    facecolor="#141824", edgecolor="#232838", labelcolor="#E6E9F2")
    leg.get_title().set_color("#8A91A8")

    os.makedirs(FIG_DIR, exist_ok=True)
    out = os.path.join(FIG_DIR, "{}_alphafold.png".format(acc))
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)   # free the figure — important when rendering many in one process
    return out


def main():
    acc = sys.argv[1] if len(sys.argv) > 1 else "P10155"
    label = sys.argv[2] if len(sys.argv) > 2 else "Ro60"
    out = render(acc, label)
    print("Saved: {}".format(out) if out else "No structure file for {}".format(acc))


if __name__ == "__main__":
    main()
