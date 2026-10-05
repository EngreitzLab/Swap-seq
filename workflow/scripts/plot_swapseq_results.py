#! /usr/bin/env python

"""Swap-seq screen-result plots.

  volcano.pdf                          per-edit effect vs BH significance, 95% CI bars
  horizontal_barplots/<element>.pdf    every edit tiling one element; red if BH < threshold
"""

import argparse
import os
import re

import matplotlib
matplotlib.use("Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import t

from plot_style import register_arial

EDIT_INDEX = "Swap-seq edit index (per edit)"
EDIT_NAME = "Swap-seq edit name (per edit)"
CATEGORY = "Design categories"
EFFECT = "mean_effect_percent"

COLOR_SIG = "#FF2426"
COLOR_INSIG = "gray"

CATEGORY_COLORS = {
    "PPIF promoter": "#90278E",
    "PPIF splice junction": "#FFC900",
    "CRISPRi hits": "#603913",
    "candidate enhancers": "#A97C50",
    "CTCF ChIP-seq peak": "#2BB673",
    "Other candidate elements": "#3386C9",
    "Negative control": "#58595B",
}
LEGEND_ORDER = ["PPIF promoter", "PPIF splice junction", "CRISPRi hits",
                "candidate enhancers", "CTCF ChIP-seq peak",
                "Other candidate elements", "Negative control"]
## Background to foreground.
DRAW_ORDER = ["Negative control", "Other candidate elements", "PPIF promoter",
              "PPIF splice junction", "CTCF ChIP-seq peak", "candidate enhancers",
              "CRISPRi hits"]

## design categories
CATEGORY_COLLAPSE = {
    "PPIF promoter": "PPIF promoter",
    "PPIF splice junction": "PPIF splice junction",
    "CRISPRi hits": "CRISPRi hits",
    "Candidate enhancer": "candidate enhancers",
    "CTCF ChIP-seq peaks": "CTCF ChIP-seq peak",
    "Other candidate elements": "Other candidate elements",
    "Negative controls": "Negative control",
}
DROP_CATEGORIES = {"Nearby-gene promoter"} # promoter controls at other genes

ELEMENT_ALIASES = [
    ("PPIFprom", "chr10:81106774-81107400_PPIF_promoter"),
    ("PPIFenh", "chr10:81045660-81046982_Nasser_enhancer_Peak88"),
]

MIN_EDITS_PER_ELEMENT = 2

BAR_ROW_PITCH = 0.5
BAR_HEIGHT = 0.62


def style():
    """Arial, ticks style, editable PDF text."""
    mpl.rcParams.update({
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 1.0,
    })
    register_arial()
    sns.set_style("ticks")
    sns.set_context("talk", font_scale=1)
    # seaborn's context/style calls reset the font family, so re-assert it afterwards.
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--edit-effects", required=True,
                        help="per-edit table (SwapseqEditEffects.tsv)")
    parser.add_argument("--replicate-effects", required=True,
                        help="per-edit x replicate table (SwapseqEditReplicateEffects.tsv)")
    parser.add_argument("--outdir", required=True, help="directory to write into")
    parser.add_argument("--bh-threshold", type=float, default=0.001,
                        help="threshold at which an edit is called significant")
    return parser.parse_args()


def bar_color(bh, threshold):
    return COLOR_SIG if np.isfinite(bh) and bh < threshold else COLOR_INSIG


def ci_95(sd, n):
    """Half-width of the 95% CI of the mean; 0 when n < 2. Scalars or Series."""
    n = np.asarray(n, dtype=float)
    sd = np.asarray(sd, dtype=float)
    half = np.where(n > 1, t.ppf(0.975, np.maximum(n - 1, 1)) * sd / np.sqrt(n), 0.0)
    return half if half.ndim else float(half)


def element_of(edit_name):
    """The element an edit tiles: `<element>_<editIndex>`, or via ELEMENT_ALIASES."""
    name = str(edit_name)
    for token, element in ELEMENT_ALIASES:
        if token in name:
            return element
    stripped = re.sub(r"(_newOverlap-\d+)?_\d+$", "", name)
    return stripped or name


def plot_volcano(edits, outfile, bh_threshold):
    data = edits.dropna(subset=[EFFECT, "BH_p_value", CATEGORY]).copy()
    data = data[~data[CATEGORY].isin(DROP_CATEGORIES)]
    data["class"] = data[CATEGORY].map(CATEGORY_COLLAPSE).fillna("Other candidate elements")
    # BH of exactly 0 is infinite on a log axis; floor at the smallest non-zero value
    positive = data.loc[data["BH_p_value"] > 0, "BH_p_value"]
    floor = positive.min() if len(positive) else 1e-300
    data["-log10p"] = -np.log10(data["BH_p_value"].clip(lower=floor))
    data["sig"] = data["BH_p_value"] < bh_threshold

    data["CI_95"] = ci_95(data["SD"], data["n_replicates"])

    x_max = np.abs(data[EFFECT]).max() * 1.1
    fig, ax = plt.subplots(figsize=(9, 6.5))

    for depth, name in enumerate(DRAW_ORDER):
        color = CATEGORY_COLORS[name]
        group = data[data["class"] == name]
        insignificant, significant = group[~group["sig"]], group[group["sig"]]
        if len(insignificant):
            ax.errorbar(insignificant[EFFECT], insignificant["-log10p"],
                        xerr=insignificant["CI_95"], fmt="o", color=color, ecolor=color,
                        elinewidth=0.7, capsize=0, zorder=depth, markersize=7, alpha=0.25,
                        markeredgewidth=0)
        if len(significant):
            ax.errorbar(significant[EFFECT], significant["-log10p"],
                        xerr=significant["CI_95"], fmt="o", color=color, ecolor=color,
                        elinewidth=1.1, capsize=0, zorder=depth + 10, markersize=9,
                        alpha=0.9, markeredgewidth=0)

    ax.axhline(-np.log10(bh_threshold), color="black", linestyle="--", linewidth=1)
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set_xlim(-x_max, x_max)
    ax.set_xlabel("Effect Size (%)", fontsize=14)
    ax.set_ylabel("-log10(BH p-value)", fontsize=14)
    ax.set_title(f"Swap-seq per-edit effects  ({len(data)} edits, "
                 f"{int(data['sig'].sum())} significant at BH < {bh_threshold:g})",
                 fontsize=12, pad=10)

    # Explicit handles so every category appears, including any with no hits.
    handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=9,
                      markeredgewidth=0, color=CATEGORY_COLORS[name], label=name)
               for name in LEGEND_ORDER]
    ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=10,
              borderaxespad=0, frameon=False)

    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  wrote {outfile} ({len(data)} edits, {int(data['sig'].sum())} significant)")


def plot_element_barplot(name, edits, reps, outfile, bh_threshold):
    """Horizontal barplot of every edit tiling one element."""
    order = edits.sort_values(EFFECT, ascending=True)
    labels, means, errors, colors, dots = [], [], [], [], []
    for _, row in order.iterrows():
        values = reps.loc[reps[EDIT_INDEX] == row[EDIT_INDEX], EFFECT].values
        labels.append(row[EDIT_NAME])
        means.append(row[EFFECT])
        # CI from the tested SD, so the bars agree with the volcano and the p-value
        errors.append(ci_95(row["SD"], row["n_replicates"]))
        colors.append(bar_color(row["BH_p_value"], bh_threshold))
        dots.append(values)

    fig, ax = plt.subplots(figsize=(10, max(2.5, len(labels) * BAR_ROW_PITCH + 1.4)))
    for i, (mean, error, color) in enumerate(zip(means, errors, colors)):
        ax.barh(i, mean, height=BAR_HEIGHT, xerr=error if error > 0 else None, color=color,
                edgecolor="none", linewidth=0,
                error_kw={"elinewidth": 1.5, "capsize": 0, "color": "black"})
    for i, values in enumerate(dots):
        if len(values):
            ax.scatter(values, [i] * len(values), color="black", alpha=0.25,
                       edgecolors="none", s=12, zorder=5)

    ax.set_yticks(range(len(labels)))
    ## Explicit sizes: the seaborn "talk" context would swamp the panel with long names.
    ax.set_yticklabels(labels, fontsize=9)
    ax.tick_params(axis="x", labelsize=10)
    ax.set_ylim(len(labels) - 0.5, -0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.5)
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Average Percent Effects (%)", fontsize=12)
    ax.set_ylabel("Variant Name", fontsize=12)
    ax.set_title(name, fontsize=11, pad=10)

    ax.legend(handles=[Patch(facecolor=COLOR_SIG, edgecolor="none",
                             label=f"BH < {bh_threshold:g}"),
                       Patch(facecolor=COLOR_INSIG, edgecolor="none",
                             label=f"BH $\\geq$ {bh_threshold:g}")],
              loc="best", frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig(outfile, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)


def safe_filename(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")[:120]


def plot_all_elements(edits, reps, outdir, bh_threshold):
    os.makedirs(outdir, exist_ok=True)
    real = edits[edits[EDIT_INDEX] != 0].copy()
    real["element"] = real[EDIT_NAME].map(element_of)

    written, pooled = 0, []
    for element, group in real.groupby("element"):
        if len(group) < MIN_EDITS_PER_ELEMENT:
            pooled.append(group)
            continue
        plot_element_barplot(element, group, reps,
                             os.path.join(outdir, safe_filename(element) + ".pdf"),
                             bh_threshold)
        written += 1
    if pooled:
        group = pd.concat(pooled)
        plot_element_barplot(f"Single-edit elements (n={len(group)})", group, reps,
                             os.path.join(outdir, "single_edit_elements.pdf"), bh_threshold)
        written += 1
        print(f"    single-edit elements: {sorted(group[EDIT_NAME])}")
    print(f"  wrote {written} barplots into {outdir}")


def main():
    args = parse_args()
    style()
    os.makedirs(args.outdir, exist_ok=True)

    edits = pd.read_table(args.edit_effects)
    reps = pd.read_table(args.replicate_effects)
    if "Excluded" in reps.columns:
        reps = reps[~reps["Excluded"].astype(str).str.upper().isin({"TRUE", "1"})]

    plot_volcano(edits, os.path.join(args.outdir, "volcano.pdf"), args.bh_threshold)
    plot_all_elements(edits, reps, os.path.join(args.outdir, "horizontal_barplots"),
                      args.bh_threshold)


if __name__ == "__main__":
    main()
