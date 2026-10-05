#! /usr/bin/env python

"""Swap-seq QC plots, drawn from the pipeline's own Swap-seq tables.

  replicate_correlation_biorep_{perEdit,perPegRNA}.pdf   BioRep 1 vs BioRep 2 scatter
  replicate_correlation_ffreps_{perEdit,perPegRNA}.pdf   FlowFISH replicate scatter
  replicate_correlation_{all,filtered}.pdf               FFrep correlation heatmap
  barcode_consistency_min{2,3,4}.pdf                     split-half barcode reproducibility
  pool_allele_frequency_by_biorep.pdf                    edited-allele fraction per BioRep
  pool_allele_frequency_distribution_combined.pdf        per-barcode allele frequency
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import combinations
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats

from plot_style import register_arial

## Column names as written by swapseq_edit_effects.py.
EDIT_INDEX = "Swap-seq edit index (per edit)"
EDIT_NAME = "Swap-seq edit name (per edit)"
BARCODE_NAME = "Deletion Pair Name (per barcode)"
REP_EFFECT = "mean_effect_percent"
PEG_EFFECT = "Swap-seq final effect size percent"
SPIKEIN_FREQ = "Swap-seq freq"

SPLIT_SEED = 42

BAR_COLOR = "#808080"
ALLELE_FREQ_THRESHOLD = 1e-5
SORT_BIN_COLUMNS = ["A", "B", "C", "D"]


def style():
    """House style for these figures: Arial, no top/right spines, editable PDF text."""
    register_arial()
    matplotlib.rcParams.update({
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 150,
    })


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--replicate-effects", required=True,
                        help="per-edit x replicate table (SwapseqEditReplicateEffects.tsv)")
    parser.add_argument("--pegrna-effects", required=True,
                        help="per-pegRNA table (SwapseqPegRNAEffects.tsv)")
    parser.add_argument("--pegrna-effects-unfiltered", default="",
                        help="per-pegRNA table before the min-allele-frequency cut "
                             "(SwapseqPegRNAEffects.unfiltered.tsv). Required for the "
                             "pool-editing panels, which report how many barcodes clear "
                             "that threshold.")
    parser.add_argument("--outdir", required=True, help="directory to write the PDFs into")
    return parser.parse_args()


def replicate_label(biorep, ffrep):
    return f"BioRep{biorep} FF{ffrep}"


DEFAULT_LIMITS = (-110, 100)


def scatter_limits(values, default=DEFAULT_LIMITS):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not len(finite):
        return default
    span = default[1] - default[0]
    low = min(default[0], finite.min() - 0.02 * span)
    high = max(default[1], finite.max() + 0.02 * span)
    return (low, high)


def limits_widened(limits, default=DEFAULT_LIMITS, tolerance=0.2):
    span = default[1] - default[0]
    return ((default[0] - limits[0]) > tolerance * span
            or (limits[1] - default[1]) > tolerance * span)


def excluded_replicate_pairs(reps):
    """The (BioRep, FFRep) pairs the Swap-seq step flagged as excluded."""
    if "Excluded" not in reps.columns:
        return set()
    flagged = reps[reps["Excluded"].astype(str).str.upper().isin({"TRUE", "1"})]
    return set(zip(flagged["BioRep"], flagged["FFRep"]))


def drop_excluded(table, excluded_pairs):
    if not excluded_pairs:
        return table
    keep = [(b, f) not in excluded_pairs
            for b, f in zip(table["BioRep"], table["FFRep"])]
    return table[np.array(keep, dtype=bool)]


def plot_replicate_correlation(reps, outfile, title):
    """Heatmap of per-edit effect correlation between every pair of replicates."""
    wide = reps.pivot_table(index=EDIT_INDEX, columns=["BioRep", "FFRep"], values=REP_EFFECT)
    if wide.shape[1] < 2:
        print(f"  skipping {outfile}: fewer than two replicates")
        return
    # Pairwise Pearson over edits measured in both replicates of each pair.
    corr = wide.corr(method="pearson", min_periods=3)
    labels = [replicate_label(b, f) for b, f in corr.columns]

    fig, ax = plt.subplots(figsize=(1 + 0.55 * len(labels), 0.9 + 0.55 * len(labels)))
    im = ax.imshow(corr.values, cmap="RdYlBu_r", vmin=0, vmax=1)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    for i in range(len(labels)):
        for j in range(len(labels)):
            value = corr.values[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=5,
                        color="white" if value > 0.75 or value < 0.25 else "black")
    mean_off_diagonal = np.nanmean(corr.values[~np.eye(len(labels), dtype=bool)])
    ax.set_title(f"{title}\nPearson r of per-edit effect between replicates "
                 f"(n={wide.shape[0]} edits, mean off-diagonal r={mean_off_diagonal:.2f})",
                 fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.7, label="Pearson r")
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    print(f"  wrote {outfile} (mean off-diagonal r={mean_off_diagonal:.3f})")


def plot_biorep_scatter(table, key, effect, excluded_pairs, outfile, level):
    """BioRep 1 vs BioRep 2 scatter."""
    data = drop_excluded(table, excluded_pairs)
    per_biorep = data.groupby([key, "BioRep"])[effect].mean().reset_index()
    wide = per_biorep.pivot(index=key, columns="BioRep", values=effect)
    bioreps = sorted(wide.columns)
    if len(bioreps) < 2:
        print(f"  skipping {outfile}: only one biological replicate")
        return
    paired = wide[bioreps[:2]].dropna()
    if len(paired) < 3:
        print(f"  skipping {outfile}: only {len(paired)} paired points")
        return
    x, y = paired[bioreps[0]], paired[bioreps[1]]
    r = stats.pearsonr(x, y)[0]

    fig, ax = plt.subplots(figsize=(5.6, 5.6))
    limits = scatter_limits(np.concatenate([x.values, y.values]))
    ax.plot(limits, limits, "k--", linewidth=1.2, zorder=1)
    ax.scatter(x, y, s=45, color="#4B85C5", edgecolor="black", linewidth=0.4, zorder=2)
    ax.text(0.05, 0.93, f"Pearson's r = {r:.2f}\nn = {len(paired)}",
            transform=ax.transAxes, fontsize=12, va="top")
    ax.set_xlim(*limits)
    ax.set_ylim(*limits)
    ax.set_xlabel(f"BioRep {bioreps[0]} effect (%)", fontsize=13)
    ax.set_ylabel(f"BioRep {bioreps[1]} effect (%)", fontsize=13)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"Biological replicate correlation, {level}\n"
                 "allele-frequency filtered, poorly correlating sorts excluded",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(outfile)
    plt.close(fig)
    print(f"  wrote {outfile} ({level}: n={len(paired)}, r={r:.3f})")


def plot_ffrep_scatter_pages(table, key, effect, excluded_pairs, outfile, level):
    """One scatter per FlowFISH replicate pair."""
    data = table.copy()
    data["Replicate"] = [replicate_label(b, f)
                         for b, f in zip(data["BioRep"], data["FFRep"])]
    excluded = {replicate_label(b, f) for b, f in excluded_pairs}

    wide = data.pivot_table(index=key, columns="Replicate", values=effect)
    ## Numeric order, not lexicographic: "BioRep1 FF10" must not sort before "FF2".
    order = {replicate_label(b, f): (b, f)
             for b, f in zip(data["BioRep"], data["FFRep"])}
    replicates = sorted(wide.columns, key=lambda name: order[name])
    if len(replicates) < 2:
        print(f"  skipping {outfile}: fewer than two replicates")
        return

    pages = []
    for a, b in combinations(replicates, 2):
        paired = wide[[a, b]].dropna()
        if len(paired) < 3:
            continue
        r = stats.pearsonr(paired[a], paired[b])[0]
        pages.append((r, a, b, paired))

    ## Limits per page: one bad sort should not stretch the axes of every other page.
    with PdfPages(outfile) as pdf:
        for r, a, b, paired in pages:
            limits = scatter_limits(paired.to_numpy().ravel())
            expanded = limits_widened(limits)
            fig, ax = plt.subplots(figsize=(5.2, 5.2))
            ax.plot(limits, limits, "k--", linewidth=1.2, zorder=1)
            ax.scatter(paired[a], paired[b], s=28, color="#4B85C5", edgecolor="black",
                       linewidth=0.3, zorder=2)
            ax.text(0.05, 0.93, f"Pearson's r = {r:.2f}\nn = {len(paired)}",
                    transform=ax.transAxes, fontsize=10, va="top")
            ax.set_xlim(*limits)
            ax.set_ylim(*limits)
            ax.set_xlabel(f"{a} effect (%)", fontsize=11)
            ax.set_ylabel(f"{b} effect (%)", fontsize=11)
            ax.set_aspect("equal", adjustable="box")
            tags = [name for name in (a, b) if name in excluded]
            note = f"  [excluded: {', '.join(tags)}]" if tags else ""
            if expanded:
                note += f"  [axes widened to {limits[0]:g}, {limits[1]:g}]"
            ax.set_title(f"{a}  vs  {b}{note}\n{level}", fontsize=10)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    if pages:
        worst = min(pages, key=lambda page: page[0])
        print(f"  wrote {outfile} ({level}: {len(pages)} pages in replicate order; "
              f"lowest pair {worst[1]} vs {worst[2]} r={worst[0]:.2f})")
    else:
        print(f"  wrote {outfile} (no pages)")


def plot_barcode_consistency(pegrnas, excluded_pairs, outdir, thresholds=(2, 3, 4)):
    """Split-half barcode consistency."""
    data = drop_excluded(pegrnas, excluded_pairs)
    data = data[(data[EDIT_INDEX] != 0)
                & (~data[BARCODE_NAME].str.contains("newOverlap", na=False))]

    # Collapse each barcode to a single effect across replicates, then split barcodes.
    per_barcode = (data.groupby([EDIT_INDEX, BARCODE_NAME])
                   .agg(**{PEG_EFFECT: (PEG_EFFECT, "mean"),
                           "_n_biorep1": ("BioRep", lambda s: (s == 1).sum()),
                           "_n_biorep2": ("BioRep", lambda s: (s == 2).sum())})
                   .reset_index())
    before = len(per_barcode)
    per_barcode = per_barcode[(per_barcode["_n_biorep1"] > 0)
                              & (per_barcode["_n_biorep2"] > 0)]
    print(f"  barcode consistency: {len(per_barcode)} barcodes "
          f"({before - len(per_barcode)} dropped for missing a biological replicate)")
    rng = np.random.default_rng(SPLIT_SEED)

    for minimum in thresholds:
        halves = []
        for edit, group in per_barcode.groupby(EDIT_INDEX):
            if len(group) < minimum:
                continue
            order = rng.permutation(len(group))
            split = len(group) // 2
            values = group[PEG_EFFECT].values
            halves.append((values[order[:split]].mean(), values[order[split:]].mean()))
        if len(halves) < 3:
            print(f"  skipping barcode consistency min{minimum}: only {len(halves)} edits")
            continue
        x, y = np.array(halves).T
        r = stats.pearsonr(x, y)[0]

        fig, ax = plt.subplots(figsize=(5.6, 5.6))
        limits = scatter_limits(np.concatenate([x, y]))
        ax.plot(limits, limits, "k--", linewidth=1.2, zorder=1)
        ax.scatter(x, y, s=45, color="#4B85C5", edgecolor="black", linewidth=0.4, zorder=2)
        ax.text(0.05, 0.93, f"Pearson's r = {r:.2f}\nn = {len(halves)}",
                transform=ax.transAxes, fontsize=12, va="top")
        ax.set_xlim(*limits)
        ax.set_ylim(*limits)
        ax.set_xlabel("Barcode group 1 mean effect (%)", fontsize=13)
        ax.set_ylabel("Barcode group 2 mean effect (%)", fontsize=13)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"Split-half barcode consistency, "
                     f"edits with $\\geq${minimum} barcodes", fontsize=10)
        fig.tight_layout()
        out = os.path.join(outdir, f"barcode_consistency_min{minimum}.pdf")
        fig.savefig(out)
        plt.close(fig)
        print(f"  wrote {out} (n={len(halves)}, r={r:.3f})")


def plot_pool_allele_frequency(unfiltered, outdir):
    """The two pool-editing panels. Reads the table BEFORE the min-allele-frequency cut."""
    edits = unfiltered[unfiltered[EDIT_INDEX] != 0].copy()

    ## panel 1: edited fraction per BioRep
    totals = edits.groupby(["BioRep", "FFRep"])[SPIKEIN_FREQ].sum().reset_index()
    totals["edited_pct"] = totals[SPIKEIN_FREQ] * 100
    bioreps = sorted(totals["BioRep"].unique())
    means, cis, dots = [], [], []
    for biorep in bioreps:
        values = totals.loc[totals["BioRep"] == biorep, "edited_pct"].values
        mean = values.mean()
        half_width = (stats.t.ppf(0.975, df=len(values) - 1)
                      * values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0
        means.append(mean)
        cis.append(half_width)
        dots.append(values)
        print(f"  BioRep {biorep}: n={len(values)}, mean={mean:.2f}%, 95% CI=+/-{half_width:.2f}%")

    fig, ax = plt.subplots(figsize=(3.0, 4))
    x = np.arange(len(bioreps))
    ax.bar(x, means, yerr=cis, color=BAR_COLOR, edgecolor="none", width=0.35, capsize=4,
           error_kw={"elinewidth": 1.2, "ecolor": "black"})
    rng = np.random.default_rng(0)
    for position, values in zip(x, dots):
        jitter = rng.uniform(-0.05, 0.05, size=len(values))
        ax.scatter(position + jitter, values, color="black", s=18, alpha=0.6,
                   edgecolors="none", zorder=3)
    for position, mean, half_width in zip(x, means, cis):
        ax.text(position, mean + half_width + 0.3, f"{mean:.2f} +/- {half_width:.2f}%",
                ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels([f"BioRep {b}" for b in bioreps])
    ax.set_ylabel("Edited reads (%)", fontsize=12)
    ax.set_ylim(0, 10)
    ax.set_title("Pool editing \u2014 allele frequency", fontsize=11)
    fig.tight_layout()
    out = os.path.join(outdir, "pool_allele_frequency_by_biorep.pdf")
    fig.savefig(out, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  wrote {out}")

    ## panel 2: per-barcode allele-frequency distribution
    per_barcode = edits.groupby(BARCODE_NAME)[SPIKEIN_FREQ].mean().reset_index()
    per_barcode["af_pct"] = per_barcode[SPIKEIN_FREQ] * 100
    threshold_pct = ALLELE_FREQ_THRESHOLD * 100
    n_total = len(per_barcode)
    n_pass = int((per_barcode[SPIKEIN_FREQ] >= ALLELE_FREQ_THRESHOLD).sum())

    ## Summing A-D gives cells sorted; averaged over BioReps.
    bins_present = [c for c in SORT_BIN_COLUMNS if c in unfiltered.columns]
    cells_sorted = (unfiltered.groupby("BioRep")[bins_present].sum().sum(axis=1).mean()
                    if bins_present else float("nan"))

    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    positive = per_barcode.loc[per_barcode["af_pct"] > 0, "af_pct"]
    edges = np.logspace(np.log10(positive.min()), np.log10(per_barcode["af_pct"].max()), 60)
    ax.hist(per_barcode["af_pct"], bins=edges, color=BAR_COLOR, edgecolor="none", alpha=0.9)
    ax.set_xscale("log")
    label = f"Threshold ({threshold_pct:g}%"
    if np.isfinite(cells_sorted):
        label += f" \u2248 {ALLELE_FREQ_THRESHOLD * cells_sorted:,.0f} cells"
    ax.axvline(threshold_pct, color="red", linestyle="--", linewidth=1.5, label=label + ")")
    ax.text(0.97, 0.95,
            f"Total barcodes: {n_total:,}\nPassing threshold: {n_pass:,} "
            f"({100 * n_pass / n_total:.1f}%)",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#cccccc",
                      alpha=0.9))
    ax.set_xlabel("Mean allele frequency (%)", fontsize=12)
    ax.set_ylabel("Number of barcodes", fontsize=12)
    ax.legend(fontsize=9, loc="upper left", frameon=False)
    if np.isfinite(cells_sorted):
        secondary = ax.secondary_xaxis(
            "top", functions=(lambda v: v / 100 * cells_sorted,
                              lambda v: v / cells_sorted * 100))
        secondary.set_xlabel(f"Cells sorted per BioRep "
                             f"(mean of BioReps, {cells_sorted / 1e6:.0f}M total)",
                             fontsize=11)
    fig.tight_layout()
    out = os.path.join(outdir, "pool_allele_frequency_distribution_combined.pdf")
    fig.savefig(out, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    print(f"  wrote {out} ({n_pass}/{n_total} barcodes pass {threshold_pct:g}%; "
          f"{cells_sorted:,.0f} cells sorted per BioRep)")


def main():
    args = parse_args()
    style()
    os.makedirs(args.outdir, exist_ok=True)

    reps = pd.read_table(args.replicate_effects)
    pegrnas = pd.read_table(args.pegrna_effects)

    excluded_pairs = excluded_replicate_pairs(reps)

    ## Both levels: per edit and per pegRNA.
    levels = [
        ("perEdit", reps, EDIT_INDEX, REP_EFFECT, "one dot per edit"),
        ("perPegRNA", pegrnas[pegrnas[EDIT_INDEX] != 0], BARCODE_NAME, PEG_EFFECT,
         "one dot per pegRNA"),
    ]
    for suffix, table, key, effect, level in levels:
        plot_biorep_scatter(
            table, key, effect, excluded_pairs,
            os.path.join(args.outdir, f"replicate_correlation_biorep_{suffix}.pdf"), level)
        plot_ffrep_scatter_pages(
            table, key, effect, excluded_pairs,
            os.path.join(args.outdir, f"replicate_correlation_ffreps_{suffix}.pdf"), level)

    plot_replicate_correlation(
        reps, os.path.join(args.outdir, "replicate_correlation_all.pdf"),
        "All replicates")
    retained = reps[~reps["Excluded"].astype(str).str.upper().isin({"TRUE", "1"})]
    if len(retained) < len(reps):
        plot_replicate_correlation(
            retained, os.path.join(args.outdir, "replicate_correlation_filtered.pdf"),
            "Retained replicates (exclude_replicates applied)")
    else:
        print("  no replicates excluded; skipping the filtered heatmap")

    plot_barcode_consistency(pegrnas, excluded_pairs, args.outdir)

    if args.pegrna_effects_unfiltered:
        plot_pool_allele_frequency(pd.read_table(args.pegrna_effects_unfiltered),
                                   args.outdir)
    else:
        print("  no --pegrna-effects-unfiltered given; skipping the pool-editing panels")


if __name__ == "__main__":
    main()
