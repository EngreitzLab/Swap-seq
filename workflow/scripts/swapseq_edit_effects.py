#! /usr/bin/env python

## Swap-seq per-edit effect sizes.
##
## Takes the MLE allelic effects table and produces one row per edit, reading
## `effect_size_scaled_qpcr` (the effect before any heterozygosity adjustment):
##   1. Anchor allele frequencies on the barcode-only spike-ins, which are a known
##      fraction of all alleles
##   2. Heterozygous correction  effect = 1 + 2 * (effect_size_scaled_qpcr - 1)/(1 + g).
##   3. Aggregate barcodes -> edits, filter, one-sample t-test, Benjamini-Hochberg.

import argparse
import re

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, ttest_1samp

CONTROL_PREFIX = "NO_PegRNA_Negative_Control_Barcode-"
REFERENCE_ID = "Swap-seq_Barcode:InferredReference"
SPIKEIN_COLUMN = "SpikeInControl"
CATEGORY_COLUMN = "DesignCategory"
CONTROL_CATEGORY = "Negative controls"

PEGRNA_DROP_COLUMNS = ["ExperimentIDReplicates", "guide", "Location", "Variant"]

## Output headers
PEGRNA_HEADERS = {
    "PegRNAIndex": "Swap-seq pegRNA index (per barcode)",
    "EditIndex": "Swap-seq edit index (per edit)",
    "VariantID": "Deletion Pair Name (per barcode)",
    "EditName": "Swap-seq edit name (per edit)",
    "DesignCategory": "Design categories",
    "anchored_freq": "Swap-seq freq",
    "editing_efficiency": "Swap-seq freq guide adjusted",
    "effect_ratio": "Swap-seq final effect size",
    "effect_percent": "Swap-seq final effect size percent",
}
EDIT_HEADERS = {
    "EditIndex": "Swap-seq edit index (per edit)",
    "EditName": "Swap-seq edit name (per edit)",
    "DesignCategory": "Design categories",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("effects_table",
                        help="V1 allelic effects table from the pipeline (.tsv or .tsv.gz)")
    parser.add_argument("output_file", help="per-edit effect table to write (.tsv)")
    parser.add_argument("--pegrna-output", default="",
                        help="if given, also write the per-pegRNA (per-barcode x "
                             "replicate) effects here. These are the same Swap-seq "
                             "adjusted values the per-edit table aggregates, before "
                             "barcodes are averaged into edits.")
    parser.add_argument("--unfiltered-pegrna-output", default="",
                        help="if given, also write the per-pegRNA table BEFORE the "
                             "--min-allele-freq cut. The QC pool-editing panels need it: "
                             "the point of those plots is how many designed barcodes clear "
                             "the frequency threshold, which the filtered table cannot show.")
    parser.add_argument("--replicate-output", default="",
                        help="if given, also write per-edit effects broken out by "
                             "replicate (edit x BioRep x FFRep). The level the QC and "
                             "per-locus plots read: barcodes averaged within an edit, "
                             "replicates not yet collapsed.")
    parser.add_argument("--variant-info", default="",
                        help="Variant Table. If it has a %(const)s column, the alleles "
                             "marked TRUE there are treated as barcode-only spike-in "
                             "controls: they anchor the frequency scale and receive no "
                             "heterozygous correction. Without this the script falls back "
                             "to matching VariantIDs against a name prefix."
                             % {"const": SPIKEIN_COLUMN})
    parser.add_argument("--spikein-allele-fraction", type=float, default=0.075,
                        help="fraction of ALL alleles carrying a barcode-only control "
                             "cassette. Controls are spiked into the plasmid library at a "
                             "known rate and land on one of two alleles, so this is "
                             "spike_rate x copies_per_cell x 0.5 (default: 0.075, i.e. "
                             "15%% spike-in and single-copy infection)")
    parser.add_argument("--min-allele-freq", type=float, default=1e-5,
                        help="drop barcode x replicate rows below this anchored frequency")
    parser.add_argument("--min-barcodes-test", type=int, default=2,
                        help="barcodes required to test an edit")
    parser.add_argument("--min-barcodes-report", type=int, default=3,
                        help="barcodes required to report a non-significant edit")
    parser.add_argument("--bh-threshold", type=float, default=0.001,
                        help="Benjamini-Hochberg threshold for calling significance")
    parser.add_argument("--exclude-replicates", default="",
                        help="comma-separated BioRep:FFRep pairs to drop everywhere, "
                             "e.g. '1:4,1:5' for sorts that correlate poorly with the rest")
    return parser.parse_args()


def spikein_control_ids(variant_info):
    """VariantIDs of the barcode-only spike-ins. None if the column is absent."""
    if not variant_info:
        return None
    variants = pd.read_table(variant_info)
    if SPIKEIN_COLUMN not in variants.columns:
        return None
    flag = variants[SPIKEIN_COLUMN].astype(str).str.strip().str.upper()
    return set(variants.loc[flag.isin({"TRUE", "T", "1", "YES"}), "VariantID"])


def design_categories(variant_info):
    """Edit index -> design category. Index 0 = controls."""
    if not variant_info:
        return {}
    variants = pd.read_table(variant_info)
    if CATEGORY_COLUMN not in variants.columns:
        return {}
    categories = {}
    for variant_id, category in zip(variants["VariantID"],
                                    variants[CATEGORY_COLUMN].astype(str).str.strip()):
        categories.setdefault(edit_index(variant_id), category)
    return categories


def pegrna_index(variant_id):
    """The last integer of a VariantID = the pegRNA index. 0 if absent."""
    match = re.search(r"_(\d+)$", variant_id)
    return int(match.group(1)) if match else 0


def parse_excluded(spec):
    """for example, if 1:4 and 1:5 are excluded, turn '1:4,1:5' into {(1, 4), (1, 5)}."""
    pairs = set()
    for item in filter(None, (s.strip() for s in spec.split(","))):
        biorep, ffrep = item.split(":")
        pairs.add((int(biorep), int(ffrep)))
    return pairs


def edit_index(variant_id):
    """The integer before `_Barcode` = the edit index. 0 for controls and the reference."""
    match = re.search(r"_(\d+)_Barcode", variant_id)
    return int(match.group(1)) if match else 0


def edit_name(variant_id):
    """Collapse a per-barcode VariantID to its edit name by stripping the suffix."""
    if edit_index(variant_id) == 0:
        return variant_id
    if "newOverlap" in variant_id:
        return re.sub(r"_\d+$", "", re.sub(r"_Barcode-\d+", "", variant_id))
    return re.sub(r"_Barcode.*$", "", variant_id)


def anchor_and_correct(alleles, spikein_fraction):
    """ The controls are a known fraction of all alleles, so scaling by
    spikein_fraction / sum(control freqs) puts every frequency on that scale. Dividing by
    the barcode's plasmid frequency then gives the editing efficiency.
    """
    alleles = alleles.copy()

    # The reference is itself a spiked control cassette, so it belongs in the denominator.
    control_freq = (alleles[alleles["IsSpikeIn"]]
                    .groupby(["BioRep", "FFRep"])["freq"].sum())
    missing = set(zip(alleles["BioRep"], alleles["FFRep"])) - set(control_freq.index)
    if missing:
        raise SystemExit(f"no barcode-only control alleles in replicate(s) {sorted(missing)}; "
                         f"cannot anchor the frequency scale")

    scale = np.array([spikein_fraction / control_freq[(b, f)]
                      for b, f in zip(alleles["BioRep"], alleles["FFRep"])])
    alleles["anchored_freq"] = alleles["freq"].values * scale
    alleles["editing_efficiency"] = np.minimum(
        alleles["anchored_freq"] / alleles["guide_freq"], 1.0)
    alleles.loc[alleles["IsSpikeIn"], "editing_efficiency"] = 1.0

    # 2/(1 + g) undoes the dilution from measuring a heterozygous edit.
    alleles["effect_ratio"] = 1 + (alleles["effect_size_scaled_qpcr"] - 1) * 2 / (
        1 + alleles["editing_efficiency"])
    alleles["effect_percent"] = (alleles["effect_ratio"] - 1) * 100
    return alleles


def barcode_support_filter(per_rep, minimum):
    """Keep rows with >= `minimum` barcodes. """
    exempt = (per_rep["IsSpikeIn"]
              | per_rep["EditName"].str.contains("newOverlap", na=False))
    return per_rep[exempt | (per_rep["n_barcodes"] >= minimum)]


def has_multiple_bioreps(group):
    return group["BioRep"].nunique() >= 2


def test_effects(per_rep, excluded, min_barcodes):
    """One-sample t-test and BH correction."""
    supported = barcode_support_filter(per_rep, min_barcodes)
    is_control = supported["IsSpikeIn"]

    def drop_excluded(df):
        mask = [(b, f) in excluded for b, f in zip(df["BioRep"], df["FFRep"])]
        return df[~np.array(mask, dtype=bool)] if len(df) else df

    # Average barcodes within a replicate first.
    edit_reps = (supported[~is_control]
                 .groupby(["EditIndex", "BioRep", "FFRep"])["effect_percent"]
                 .mean().reset_index())
    edit_reps = drop_excluded(edit_reps)
    edit_reps = edit_reps[edit_reps["EditIndex"] != 0]
    edit_reps = edit_reps.groupby("EditIndex").filter(has_multiple_bioreps)

    control_reps = drop_excluded(supported[is_control]).groupby("EditName").filter(
        has_multiple_bioreps)
    name_of = (supported.groupby("EditIndex")["EditName"]
               .agg(lambda names: min(names, key=len)).to_dict())

    def summarize(group):
        """Test one edit's replicate effects against no effect."""
        values = group["effect_percent"].values
        sd = values.std(ddof=1)
        return {"mean_effect_percent": values.mean(),
                "n_replicates": len(values),
                "n_biorep1": int((group["BioRep"] == 1).sum()),
                "n_biorep2": int((group["BioRep"] == 2).sum()),
                "p_value": ttest_1samp(values, 0)[1],
                "SD": sd, "Variance": sd ** 2}

    rows = []
    for key, group in edit_reps.groupby("EditIndex"):
        rows.append({"EditIndex": key, "EditName": name_of.get(key, ""), "kind": "edit",
                     **summarize(group)})
    for key, group in control_reps.groupby("EditName"):
        rows.append({"EditIndex": 0, "EditName": key, "kind": "control",
                     **summarize(group)})

    results = pd.DataFrame(rows)
    results["BH_p_value"] = false_discovery_control(results["p_value"])
    return results


def apply_display_filter(results, per_rep, excluded, min_barcodes, bh_threshold):
    """Drop non-significant edits with too few barcodes."""
    well_supported = barcode_support_filter(per_rep, min_barcodes)
    well_supported = well_supported[~well_supported["IsSpikeIn"]]
    mask = [(b, f) in excluded for b, f in zip(well_supported["BioRep"],
                                               well_supported["FFRep"])]
    well_supported = well_supported[~np.array(mask, dtype=bool)]
    keep = set(well_supported.groupby("EditIndex").filter(has_multiple_bioreps)["EditIndex"])

    significant = set(results.loc[(results["kind"] == "edit")
                                 & (results["BH_p_value"] < bh_threshold), "EditIndex"])
    keep |= significant
    return results[(results["kind"] == "control") | results["EditIndex"].isin(keep)]


def write_replicate_table(per_rep, categories, excluded, outfile):
    """One row per edit x (BioRep, FFRep)"""
    table = per_rep[~per_rep["IsSpikeIn"]].copy()
    table = table[table["EditIndex"] != 0]

    # Collapse to the parent name (the shortest), as the per-edit test does.
    grouped = (table.groupby(["EditIndex", "BioRep", "FFRep"])
               .agg(mean_effect_percent=("effect_percent", "mean"),
                    n_barcodes=("n_barcodes", "sum"),
                    n_edit_names=("EditName", "nunique"),
                    EditName=("EditName", lambda names: min(names, key=len)))
               .reset_index())
    grouped["DesignCategory"] = grouped["EditIndex"].map(categories).fillna("")
    grouped["Excluded"] = [(b, f) in excluded
                           for b, f in zip(grouped["BioRep"], grouped["FFRep"])]

    columns = ["EditIndex", "EditName", "DesignCategory", "BioRep", "FFRep",
               "mean_effect_percent", "n_barcodes", "n_edit_names", "Excluded"]
    (grouped[columns]
     .rename(columns=EDIT_HEADERS)
     .sort_values(list(EDIT_HEADERS.get(c, c) for c in ["EditIndex", "BioRep", "FFRep"]))
     .to_csv(outfile, sep="\t", index=False))
    print(f"wrote {outfile}")


def write_pegrna_table(filtered, outfile):
    """One row per barcode x replicate."""
    filtered = filtered.copy()

    # The column is "(per edit)", so it holds the edit's name, not the barcode's; these
    # differ for `newOverlap`. Index 0 (controls, reference) keeps its own name.
    parent = (filtered[filtered["EditIndex"] != 0]
              .groupby("EditIndex")["EditName"]
              .agg(lambda names: min(names, key=len)))
    is_edit = filtered["EditIndex"] != 0
    filtered.loc[is_edit, "EditName"] = filtered.loc[is_edit, "EditIndex"].map(parent)

    lead = ["PegRNAIndex", "EditIndex", "VariantID", "EditName", "DesignCategory"]
    tail = ["anchored_freq", "editing_efficiency", "effect_ratio", "effect_percent"]
    middle = [c for c in filtered.columns
              if c not in lead + tail + PEGRNA_DROP_COLUMNS
              # internal bookkeeping, not part of the published table
              and c not in ("IsSpikeIn",)]
    table = filtered[lead + middle + tail].rename(columns=PEGRNA_HEADERS)
    table.sort_values(list(table.columns[:3])).to_csv(outfile, sep="\t", index=False)
    print(f"wrote {outfile}")


def main():
    args = parse_args()
    excluded = parse_excluded(args.exclude_replicates)

    alleles = pd.read_table(args.effects_table)
    alleles["EditIndex"] = alleles["VariantID"].map(edit_index)
    alleles["EditName"] = alleles["VariantID"].map(edit_name)
    alleles["PegRNAIndex"] = alleles["VariantID"].map(pegrna_index)

    categories = design_categories(args.variant_info)
    if categories:
        print(f"  design categories: {len(set(categories.values()))} distinct, "
              f"over {len(categories)} edits")
    else:
        print(f"  design categories: none (no {CATEGORY_COLUMN} column supplied)")
    alleles["DesignCategory"] = alleles["EditIndex"].map(categories).fillna("")
    # The reference is not in the Variant Table but is a control cassette.
    if categories:
        alleles.loc[alleles["VariantID"] == REFERENCE_ID, "DesignCategory"] = (
            categories.get(0, CONTROL_CATEGORY))

    # Declared in the Variant Table when possible, else matched on the name prefix.
    declared = spikein_control_ids(args.variant_info)
    if declared is None:
        print(f"  spike-in controls: from name prefix {CONTROL_PREFIX!r} "
              f"(no {SPIKEIN_COLUMN} column supplied)")
        alleles["IsSpikeIn"] = alleles["VariantID"].str.startswith(CONTROL_PREFIX)
    else:
        print(f"  spike-in controls: {len(declared)} declared in {args.variant_info}")
        alleles["IsSpikeIn"] = alleles["VariantID"].isin(declared)

    # check before adding the reference, which would otherwise hide an empty set
    if not alleles["IsSpikeIn"].any():
        raise SystemExit("no spike-in control alleles found; cannot anchor the frequency "
                         f"scale. Mark them TRUE in a {SPIKEIN_COLUMN} column of the "
                         "Variant Table.")
    alleles.loc[alleles["VariantID"] == REFERENCE_ID, "IsSpikeIn"] = True

    alleles = anchor_and_correct(alleles, args.spikein_allele_fraction)

    filtered = alleles[alleles["anchored_freq"] >= args.min_allele_freq]
    kept = filtered[filtered["VariantID"] != REFERENCE_ID]
    per_rep = (kept.groupby(["EditName", "EditIndex", "IsSpikeIn", "BioRep", "FFRep"])
               .agg(effect_percent=("effect_percent", "mean"),
                    n_barcodes=("VariantID", "count"))
               .reset_index())

    results = test_effects(per_rep, excluded, args.min_barcodes_test)
    reported = apply_display_filter(results, per_rep, excluded,
                                    args.min_barcodes_report, args.bh_threshold)

    counts = {
        "rows": len(filtered),
        "barcodes": filtered["VariantID"].nunique(),
        "tested_edits": int((results["kind"] == "edit").sum()),
        "tested_controls": int((results["kind"] == "control").sum()),
        "reported_edits": int((reported["kind"] == "edit").sum()),
        "significant": int(((reported["kind"] == "edit")
                            & (reported["BH_p_value"] < args.bh_threshold)).sum()),
    }
    for key, value in counts.items():
        print(f"  {key:<18} {value}")

    if args.unfiltered_pegrna_output:
        write_pegrna_table(alleles, args.unfiltered_pegrna_output)
    if args.pegrna_output:
        write_pegrna_table(filtered, args.pegrna_output)
    if args.replicate_output:
        write_replicate_table(per_rep, categories, excluded, args.replicate_output)

    # One edit has one category.
    by_index = (kept[kept["EditIndex"] != 0].groupby("EditIndex")["DesignCategory"]
                .first().to_dict())
    by_name = kept.groupby("EditName")["DesignCategory"].first().to_dict()
    reported = reported.copy()
    reported["DesignCategory"] = [
        by_name.get(name, "") if index == 0 else by_index.get(index, "")
        for index, name in zip(reported["EditIndex"], reported["EditName"])]
    significant = f"Significant (BH<{args.bh_threshold:g})"
    reported[significant] = reported["BH_p_value"] < args.bh_threshold

    columns = ["EditIndex", "EditName", "mean_effect_percent", "n_replicates",
               "n_biorep1", "n_biorep2", "p_value", "SD", "Variance", "BH_p_value",
               significant, "DesignCategory"]
    (reported.sort_values(["kind", "EditIndex", "EditName"])[columns]
     .rename(columns=EDIT_HEADERS)
     .to_csv(args.output_file, sep="\t", index=False))
    print(f"wrote {args.output_file}")


if __name__ == "__main__":
    main()
