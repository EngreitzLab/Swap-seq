# Swap-seq analysis pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.0-brightgreen.svg)](https://snakemake.bitbucket.io)

Swap-seq measures the effect of many programmed deletions on a gene's expression at once.
Pooled twin prime editing removes target sequences and installs a short barcode cassette in its place.
Cells are sorted on the target gene's RNA-FISH signal and gDNA is extracted from each sorted bin.
Barcode cassettes are amplified and sequenced. Editing rates are
quantified with [CRISPResso2](https://github.com/pinellolab/CRISPResso2), effect sizes are
inferred from each barcode's distribution across FACS bins by maximum likelihood, and
significance is computed with 1-sample t-tests corrected for multiple testing.

The workflow is a fork of
[Variant-EFFECTS](https://github.com/EngreitzLab/Variant-EFFECTS). Read alignment, barcode counting and the maximum-likelihood effect-size
estimator are unchanged; what differs is described under
[Differences from Variant-EFFECTS](#differences-from-variant-effects).

## Contents

- [Usage](#usage)
  - [Input 1: Sample List](#input-1-sample-list)
  - [Input 2: Amplicon Table](#input-2-amplicon-table)
  - [Input 3: Sorting parameters files](#input-3-sorting-parameters-files)
  - [Input 4: Variant Table](#input-4-variant-table)
  - [Input 5: Guide Counts Table](#input-5-guide-counts-table)
- [Workflow](#workflow)
  - [Step 0: Demultiplex](#step-0-demultiplex)
  - [Step 1: Clone this github repository](#step-1-clone-this-github-repository)
  - [Step 2: Install conda environment](#step-2-install-conda-environment)
  - [Step 3: Trim reads to the cassette](#step-3-trim-reads-to-the-cassette)
  - [Step 4: Configure workflow](#step-4-configure-workflow)
  - [Step 5: Execute workflow](#step-5-execute-workflow)
- [Differences from Variant-EFFECTS](#differences-from-variant-effects)
- [Outputs](#outputs)
  - [Structure](#structure)
  - [Tables in `summary/`](#tables-in-summary)
  - [Key files](#key-files)
  - [Plots](#plots)
- [License](#license)

## Usage

There are several required inputs prior to executing this pipeline. For ease
of use, generate a subdirectory `config/` in the directory you are performing the data
analysis. Generate and place the following documents inside it.

A complete, working set of all five inputs used for the *PPIF* screen is in
[`example/`](example/), together with [`example/log.sh`](example/log.sh) recording every
command that was run, and the demultiplexing Sample Sheet from
[Step 0](#step-0-demultiplex).

### Input 1: Sample List

The Sample List lists all of the sequencing libraries that will be included in the
analysis, and describes their relationships and groupings. See
[`example/config/SampleList.tsv`](example/config/SampleList.tsv).

Required columns:

    SampleID          Unique name for each amplicon library. (e.g., BATCH-CellLine-Sample-FFRep-PCRRep-Bin)
    AmpliconID        Name of genomic amplicon contained in the library - must match corresponding AmpliconID column in the Amplicon Table (see below)
                        For Swap-seq this is a single value for the whole screen ('Swap-seq_Barcode'), because
                        every deletion is read out through the same barcode cassette amplicon.

    Batch             Batch ID used to identify the appropriate FACS sort params file (config['sortparamsdir']/{Batch}_{SampleNumber}.csv)
    SampleNumber      FlowFISH sample number - used to identify the appropriate FACS sort params file (config['sortparamsdir']/{Batch}_{SampleNumber}.csv)

    Bin               Name of a FACS-sorted bin (e.g.: A B C D). 'All' for FlowFISH-input edited samples. 'Neg' or blank if not applicable
    PCRRep            PCR replicate number or name
    ControlForAmplicon TRUE or FALSE. Set to TRUE for unedited samples that will be used to evaluate background sequencing/PCR error rate

    [Experiment Keys] Provide any number of additional columns (e.g., CellLine, Guides, TestProbe) that distinguish different samples.
                        Key columns are defined as such by the 'experiment_keycols' parameter in the config file.
    [Replicate Keys]  Provide any number of additional columns (e.g., FlowFISHRep) that distinguish different experimental replicates.
                        Replicate columns are defined by the 'replicate_keycols' parameter in the config file.
                        PCR replicate counts for each unique replicate key will be summed at the level of this replicate ID.

Note: `Bin` must NOT be listed in `replicate_keycols`, and amplicon names must not contain
spaces.

The *PPIF* screen used 2 BioReps x 7 FlowFISH reps x (4 sorted bins + 1 unsorted input) x
16 PCR reps = **924** libraries: `Bin` is `A`/`B`/`C`/`D` for the sorted bins and `All`
for the unsorted input.

### Input 2: Amplicon Table

The Amplicon Table lists details for the genomic PCR amplicons used in the experiment.
Information is pulled into the Sample List by the `AmpliconID` column. See
[`example/config/AmpliconInfo.tsv`](example/config/AmpliconInfo.tsv).

Required columns:

    AmpliconID                  User-defined name of the amplicon.
    AmpliconSeq                 Full sequence to align to.
    QuantificationWindowStart   Input for CRISPResso. Zero-based coordinate for quantifying reference allele.
    QuantificationWindowEnd     Input for CRISPResso. Zero-based coordinate for quantifying reference allele.
    ReferenceErrorThreshold     Integer indicating how many errors (mismatch/insertion/deletion) are tolerable when inferring the reference allele.

### Input 3: Sorting parameters files

This file lists statistics and values derived from the FACS sort for each sample. Files
must be named `{Batch}_{SampleNumber}.csv` and placed in the `config['sortparamsdir']`
directory. See [`example/sortParams/`](example/sortParams/) — one file per replicate.

Required columns:

    Name               Sorting gate name on the cytometer
    Barcode            Name of the sorted bin, needs to match "Bin" column in the Sample List
    Count              Number of cells sorted into this bin
    Mean               Mean fluorescence values of cells sorted into this bin
    Min                Minimum fluorescence value sorted into this bin (e.g., left edge of the gate)
    Max                Maximum fluorescence value sorted into this bin (e.g., right edge of the gate)

These tables were made from the sorting result output files of the FACS machine.

The `Barcode` column must also contain a row named `Total` for the whole population the
bins were drawn from — the pipeline stops with an error without it.

### Input 4: Variant Table

The Variant Table lists details for all the Swap-seq edits in the experiment. See
[`example/config/VariantList.tsv`](example/config/VariantList.tsv).

Required columns:

    AmpliconID         User-defined name of the amplicon that matches AmpliconID in the provided amplicon table
    VariantID          Unique readable name of the edit
    MappingSequence    Barcode cassette sequence corresponding to the specific Swap-seq edit; used for alignment
    RefAllele          TRUE/FALSE if this is (one of) the reference alleles. Always FALSE for
                         Swap-seq: the reference is the one barcode-only cassette held out of this
                         table, which the pipeline itself emits as "Swap-seq_Barcode:InferredReference"
                         and flags as the reference. The column must be present either way
    SpikeInControl     TRUE/FALSE. TRUE marks a barcode-only control cassette: no pegRNA, so no edit
                         is made. These were used as negative controls to ensure the effect sizes are calibrated and
                         for anchoring allele frequency calculating. See Method --> "Analysis of Swap-seq screens section" of the manuscript
    DesignCategory     What kind of element the edit targets: PPIF promoter, PPIF splice junction,
                         CRISPRi hits, Candidate enhancer, CTCF ChIP-seq peaks, Other candidate elements,
                         Nearby-gene promoter, Negative controls. Carried through to both output tables.
                         Optional -- the category columns come out blank if it is absent

For Swap-seq one row is one **barcoded construct**, and `MappingSequence` is the **full
71 bp cassette**. As an example, the current `VariantID` encodes the edit:

```
chr10:81080709-81081391_Nasser_enhancer_Peak96_77_Barcode-1_305

Element:                                 "chr10:81080709-81081391_Nasser_enhancer_Peak96"
Index of the edit:                       "77"
Index of barcode for this edit (1 to 4): "Barcode-1"
Index of the pegRNA:                     "305"
```

Note: `MappingSequence` must be UPPER CASE to match CRISPResso2's output. 

### Input 5: Guide Counts Table

Counts of each pegRNA in the plasmid library, used to normalize allele frequency ("what's the proportion of all sequenced alleles that align to edit X") 
by how abundant that pegRNA was in the library, which gives us the per pegRNA editing efficiency ("what's the portion of all sequenced alleles from cells that received pegRNA X that align to edit X"). 
Editing efficiencies are used for heterozygous adjustment during the effect size calculation.

```
43936   *
25697   chr10:80828483-80828842_ZMIZ1_promoter_1_Barcode-1_1
25245   chr10:80828483-80828842_ZMIZ1_promoter_1_Barcode-2_2
...
```

See
[`example/config/SC_Swap-seq_TDL_plasmid_library_PCRRep1_S1_read1.count.txt`](example/config/SC_Swap-seq_TDL_plasmid_library_PCRRep1_S1_read1.count.txt).


## Workflow

### Step 0: Demultiplex

Demultiplexing is not part of the Snakemake workflow. See
[`example/demultiplex.batch`](example/demultiplex.batch) for the `bcl2fastq` call used:

    bcl2fastq --runfolder-dir $BCLDIR --output-dir $PROJECT/fastq/ \
      --sample-sheet $PROJECT/config/SampleSheet.csv \
      --use-bases-mask Y*,I8N2,I8N2,Y* --no-lane-splitting \
      --create-fastq-for-index-reads --barcode-mismatches 0


The **Sample Sheet** [`example/config/SampleSheet.csv`](example/config/SampleSheet.csv)
is the bcl2fastq demultiplexing sheet, and is a different file from the **Sample List** in [Input 1](#input-1-sample-list)

### Step 1: Clone this github repository

[Clone](https://help.github.com/en/articles/cloning-a-repository) this to your local system
or server where you want to perform the data analysis.

### Step 2: Install conda environment

Install the "VFFenv" conda environment using
[conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or
[mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).

    mamba env create -f envs/VFFenv_release_240702.yml
    #or
    conda env create --file envs/VFFenv_dev_240702.yml

### Step 3: Trim reads to the cassette

    bash Swap-seq/example/trim_fastq.sh

It reads `fastq/*.fastq.gz` and writes `trimed_fastq/` — so run it from your analysis directory,
the one holding `config/`. Later in the config file, set the path to fastq files as `trimed_fastq`:
[`example/config/config.json`](example/config/config.json)

### Step 4: Configure workflow

Copy [`example/config/config.json`](example/config/config.json) to your `config/` folder
and edit the fields to point to the right files.

The fields, in the order they appear in that file:

    sample_list                          Path to the Sample List
    amplicon_info                        Path to the Amplicon Table
    variant_info                         Path to the Variant Table
    sortparamsdir                        Path to the directory of FACS sort parameters files
    fastqdir                             Path to the trimmed FASTQ directory (see Step 3), not the raw one
    codedir                              Path to this code (e.g. "Swap-seq/")
    replicate_keycols                    Comma-separated Sample List columns marking experimental replicates
    experiment_keycols                   Comma-separated Sample List columns marking different experiments
    crispresso_min_average_read_quality  Passed to CRISPResso2 as -q; reads averaging below this are discarded
    crispresso_min_single_bp_quality     Passed to CRISPResso2 as -s; 0 disables per-base filtering
    genotyping_only                      "True" to stop after quantifying variants, skipping effect sizes
    reps                                 Number of PCR replicates
    pooled                               "True" for a pool of variants, "False" for a single variant
    ff_tss_guide_kd                      Scaling factor between FlowFISH and qPCR measurements of a total
                                           knockdown of the gene of interest, to correct for non-specific
                                           FlowFISH probe binding. Set to 1 if unknown
    qpcr_tss_guide_kd                    See above. Set to 1 if unknown
    guide_counts_file                    Path to the Guide Counts Table
    spikein_plasmid_fraction             Fraction of the plasmid library that is barcode-only spike-in
                                           cassettes. This value is used for adjusting allele frequencies (see Method: Analysis of Swap-seq screens of the manuscript)
                                           Under low MOI lentiviral infection, we assume the alleles carrying the control barcodes will be half of this ratio
                                           as lentivirus integrates in a hemizygous way. Default is 0.15 (i.e. 15% of the library)
    exclude_replicates                   Set BioRep:FFRep pairs to drop for further analysis if, for example, there are known technical issues with some replicates.
                                         BioRep 1 FlowFISH reps 4 and 5 were dropped in the screen due to known experimental error, which we set as "1:4,1:5" here.

Optional and not in the example:

    single_end                           "True" for single-end rather than paired-end reads. Defaults to
                                           False when omitted

### Step 5: Execute workflow

Activate the conda environment:

    conda activate VFFenv

Test your configuration by performing a dry-run via

    snakemake -s Swap-seq/workflow/Snakefile --configfile config/config.json -n

Execute the workflow locally via

    snakemake -s Swap-seq/workflow/Snakefile --configfile config/config.json --cores $N

Or run it in a cluster environment via

```
snakemake \
  -s Swap-seq/workflow/Snakefile \
  --configfile config/config.json \
  --cores 1 \
  --jobs 50 \
  --max-jobs-per-second 2 \
  --max-status-checks-per-second 2 \
  --restart-times 1 \
  -k \
  --cluster "sbatch -n 1 -c 1 --mem 4G -t 1:00:00 -p owners,engreitz,normal -J VFF_{rule} -o log/{rule}_{wildcards} -e log/{rule}_{wildcards}"
```

# Differences from Variant-EFFECTS

| | Variant-EFFECTS | Swap-seq |
|---|---|---|
| Readout | the edited locus | a barcode cassette |
| `MappingSequence` | variant plus 3–5 bp of context | the full 71 bp cassette |
| CRISPResso2 alignment | default gap penalty | `--needleman_wunsch_gap_open -200` to prevent the software from inserting gaps in the barcode cassettes and cause misalignment |
| Editing efficiency | measured at the edited locus | derived from non-editing barcode-only cassettes |
| Wild-type drag correction | `V2_effect_size_adjustment.py` The WT baseline of Variant-EFFECTS is measured by the distribution of unedited alleles, which are from two sources -- homozygous, WT cells and heterozygous, edited cells. Therefore, the WT baseline is "dragged" by the heterozygous cells and requires additional corrections.| **removed.** The WT of Swap-seq is from non-editing barcode cassettes. WT cells or unedited alleles in heteterozygous, edited cells are not sequenced because they don't have a barcode cassettes and therefore do not contribute to WT distribution. For details, see Method: Swap-seq PPIF tiling deletion screen design of the manuscript.|

## Outputs

The *PPIF* screen's outputs are in the repo as examples:
[`example/result_tables/`](example/result_tables/) for the two key files and
[`example/result_plots/`](example/result_plots/) for the plots.

### Structure

All result files are under the `/results/` directory:

    summary/                        Key output of the pipeline: tables with results and plots.
    summary/QC/                     Replicate correlations, barcode consistency, edited-allele fraction.
    summary/Screen_results/         The per-edit volcano and barplots per targeted element.

    byPCRRep/                       Contains variant count information, mle logs, mle outputs, and PDFs quantifying
                                      effects within a given PCR replicate.

    byExperimentRepCorFilter/       Contains the same as above but the analysis is performed where PCR reps are
                                      aggregated by individual FlowFISH samples. Imposes a PCR correlation
                                      filter at r = 0.8, set by --minCorrelation in
                                      workflow/scripts/GetPCRReplicateCorrelation.R
    byExperimentRep/                See above minus the correlation filter

    aligned/                        Bam files for every sample. Includes unaligned fastq files for troubleshooting
    crispresso/                     CRISPResso2 output files
    variantCounts/                  Count files for both reference alleles and all the variants analyzed for every fastq.

### Tables in `summary/`

| File | Contents |
|---|---|
| `VariantCounts.flat.tsv.gz`, `VariantCounts.matrix.tsv.gz` | barcode counts per sample|
| `alignment.counts.tsv` | bowtie2 alignment counts |
| `PCRReplicateCorrelations.tsv`, `PCRReplicateCorrelations.LowQualSamples.tsv`, `VariationVsAlleleFrequency.tsv` | PCR-replicate correlations and variance vs allele frequency |
| `AllelicEffects.byPCRRep.ExperimentIDPCRRep.flat.tsv.gz` | MLE effects per PCR replicate |
| `AllelicEffects.byExperimentRep.ExperimentIDReplicates.flat.tsv.gz` | MLE effects per pegRNA, before the Swap-seq adjustment — the input to the Swap-seq step |
| `SwapseqPegRNAEffects.tsv` | per-pegRNA effects after the Swap-seq adjustment, one row per barcode x replicate |
| `SwapseqPegRNAEffects.unfiltered.tsv` | same as above but before applying the `min_allele_freq` threshold |
| `SwapseqEditReplicateEffects.tsv` | per-edit effects separated by replicate, barcodes averaged within an edit |
| `SwapseqEditEffects.tsv` | per-edit effect sizes and statistics |
| `QC/DesiredVariants.PCRReplicateCorrelations.tsv`, `.PCRReplicateVariantCV.tsv`, `.RData` | per-barcode frequency summaries used for the QC panels |

### Key files

| File | Contents |
|---|---|
| [`summary/SwapseqPegRNAEffects.tsv`](example/result_tables/SwapseqPegRNAEffects.tsv) | Per-pegRNA effect sizes, one row per barcode x replicate |
| [`summary/SwapseqEditEffects.tsv`](example/result_tables/SwapseqEditEffects.tsv) | Per-edit effect sizes and statistics |

The MLE reports an effect per **pegRNA** (equivalently, per barcode: one barcode cassette
per pegRNA). The `swapseq_edit_effects` rule applies the Swap-seq adjustment at the pegRNA level
— anchoring allele frequencies on the barcode-only spike-ins, converting them to editing
efficiency, then applying the heterozygous correction — and writes the result to
`SwapseqPegRNAEffects.tsv`. Then barcodes are averaged within
an edit and each edit is tested against no effect with a one-sample t-test and Benjamini-Hochberg correction, giving
`SwapseqEditEffects.tsv`.

### Plots

Drawn from the Swap-seq tables above. The *PPIF* versions of all of these are in
[`example/result_plots/`](example/result_plots/).

The figures ask for Arial, which matplotlib does not ship. Point `SWAPSEQ_ARIAL_TTF` at an
Arial `.ttf` to match the manuscript; otherwise they fall back to the default sans-serif.

| Plot in `summary/QC/` | Contents |
|---|---|
| `replicate_correlation_biorep_perEdit.pdf` | BioRep 1 vs BioRep 2 scatter, **one dot per edit**. Allele-frequency filtered, BR1 FFrep4 and 5 excluded |
| `replicate_correlation_biorep_perPegRNA.pdf` | the same, one dot per pegRNA |
| `replicate_correlation_ffreps_perEdit.pdf`, `..._perPegRNA.pdf` | one page per pair of FlowFISH replicates |
| `replicate_correlation_all.pdf`, `..._filtered.pdf` | the same across every replicate pair, as a heatmap |
| `barcode_consistency_min{2,3,4}.pdf` | test barcode consistency by randomly splitting the barcodes to two groups and plotting correlation of the mean effects of each group |
| `pool_allele_frequency_by_biorep.pdf` | total edited-allele fraction per biological replicate, mean over FlowFISH reps, 95% CI, dots = the reps |
| `pool_allele_frequency_distribution_combined.pdf` | distribution of each barcode's mean allele frequency with the threshold marked, and a second x axis giving the equivalent cells sorted per BioRep |
| `DesiredVariants.replicateCorrelations.pdf` | barcode-frequency correlation between PCR replicates of the same FlowFISH rep and bin, plus variability vs abundance |

| Plot in `summary/Screen_results/` | Contents |
|---|---|
| `volcano.pdf` | per-edit effect against BH-adjusted significance with 95% CI error bars, coloured by design category |
| `horizontal_barplots/<element>.pdf` | every edit tiling one element: bar = mean over replicates, error = 95% CI, dots = replicates; red when BH < 0.001, grey otherwise|

## License

GPL-3.0, see [LICENSE](LICENSE).
