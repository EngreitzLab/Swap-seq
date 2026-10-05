## Plotting.
##
## Outputs are split in two:
##   results/summary/QC/             
##   results/summary/Screen_results/


############################################################
## QC



## Barcode-frequency correlation between PCR replicates, plus variability vs abundance.
rule plot_genotyping_stats:
	input:
		variantCounts="results/summary/VariantCounts.flat.tsv.gz",
		samplesheet="SampleList.snakemake.tsv"
	output:
		"results/summary/QC/DesiredVariants.RData"
	params:
		codedir=config['codedir']
	shell:
		"Rscript {params.codedir}/workflow/scripts/PlotVariantCounts.R --variantCounts {input.variantCounts} --samplesheet {input.samplesheet} --outbase results/summary/QC/DesiredVariants"

## Swap-seq QC: replicate-correlation heatmaps, split-half barcode consistency, and the edited-allele fraction per biological replicate.
rule plot_swapseq_qc:
	input:
		replicates='results/summary/SwapseqEditReplicateEffects.tsv',
		pegrnas='results/summary/SwapseqPegRNAEffects.tsv',
		pegrnas_unfiltered='results/summary/SwapseqPegRNAEffects.unfiltered.tsv'
	output:
		biorep_edit='results/summary/QC/replicate_correlation_biorep_perEdit.pdf',
		biorep_pegrna='results/summary/QC/replicate_correlation_biorep_perPegRNA.pdf',
		ffreps_edit='results/summary/QC/replicate_correlation_ffreps_perEdit.pdf',
		ffreps_pegrna='results/summary/QC/replicate_correlation_ffreps_perPegRNA.pdf',
		correlation='results/summary/QC/replicate_correlation_all.pdf',
		consistency='results/summary/QC/barcode_consistency_min3.pdf',
		pool_biorep='results/summary/QC/pool_allele_frequency_by_biorep.pdf',
		pool_dist='results/summary/QC/pool_allele_frequency_distribution_combined.pdf'
	params:
		codedir=config['codedir']
	shell:
		"""
		bash -c '
			. $HOME/.bashrc
			conda activate VFFenv
			python {params.codedir}/workflow/scripts/plot_swapseq_qc.py \
				--replicate-effects {input.replicates} \
				--pegrna-effects {input.pegrnas} \
				--pegrna-effects-unfiltered {input.pegrnas_unfiltered} \
				--outdir results/summary/QC'
		"""


############################################################
## Screen results

## Volcano plot, plus barplots where each bar is an edit with the element and each dot a replicate.
rule plot_swapseq_results:
	input:
		edits='results/summary/SwapseqEditEffects.tsv',
		replicates='results/summary/SwapseqEditReplicateEffects.tsv'
	output:
		volcano='results/summary/Screen_results/volcano.pdf'
	params:
		codedir=config['codedir'],
		bh=config.get('bh_threshold', 0.001)
	shell:
		"""
		bash -c '
			. $HOME/.bashrc
			conda activate VFFenv
			python {params.codedir}/workflow/scripts/plot_swapseq_results.py \
				--edit-effects {input.edits} \
				--replicate-effects {input.replicates} \
				--bh-threshold {params.bh} \
				--outdir results/summary/Screen_results'
		"""


############################################################
## Per-experiment MLE effect plots

rule plot_allelic_effect_sizes:
	input:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref.txt'
	output:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref.pdf'
	params:
		codedir=config['codedir']
	shell:
		"""
		Rscript {params.codedir}/workflow/scripts/PlotMleVariantEffects.R --mleEffects {input} --outfile {output}
		"""

rule plot_allelic_effect_sizes_ignoreInputBin:
	input:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref_ignoreInputBin.txt'
	output:
		'results/{replicateDirectory}/{ExperimentIDReplicates}.effects_vs_ref_ignoreInputBin.pdf'
	params:
		codedir=config['codedir']
	shell:
		"""
		Rscript {params.codedir}/workflow/scripts/PlotMleVariantEffects.R --mleEffects {input} --outfile {output}
		"""
