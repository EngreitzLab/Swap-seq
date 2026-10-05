# Swap-seq per-edit effect sizes.
#
# The MLE reports an effect per BARCODE. This turns that into an effect per EDIT:
# anchor allele frequencies on the barcode-only spike-in cassettes, convert to editing
# efficiency, apply the heterozygous correction, aggregate barcodes into edits, then
# one-sample t-test and Benjamini-Hochberg.


rule swapseq_edit_effects:
	input:
		flat='results/summary/AllelicEffects.byExperimentRep.ExperimentIDReplicates.flat.tsv.gz'
	output:
		edits='results/summary/SwapseqEditEffects.tsv',
		pegrnas='results/summary/SwapseqPegRNAEffects.tsv',
		pegrnas_unfiltered='results/summary/SwapseqPegRNAEffects.unfiltered.tsv',
		replicates='results/summary/SwapseqEditReplicateEffects.tsv'
	params:
		codedir=config['codedir'],
		variantInfo=config['variant_info'],
		spikein=0.5 * config.get('spikein_plasmid_fraction', 0.15), # at low MOI, each cells get one lentivirus which integrates in a hemizygious way
		exclude=config.get('exclude_replicates', '') # BioRep:FFRep pairs to drop
	shell:
		"""
		bash -c '
			. $HOME/.bashrc
			conda activate VFFenv
			python {params.codedir}/workflow/scripts/swapseq_edit_effects.py \
				{input.flat} \
				{output.edits} \
				--pegrna-output {output.pegrnas} \
				--unfiltered-pegrna-output {output.pegrnas_unfiltered} \
				--replicate-output {output.replicates} \
				--variant-info {params.variantInfo} \
				--spikein-allele-fraction {params.spikein} \
				--exclude-replicates "{params.exclude}"'
		"""
