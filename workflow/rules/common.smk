from snakemake.utils import validate
import pandas as pd
import os
import glob

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
# singularity: "docker://continuumio/miniconda3"



###########################################################################################
##### load config and sample sheets #####

# configfile: "config/config.yaml"   ## Read from command line instead
# validate(config, schema="../schemas/config.schema.yaml")


def find_fastq_files(samplesheet, fastqdir):
	## Adds columns 'fastqR1' and 'fastqR2' to the sample sheet, only if they do not already exist

	if single_end:
		reads = ["1"]
	else:
		reads = ["1","2"]

	for read in reads:
		colName = 'fastqR' + read
		if not colName in samplesheet.columns:
			samplesheet[colName] = ""
			for i in samplesheet.index:
				currSample = samplesheet.at[i,'SampleID']
				## Try looking for the bcl2fastq output file format
				file = glob.glob("{}_*_R{}_*fastq.gz".format(os.path.join(fastqdir, currSample), read))
				if len(file) == 0: ## Try looking for the barcode-splitter output format
					file = glob.glob("{}-read-{}.fastq.gz".format(os.path.join(fastqdir, currSample), read))
				if len(file) > 1:
					raise ValueError("Found more than one FASTQ file for sample :" + currSample)
				elif len(file) == 0:
					print("Warning: Could not find FASTQ file for read " + read + " and sample: " + currSample)
				elif len(file) == 1:
					samplesheet.at[i,colName] = file[0]


	return samplesheet


def find_sort_params_files(samplesheet):
	## If the user provided 'sortParamsFile' in the samplesheet, use that. Otherwise, look in the sortParams directory
	if not 'sortParamsFile' in samplesheet.columns:
		samplesheet['sortParamsFile'] = [os.path.join(config['sortparamsdir'], str(row['Batch']) + "_" + str(row['SampleNumber']) + ".csv") for idx, row in samplesheet.iterrows()]
		samplesheet.loc[samplesheet['SampleNumber'].isnull(),'sortParamsFile'] = ""
	return samplesheet


def add_experiment_names(samplesheet):
	if ('ExperimentIDReplicates' in samplesheet.columns) or ('ExperimentID' in samplesheet.columns) or ('ExperimentIDPCRRep' in samplesheet.columns):
		print("Warning: ExperimentID columns found and will be overwritten in the sample sheet")

	## Experiments at the level of PCR replicates
	cols_ExperimentIDPCRRep = keyCols + repCols + ['PCRRep']
	cols_ExperimentIDReplicates = keyCols + repCols

	s = samplesheet[cols_ExperimentIDPCRRep].drop_duplicates()
	s['ExperimentIDPCRRep'] = ['-'.join([str(v) for v in list(row.values)]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	## Experiments at the level of specified replicate columns, before spike correction
	s = samplesheet[cols_ExperimentIDReplicates].drop_duplicates()
	s['ExperimentIDReplicates'] = ['-'.join([str(v) for v in list(row.values)]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	## Experiments at the level of experiments (combined across replicates)
	s = samplesheet[keyCols].drop_duplicates()
	s['ExperimentID'] = ['-'.join([str(v) for v in list(row.values)]) for index,row in s.iterrows()]
	samplesheet = samplesheet.merge(s)

	return(samplesheet)


def add_outputs(samplesheet):
	samplesheet['CRISPRessoDir'] = ['results/crispresso/CRISPResso_on_{SampleID}/'.format(SampleID=row['SampleID']) for index, row in samplesheet.iterrows()]
	samplesheet['variantCountFile'] = ['results/variantCounts/{SampleID}.variantCounts.txt'.format(SampleID=row['SampleID']) for index, row in samplesheet.iterrows()]
	samplesheet['referenceAlleleFile'] = ['results/variantCounts/{SampleID}.referenceAlleles.txt'.format(SampleID=row['SampleID']) for index, row in samplesheet.iterrows()]
	if not genotyping_only:
		samplesheet['ExperimentIDPCRRep_BinCounts'] = ['results/byPCRRep/{}.bin_counts.txt'.format(e) for e in samplesheet['ExperimentIDPCRRep']]
		samplesheet['ExperimentIDReplicates_BinCounts'] = ['results/byExperimentRep/{}.bin_counts.txt'.format(e) for e in samplesheet['ExperimentIDReplicates']]
	return samplesheet


def validate_sample_sheet(samplesheet):
	print("Validating the Sample Sheet ...\n")

	if not samplesheet['SampleID'].is_unique:
		raise ValueError("SampleID column in samplesheet must not contain duplicates.")

	for col in requiredCols:
		if not col in samplesheet.columns:
			raise ValueError("Missing required column in sample sheet: " + col)

	for col in keyCols:
		if not col in samplesheet.columns:
			raise ValueError("Missing column in sample sheet that is provided in experiment_keycols in the config file: " + col)

	for col in repCols:
		if not col in samplesheet.columns:
			raise ValueError("Missing column in sample sheet that is provided in replicate_keycols in the config file: " + col)

	print("Found all Experiment Key columns. Generating comparisons for the following experiments:")
	print('\t'.join(keyCols))
	for index, row in samplesheet[keyCols].drop_duplicates().iterrows():
		print('\t'.join([str(v) for v in row.values]))


def load_sample_sheet(samplesheetFile, ampliconInfoFile, idcol='AmpliconID'):
	samplesheet = pd.read_table(samplesheetFile, dtype=str)
	samplesheet.dropna(how='all', inplace=True)
	validate_sample_sheet(samplesheet)
	samplesheet.index = samplesheet['SampleID']  ## Requires that SampleID is unique

	if not set(ampliconRequiredCols).issubset(samplesheet.columns):
		amplicons = pd.read_table(ampliconInfoFile)
		if not set(samplesheet[idcol]).issubset(amplicons[idcol]):
			raise ValueError("Some AmpliconIDs in the sample sheet are not specified in the amplicon info table.")
		if not set(ampliconRequiredCols).issubset(amplicons.columns):
			raise ValueError("Amplicon info file must contain AmpliconID AmpliconSeq GuideSpacer")
		samplesheet = samplesheet.merge(amplicons[ampliconRequiredCols])
		if not set(ampliconRequiredCols).issubset(samplesheet.columns):
			raise ValueError("Failed to merge samplesheet and amplicon info file.")

	samplesheet = find_fastq_files(samplesheet, fastqdir)
	if not genotyping_only:
		samplesheet = find_sort_params_files(samplesheet)
	samplesheet = add_experiment_names(samplesheet)
	samplesheet = add_outputs(samplesheet)
	samplesheet.index = samplesheet['SampleID']

	return samplesheet


def get_bin_list():
	binList = samplesheet['Bin'].drop_duplicates()
	if not genotyping_only:
		if  "All" not in binList.tolist():
			print("\nWARNING: Did not find any entries with Bin == 'All' (unsorted edited cells input into FlowFISH). Was this intended, or was 'All' mispelled?\n\n")
		if "Neg" not in binList.tolist():
			print("\nWARNING: Did not find any entries with Bin == 'Neg' (unedited cells used to assess sequencing error rate). Was this intended, or was 'Neg' mispelled?\n\n")
	binList = binList[(binList != "All") & (binList != "Neg") & (binList.notnull())]
	binList = [str(b) for b in list(binList)]
	binList.sort()
	print("Processing unique bins: " + ' '.join(binList))
	return(binList)

# global variables
genotyping_only = ('genotyping_only' in config) and (config['genotyping_only'].lower() == 'true')
requiredCols = ['SampleID','AmpliconID','Bin','PCRRep','ControlForAmplicon']

single_end = ('single_end' in config) and (config['single_end'].lower() == 'true')

ampliconRequiredCols = ['AmpliconID','AmpliconSeq','QuantificationWindowStart','QuantificationWindowEnd', 'ReferenceErrorThreshold']  ## To do:  Allow specifying crispresso quantification window for different amplicons
keyCols = config['experiment_keycols'].split(',')
repCols = config['replicate_keycols'].split(',')
codedir = config['codedir']
fastqdir = config['fastqdir']
sortparamsdir = config['sortparamsdir'] if not genotyping_only else None

## 'sample_list', not 'sample_sheet': in Swap-seq the Sample Sheet is the bcl2fastq
## demultiplexing sheet, a different file. Configs written against the old key name would
## otherwise fail with a bare KeyError, so say what to change.
if 'sample_list' not in config and 'sample_sheet' in config:
	raise ValueError("config key 'sample_sheet' was renamed to 'sample_list', to keep it "
	                 "distinct from the bcl2fastq Sample Sheet used for demultiplexing. "
	                 "Rename it in your config file.")

samplesheet = load_sample_sheet(config['sample_list'], config['amplicon_info'])

# only rewrite when changed: rules take this as an input
sampleListSnapshot = samplesheet.to_csv(index=False, header=True, sep='\t')
sampleListFile = "SampleList.snakemake.tsv"
if not os.path.exists(sampleListFile) or open(sampleListFile).read() != sampleListSnapshot:
	with open(sampleListFile, 'w') as fh:
		fh.write(sampleListSnapshot)

binList = get_bin_list()


#######################################################################################
####### helpers ###########

def all_input(wildcards):

	wanted_input = []

	## CRISPResso output:
	wanted_input.extend(list(samplesheet['CRISPRessoDir'].unique()))
	wanted_input.append("results/crispresso/CRISPRessoAggregate_on_Aggregate/")
	wanted_input.append("results/summary/VariantCounts.flat.tsv.gz")
	wanted_input.append("results/summary/VariantCounts.matrix.tsv.gz")
	wanted_input.append("results/summary/QC/DesiredVariants.RData")

	## Bowtie2 alignments:
	wanted_input.extend(
	 	['results/aligned/{s}/{s}.bam'.format(s=s) for s in samplesheet['SampleID'].unique()]
	  )
	wanted_input.append("results/summary/alignment.counts.tsv")
	
	## Variant counts:
	wanted_input.extend(list(samplesheet['variantCountFile'].unique()))

	## At what point do we merge in the spike-in data?

	if not genotyping_only:

		## Output files for PCR replicates (before merging spike-in data)
		wanted_input.extend(list(samplesheet['ExperimentIDPCRRep_BinCounts'].unique()))
		wanted_input.extend([
			'results/byPCRRep/{}.effects_vs_ref.pdf'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDPCRRep'].unique()
		])
		wanted_input.extend([
			'results/byPCRRep/{}.effects_vs_ref_ignoreInputBin.pdf'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDPCRRep'].unique()
		])
		wanted_input.append('results/summary/AllelicEffects.byPCRRep.ExperimentIDPCRRep.flat.tsv.gz')

		## Output files for PCR replicates (after merging spike-in data) (?)
		wanted_input.extend([])

		## Output files for replicate experiments (before merging spike-in data)
		wanted_input.extend(list(samplesheet['ExperimentIDReplicates_BinCounts'].unique()))
		wanted_input.extend([
			'results/byExperimentRep/{}.effects_vs_ref.pdf'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])
		wanted_input.extend([
			'results/byExperimentRep/{}.effects_vs_ref_ignoreInputBin.pdf'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])
		wanted_input.extend([
			'results/byExperimentRepCorFilter/{}.effects_vs_ref.pdf'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])
		wanted_input.extend([
			'results/byExperimentRepCorFilter/{}.effects_vs_ref_ignoreInputBin.pdf'.format(e) for e in samplesheet.loc[samplesheet['Bin'].isin(binList)]['ExperimentIDReplicates'].unique()
		])
		wanted_input.append('results/summary/AllelicEffects.byExperimentRep.ExperimentIDReplicates.flat.tsv.gz')

		# Swap-seq effect sizes and statistics: per pegRNA, per edit x replicate, per edit
		wanted_input.append('results/summary/SwapseqPegRNAEffects.tsv')
		wanted_input.append('results/summary/SwapseqPegRNAEffects.unfiltered.tsv')
		wanted_input.append('results/summary/SwapseqEditReplicateEffects.tsv')
		wanted_input.append('results/summary/SwapseqEditEffects.tsv')

		## Swap-seq QC and screen-result figures, drawn from the tables above.
		wanted_input.extend(["results/summary/QC/replicate_correlation_biorep_perEdit.pdf",
						"results/summary/QC/replicate_correlation_biorep_perPegRNA.pdf",
						"results/summary/QC/replicate_correlation_ffreps_perEdit.pdf",
						"results/summary/QC/replicate_correlation_ffreps_perPegRNA.pdf",
						"results/summary/QC/replicate_correlation_all.pdf",
						"results/summary/QC/barcode_consistency_min3.pdf",
						"results/summary/QC/pool_allele_frequency_by_biorep.pdf",
						"results/summary/QC/pool_allele_frequency_distribution_combined.pdf"])
		wanted_input.append("results/summary/Screen_results/volcano.pdf")



	return wanted_input


