## Shawn Cai
## 2025/06/09
## Swap-seq barcode cassette library, PPIF locus, THP-1
## 2 BioReps x 7 FlowFISH reps x (4 sorted bins + input) x 16 PCR reps = 924 libraries

PROJECT=$OAK/Projects/VariantEditing/FF/250609_Swap-seq_TDL
BCLDIR=$OAK/Projects/SequencingRuns/swap-seq_TDL_novogene_250609/250607_lh00134_0730_A232C5LLT3

cd $PROJECT
git clone https://github.com/EngreitzLab/Swap-seq.git
mkdir -p $PROJECT/fastq/ $PROJECT/log/ $PROJECT/results/

## Copy this example configuration into place and edit the paths
cp -r Swap-seq/example/config    $PROJECT/config
cp -r Swap-seq/example/sortParams $PROJECT/sortParams


########################################################################
## Demultiplex
##

sbatch Swap-seq/example/demultiplex.batch


########################################################################
## Trim reads to the 71 bp cassette
##

bash Swap-seq/example/trim_fastq.sh


########################################################################
## Run the pipeline

conda activate VFFenv

## Dry run first
snakemake \
  -s Swap-seq/workflow/Snakefile \
  --configfile config/config.json -n

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
