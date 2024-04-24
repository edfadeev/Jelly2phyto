#!/bin/bash
#SBATCH --job-name=metawrap_binning
#SBATCH --cpus-per-task=32
#SBATCH --mem=200GB
#SBATCH --mail-user=ALL
#SBATCH --partition=basic
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=48:00:00

module load conda 

#load anvio
conda activate metawrap-1.3.2

GRPS=(T3J T4J T3C T4C SC SJ ML M0 T0)

for metaG in ${GRPS[${SLURM_ARRAY_TASK_ID}]};
do

#gunzip -c $WORKDIR/01_QC/*$metaG*.fastq.gz > $WORKDIR/08_BINS/fastqs/

#rename fastq files
#rename -- _R _ $WORKDIR/08_BINS/fastqs/*$metaG*.fastq

#binning
metawrap binning -a $WORKDIR/02_FASTA/${metaG}/${metaG}-contigs-prefix-formatted-only.fa \
-o $WORKDIR/08_BINS/${metaG} --concoct --metabat2 --maxbin2 -t 32 -m 200 $WORKDIR/08_BINS/fastqs/*${metaG}*.fastq

#refinment
metawrap bin_refinement -o $WORKDIR/08_BINS/${metaG}/merged -t 32 -m 200 \
-A $WORKDIR/08_BINS/${metaG}/metabat2_bins/ -B $WORKDIR/08_BINS/${metaG}/maxbin2_bins/ -C $WORKDIR/08_BINS/${metaG}/concoct_bins/ -c 70 -x 10

done
