#!/bin/bash
#
#SBATCH --job-name=metawrap_binning
#SBATCH --cpus-per-task=128
#SBATCH --mem=1000GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --partition=himem
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=168:00:00

module load conda

#load anvio
conda activate metawrap-1.3.2

#decompress fastq files
gzip -d $WORKDIR/01_QC/*.fastq.gz

rename -- _R _ $WORKDIR/01_QC/*.fastq

#binning
metawrap binning -a $WORKDIR/02_FASTA/coassembly/coassembly-contigs-prefix-formatted-only.fa \
--concoct --metabat2 --maxbin2 -t 128 -m 1000 $WORKDIR/01_QC/*.fastq \
-o $WORKDIR/08_BINS/

#refinment
metawrap bin_refinement -o $WORKDIR/08_BINS/06_MERGED/ -t 128 -m 1000 \
-A $WORKDIR/08_BINS/metabat2_bins/ -B $WORKDIR/08_BINS/maxbin2_bins/ \
-C $WORKDIR/08_BINS/concoct_bins/ -c 70 -x 10

#compress fastq files
rename -- _1 _R1 $WORKDIR/01_QC/*.fastq
rename -- _2 _R2 $WORKDIR/01_QC/*.fastq

gzip $WORKDIR/01_QC/*.fastq