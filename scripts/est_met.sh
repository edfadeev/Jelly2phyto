#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --partition=basic
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=96:00:00

module load conda
conda activate anvio-dev

#anvi-estimate-metabolism -e $WORKDIR/metagenomes.list \
#-O $WORKDIR/03_CONTIGS/ --output-modes modules,hits

#anvi-estimate-metabolism -e $WORKDIR/metagenomes.list \
#-O $WORKDIR/03_CONTIGS/Metabol_meta --matrix-format

#anvi-estimate-metabolism -M $WORKDIR/metagenomes.list \
#-O $WORKDIR/03_CONTIGS/ --output-modes modules,hits

#anvi-estimate-metabolism -M $WORKDIR/metagenomes.list \
#-O $WORKDIR/03_CONTIGS/Metabol --matrix-format

anvi-estimate-metabolism -i $WORKDIR/08_BINS/combined_bins_table.txt -O $WORKDIR/08_BINS
