#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --partition=basic
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=12:00:00

module load conda
conda activate anvio-dev

anvi-compute-functional-enrichment-across-genomes -G \
$WORKDIR/08_BINS/BINS_Jelly_Control_groups.txt \
-o $WORKDIR/08_BINS/Jelly_control_COG20_FUN_enrichment.txt \
--annotation-source COG20_FUNCTION -i $WORKDIR/08_BINS/Jelly_control_bins.txt

anvi-compute-functional-enrichment-across-genomes -G \
$WORKDIR/08_BINS/BINS_SJ_SC_groups.txt \
-o $WORKDIR/08_BINS/SJ_SC_COG20_FUN_enrichment.txt \
--annotation-source COG20_FUNCTION -i $WORKDIR/08_BINS/SJ_SC_bins.txt
