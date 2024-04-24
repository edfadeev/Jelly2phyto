#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mem=100GB
#SBATCH --partition=basic
#SBATCH --mail-user=ALL
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=12:00:00

module load conda

#activate kraken2
conda activate kraken2-2.1.3

export PATH=/home/user/fadeev/anvio-resources/KrakenTools:/lisc/user/fadeev/.local/bin:$PATH

K_dbs=/scratch/oceanography/efadeev/anvio-resources/kraken2/k2_pluspf_20240112

i=0
while read line
do
chunks[$i]="$line"
i=$((i+1))
done < $WORKDIR/07_TAXONOMY/fastq.list

for SAMPLE in ${chunks[${SLURM_ARRAY_TASK_ID}]};
do

kraken2 --threads 24 --db $K_dbs --output ${WORKDIR}/07_TAXONOMY/${SAMPLE}_kraken2.out \
--report ${WORKDIR}/07_TAXONOMY/${SAMPLE}_kraken2.report \
--paired ${WORKDIR}/01_QC/${SAMPLE}-QUALITY_PASSED_R1.fastq ${WORKDIR}/01_QC/${SAMPLE}-QUALITY_PASSED_R2.fastq

done