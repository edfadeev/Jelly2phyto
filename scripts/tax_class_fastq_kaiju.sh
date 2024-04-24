#!/bin/bash
#SBATCH --cpus-per-task=12
#SBATCH --mem=70GB
#SBATCH --partition=basic
#SBATCH --mail-user=ALL
#SBATCH --time=24:00:00
#SBATCH --output=./00_LOGS/%x-%j.out

module load conda

#load kaiju 
conda activate kaiju-1.10.1

#path to kaiju dbs
K_dbs=/lisc/scratch/oceanography/efadeev/anvio-resources/kaiju_dbs/

#generate array with protein chunk names
i=0
while read line
do
chunks[$i]="$line"
i=$((i+1))
done < $WORKDIR/07_TAXONOMY/fastq.list

for SAMPLE in ${chunks[${SLURM_ARRAY_TASK_ID}]}; do

kaiju -t $K_dbs/kaiju_db_nr_euk_2023-05-10/nodes.dmp \
-f $K_dbs/kaiju_db_nr_euk_2023-05-10/kaiju_db_nr_euk.fmi \
-i $WORKDIR/01_QC/${SAMPLE}-QUALITY_PASSED_R1.fastq.gz \
-j $WORKDIR/01_QC/${SAMPLE}-QUALITY_PASSED_R2.fastq.gz \
-o $WORKDIR/07_TAXONOMY/${SAMPLE}_nr_euk_tax.out \
-z 12 \
-v

#create krona for each sample
kaiju2krona -t $K_dbs/kaiju_db_nr_euk_2023-05-10/nodes.dmp \
-n $K_dbs/kaiju_db_nr_euk_2023-05-10/names.dmp \
-i $WORKDIR/07_TAXONOMY/${SAMPLE}_nr_euk_tax.out \
-o $WORKDIR/07_TAXONOMY/${SAMPLE}_nr_euk_tax.krona \
-u

#load Krona
conda activate krona-2.8.1

#generate html
ktImportText -o $WORKDIR/07_TAXONOMY/${SAMPLE}_nr_euk_tax.krona.html \
$WORKDIR/07_TAXONOMY/${SAMPLE}_nr_euk_tax.krona

done