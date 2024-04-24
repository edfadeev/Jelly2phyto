#!/bin/bash
#
#SBATCH --cpus-per-task=24
#SBATCH --mem=200GB
#SBATCH --partition=basic
#SBATCH --mail-user=ALL
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=48:00:00

#add kaiju to path
export PATH=/scratch/oceanography/efadeev/anvio-resources/kaiju/bin:$PATH

K_dbs=/scratch/oceanography/efadeev/anvio-resources/kaiju_dbs

  
i=0
while read line
do
chunks[$i]="$line"
i=$((i+1))
done < $WORKDIR/samples.txt

for SAMPLE in ${chunks[${SLURM_ARRAY_TASK_ID}]};
do
#file=$(basename $SAMPLE)

cp $WORKDIR/01_QC/$SAMPLE-QUALITY_PASSED_R1.fastq.gz $TMPDIR/R1.fastq.gz
cp $WORKDIR/01_QC/$SAMPLE-QUALITY_PASSED_R2.fastq.gz $TMPDIR/R2.fastq.gz

pushd $TMPDIR

#against Refseq
kaiju -t $K_dbs/kaiju_db_nr_euk_2023-05-10/nodes.dmp \
-f $K_dbs/kaiju_db_nr_euk_2023-05-10/kaiju_db_nr_euk.fmi \
-i R1.fastq.gz \
-j R2.fastq.gz
-o kaiju.out \
-z 24 \
-v

#add taxon name
kaiju-addTaxonNames -t $K_dbs/kaiju_db_nr_euk_2023-05-10/nodes.dmp \
-n $K_dbs/kaiju_db_nr_euk_2023-05-10/names.dmp \
-i kaiju.out \
-o kaiju_names.out \
-r superkingdom,phylum,order,class,family,genus,species

popd

mv $TMPDIR/kaiju.out $WORKDIR/07_TAXONOMY/$SAMPLE_kaiju.out
mv $TMPDIR/kaiju_names.out $WORKDIR/07_TAXONOMY/$SAMPLE_kaiju_names.out

done