#!/bin/bash
#
#SBATCH --cpus-per-task=6
#SBATCH --mem=50GB
#SBATCH --partition=basic
#SBATCH --mail-user=ALL
#SBATCH --time=96:00:00
#SBATCH --output=./00_LOGS/%x-%j.out

#add kaiju to path
export PATH=/home/user/fadeev/anvio-resources/kaiju/bin:$PATH

K_dbs=/home/user/fadeev/anvio-resources/kaiju_dbs/kaiju_db_progenomes_2021-03-02/

#export all the gene calls
anvi-get-sequences-for-gene-calls -c $WORKDIR/03_CONTIGS/$PROJECT-contigs.db \
-o $WORKDIR/03_CONTIGS/$PROJECT-gene_calls.fa

#add gene taxonomic assignments to genes
kaiju -t $K_dbs/nodes.dmp \
-f $K_dbs/kaiju_db_progenomes.fmi \
-i $WORKDIR/03_CONTIGS/$PROJECT-gene_calls.fa \
-o $WORKDIR/08_TAXONOMY/$PROJECT-gene_calls_progenomes_tax.out \
-z 6 \
-v

#add taxon name
kaiju-addTaxonNames -t $K_dbs/nodes.dmp \
-n $K_dbs/names.dmp \
-i $WORKDIR/08_TAXONOMY/$PROJECT-gene_calls_progenomes_tax.out \
-o $WORKDIR/08_TAXONOMY/$PROJECT-gene_calls_progenomes_names.out \
-r superkingdom,phylum,order,class,family,genus,species

#add to anvio db
anvi-import-taxonomy-for-genes -i $WORKDIR/08_TAXONOMY/$PROJECT-gene_calls_progenomes_names.out \
-c $WORKDIR/03_CONTIGS/$PROJECT-contigs.db \
-p kaiju --just-do-it
