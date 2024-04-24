#!/bin/bash
#SBATCH --job-name=merge_kaiju_reports
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=ALL
#SBATCH --partition=basic
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=12:00:00

#merge Kaiju output into a single table
module load conda
#load kaiju 
conda activate kaiju-1.10.1

#path to kaiju dbs
K_dbs=/lisc/scratch/oceanography/efadeev/anvio-resources/kaiju_dbs/

kaiju2table -t $K_dbs/kaiju_db_nr_euk_2023-05-10/nodes.dmp \
-n $K_dbs/kaiju_db_nr_euk_2023-05-10/names.dmp \
-r phylum -o $WORKDIR/07_TAXONOMY/kaiju_nr_euk_merged_Phylum.out \
-l superkingdom,phylum,class,order,family,genus,species \
$WORKDIR/07_TAXONOMY/*.out

kaiju2table -t $K_dbs/kaiju_db_nr_euk_2023-05-10/nodes.dmp \
-n $K_dbs/kaiju_db_nr_euk_2023-05-10/names.dmp \
-r class -o $WORKDIR/07_TAXONOMY/kaiju_nr_euk_merged_Class.out \
-l superkingdom,phylum,class,order,family,genus,species \
$WORKDIR/07_TAXONOMY/*.out

kaiju2table -t $K_dbs/kaiju_db_nr_euk_2023-05-10/nodes.dmp \
-n $K_dbs/kaiju_db_nr_euk_2023-05-10/names.dmp \
-r order -o $WORKDIR/07_TAXONOMY/kaiju_nr_euk_merged_Order.out \
-l superkingdom,phylum,class,order,family,genus,species \
$WORKDIR/07_TAXONOMY/*.out

kaiju2table -t $K_dbs/kaiju_db_nr_euk_2023-05-10/nodes.dmp \
-n $K_dbs/kaiju_db_nr_euk_2023-05-10/names.dmp \
-r family -o $WORKDIR/07_TAXONOMY/kaiju_nr_euk_merged_Family.out \
-l superkingdom,phylum,class,order,family,genus,species \
$WORKDIR/07_TAXONOMY/*.out

kaiju2table -t $K_dbs/kaiju_db_nr_euk_2023-05-10/nodes.dmp \
-n $K_dbs/kaiju_db_nr_euk_2023-05-10/names.dmp \
-r genus -o $WORKDIR/07_TAXONOMY/kaiju_nr_euk_merged_Genus.out \
-l superkingdom,phylum,class,order,family,genus,species \
$WORKDIR/07_TAXONOMY/*.out


