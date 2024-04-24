#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --mail-user=dr.eduard.fadeev@gmail.com
#SBATCH --partition=basic
#SBATCH --output=./00_LOGS/%x-%j.out
#SBATCH --time=24:00:00

module load conda
conda activate anvio-dev


GRPS=(SC SJ T3J T4J T3C T4C ML M T0)

for metaG in ${GRPS[${SLURM_ARRAY_TASK_ID}]};
do

###############################################################################################
#Import to anvio the results of metwrap (Completness > 70, Contamination < 10))
###############################################################################################
sed -i 's/\./_/g' $WORKDIR/08_BINS/${metaG}/merged/metawrap_70_10_bins.contigs

anvi-import-collection --collection-name ${metaG}_metawrap \
--pan-or-profile-db $WORKDIR/06_MERGED/${metaG}/PROFILE.db \
--contigs-db $WORKDIR/03_CONTIGS/${metaG}-contigs.db \
--contigs-mode $WORKDIR/08_BINS/${metaG}/merged/metawrap_70_10_bins.contigs


#generate summary of the bins 
anvi-summarize -p $WORKDIR/06_MERGED/${metaG}/PROFILE.db \
-c $WORKDIR/03_CONTIGS/${metaG}-contigs.db -C ${metaG}_metawrap \
-o $WORKDIR/09_SUMMARY/${metaG}
  
#add taxonomy to each bin for visualization
awk '{print $1,$10,$11,$12, $13, $14}' $WORKDIR/09_SUMMARY/${metaG}/bins_summary.txt | \
awk '{if (NR!=1) {print}}' - | \
awk -F' ' 'BEGIN{OFS=FS}$2 == "" {$2 = $3 = $4 = $5 = $6 ="UNKNOWN"} 1' - | \
awk -F' ' 'BEGIN{OFS=FS}$3 == "" {$3 = $4 = $5 = $6 ="UNKNOWN"} 1' - | \
awk -F' ' 'BEGIN{OFS=FS}$4 == "" {$4 = $5 = $6 ="UNKNOWN"} 1' - | \
awk -F' ' 'BEGIN{OFS=FS}$5 == "" {$5 = $6 ="UNKNOWN"} 1' - | \
awk -F' ' 'BEGIN{OFS=FS}$6 == "" {$6 ="UNKNOWN"} 1' - > $WORKDIR/09_SUMMARY/${metaG}/bins_tax.txt

#add source of each bin
awk '{print $1, $8}' $WORKDIR/08_BINS/${metaG}/merged/metawrap_70_10_bins.stats | \
awk '{if (NR!=1) {print}}' - | sed 's/bin\./bin_/' -  > $WORKDIR/09_SUMMARY/${metaG}/metawrap_70_10_bins.source

#merge tables
join -j 1 --header <(sort $WORKDIR/09_SUMMARY/${metaG}/metawrap_70_10_bins.source) <(sort $WORKDIR/09_SUMMARY/${metaG}/bins_tax.txt)| \
awk 'gsub(" ","\t")' - |
  awk 'BEGIN{print "item_name\tbinner\tClass\tOrder\tFamily\tGenus\tSpecies"}1' - > $WORKDIR/09_SUMMARY/${metaG}/bins_meta.txt

done
