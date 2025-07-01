# Anvio workflow for metagenomes assembly, annotation and binning
#login with port forwarding: ssh -L 5678:localhost:5678 slurm

################################
# Define variables
################################

#start screen
screen -S MicPir

#load modules
module load conda

#Activate conda environment:
conda activate anvio-dev

WORKDIR=/lisc/scratch/oceanography/efadeev/MicPir/MicPir_sep/ && cd $WORKDIR
SCRIPT_DIR=/lisc/scratch/oceanography/efadeev/anvio_scripts/
PROJECT=MicPir_sep

#generate manually table samples_list.txt that includes the sampled ID and the paths to the fastq files

################################
# Run metagenome assembly and annotation within Anvio
################################
mkdir $WORKDIR/00_LOGS

#generate default workflow
anvi-run-workflow -w metagenomics --get-default-config $WORKDIR/coassembly_config.json

#adjust the workflow file manualy and plot it
anvi-run-workflow -w metagenomics -c $WORKDIR/coassembly_config.json --save-workflow-graph

#run workflow
anvi-run-workflow -w metagenomics -c $WORKDIR/coassembly_config.json \
--additional-params --cluster-config $WORKDIR/cluster_config.json --jobs 200 \
--cluster "sbatch --output=$WORKDIR/00_LOGS/%x-%j.out --partition={cluster.slurm_partition} \
--mem={cluster.mem_gb}GB --cpus-per-task={cluster.nodes} --time=1-{cluster.runtime}"

################################
# Taxonomic classification of raw reads and genes
################################
mkdir 07_TAXONOMY

#kraken2 classification
#adjust array and run sbatch
sbatch --array=0-26 -D `pwd` --export=ALL,WORKDIR=$WORKDIR \
--job-name "fastq_tax_Kraken2" $SCRIPT_DIR/kraken2_class_fastq.sh

# kaiju classification using nr_euk ref database
cd $WORKDIR/01_QC
ls *R1* | sed 's/-QUALITY_PASSED_R1.fastq.gz//g' > $WORKDIR/07_TAXONOMY/fastq.list
cd $WORKDIR

#check how many chunck are there 
n_chunks=$(awk 'END { print NR }' $WORKDIR/07_TAXONOMY/fastq.list)

#run kaiju
sbatch -D `pwd` --array=0-$n_chunks --export=ALL,WORKDIR=$WORKDIR \
--job-name "fastq_tax_kaiju" $SCRIPT_DIR/tax_class_fastq_kaiju.sh

#merge Kaiju output into a single table
sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR kaiju2table.sh
  

#add gene taxonomy using kaiju
sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT \
--job-name "gene_tax_kaiju" $SCRIPT_DIR/gene_taxonomy_kaiju.sh

################################
# Metagenome binning using metawrap
################################
#the binning is conducted using the metawrap pipeline (https://github.com/bxlab/metaWRAP)

#create output directory
mkdir $WORKDIR/08_BINS/
mkdir $WORKDIR/08_BINS/fastqs
  
#TO start the pipeline we need unzipped fastq files with 1/2 ending (instead of 1/R2) 
#and the contigs of the coassembly

#rename fastq files
gzip -d $WORKDIR/01_QC/*.fastq.gz
  
rename -- _R _ $WORKDIR/01_QC/*.fastq


#start the metawrap pipeline
sbatch -D `pwd` --array=0-8 --export=ALL,WORKDIR=$WORKDIR --job-name "metawrap_binning" $REPO_DIR/scripts/metawrap_binning.sh

#summarize bins
mkdir $WORKDIR/09_SUMMARY

sbatch -D `pwd` --array=0-8 --export=ALL,WORKDIR=$WORKDIR --job-name "bins2anvio" $REPO_DIR/scripts/bins2anvio.sh

#combine into one table all the bins
touch $WORKDIR/09_SUMMARY/combined_bins_summary.tmp

for metaG in T3J T4J T3C T4C SC SJ ML M T0;
do
sed 's/bin\\./bin_/g' $WORKDIR/09_SUMMARY/$metaG/bins_summary.txt |
  awk -F' ' -v metaG=$metaG 'BEGIN{OFS="\t"};NR>1{print metaG, metaG"_"$1, $0}' - > \
  $WORKDIR/09_SUMMARY/${metaG}_bins_summary.tmp

cat $WORKDIR/09_SUMMARY/${metaG}_bins_summary.tmp >> $WORKDIR/09_SUMMARY/combined_bins_summary.tmp
done

awk 'BEGIN{OFS="\t"}; NR==1{print \
"group","db_name","bins","total_length","num_contigs",
"N50","GC_content","percent_completion",
"percent_redundancy","t_domain","t_phylum",
"t_class","t_order","t_family",
"t_genus","t_species"}; NR>1{print}' \
$WORKDIR/09_SUMMARY/combined_bins_summary.tmp > $WORKDIR/09_SUMMARY/combined_bins_summary.txt

rm $WORKDIR/09_SUMMARY/*.tmp 

###############################################################################################
#estimate metabolism for each bin
###############################################################################################
#make a combined list of all bins
for metaG in T3J T4J T3C T4C SC SJ ML M T0;
do
sed 's/\./_/g' $WORKDIR/08_BINS/$metaG/merged/metawrap_70_10_bins.stats |
  awk -F' ' -v metaG=$metaG -v wd=$WORKDIR \
'NR>1{print metaG"_"$1, $1, metaG"_metawrap", wd"/06_MERGED/"metaG"/PROFILE.db", wd"/03_CONTIGS/"metaG"-contigs.db" }' -|
  awk 'gsub(" ","\t")' -> $WORKDIR/08_BINS/$metaG/merged/internal-genomes.txt

cat $WORKDIR/08_BINS/$metaG/merged/internal-genomes.txt >> $WORKDIR/08_BINS/combined_bins.txt
done


awk 'BEGIN{OFS="\t"}; NR==1{print "name","bin_id","collection_id", "profile_db_path","contigs_db_path"}; NR>1 {print}' \
$WORKDIR/08_BINS/combined_bins.txt > $WORKDIR/08_BINS/combined_bins_table.txt

#estimate metabolism
sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR --job-name "est_met_bins" $REPO_DIR/anvio_scripts/est_met.sh

#run CAZy annotation
sbatch -D `pwd` --array=0-8 --export=ALL,WORKDIR=$WORKDIR --job-name "CAZy_annotation" $REPO_DIR/scripts/run_cazy.sh

###############################################################################################
#test Metabolic enrichment 
###############################################################################################
#Jelly-Control
#generate group-txt
grep "T3J\\|T4J\\|T3C\\|T4C" $WORKDIR/08_BINS/combined_bins_table.txt |\
sed 's/\t.*C_metawrap/\tControl/g' - |\
sed 's/\t.*J_metawrap/\tJelly/g' -|\
awk -F '\t' 'BEGIN{OFS="\t"}; {print $1, $2}' - |\
awk 'BEGIN{OFS="\t"}; NR==1{print "name","group"}; NR>1 {print}' >\
$WORKDIR/08_BINS/BINS_Jelly_Control_groups.txt

#run enrichment
anvi-compute-metabolic-enrichment -M $WORKDIR/08_BINS/BINS_modules.txt \
-G $WORKDIR/08_BINS/BINS_Jelly_Control_groups.txt \
-o $WORKDIR/08_BINS/Jelly_control_met_enrichment.txt

#SJ-SC
#generate group-txt
grep "SJ\\|SC" $WORKDIR/08_BINS/combined_bins_table.txt |\
sed 's/\t.*C_metawrap/\tSC/g' - |\
sed 's/\t.*J_metawrap/\tSJ/g' -|\
awk -F '\t' 'BEGIN{OFS="\t"}; {print $1, $2}' - |\
awk 'BEGIN{OFS="\t"}; NR==1{print "name","group"}; NR>1 {print}' >\
$WORKDIR/08_BINS/BINS_SJ_SC_groups.txt

#run enrichment
anvi-compute-metabolic-enrichment -M $WORKDIR/08_BINS/BINS_modules.txt \
-G $WORKDIR/08_BINS/BINS_SJ_SC_groups.txt \
-o $WORKDIR/08_BINS/SJ_SC_met_enrichment.txt

###############################################################################################
#test functional enrichment (COG20)
###############################################################################################
#generate combined COG20 definitions table
awk 'BEGIN{OFS="\t"; print "group","gene_callers_id","source","accession","function","e_value"}' \
> $WORKDIR/03_CONTIGS/COG20.txt

for metaG in SC SJ T3J T4J T3C T4C ML M T0;
do

anvi-export-functions -c $WORKDIR/03_CONTIGS/${metaG}-contigs.db \
-o $WORKDIR/03_CONTIGS/${metaG}-COG20.txt \
--annotation-sources COG20_CATEGORY,COG20_FUNCTION

awk -v metaG=$metaG 'BEGIN{OFS="\t"};\
NR>1{print metaG,$0}' $WORKDIR/03_CONTIGS/${metaG}-COG20.txt >> \
$WORKDIR/03_CONTIGS/COG20.txt

done

Rscript generate_COG20_def_table.R



#Jelly-Control
#generate bins table 
grep "T3J\\|T4J\\|T3C\\|T4C" $WORKDIR/08_BINS/combined_bins_table.txt | \
awk  'BEGIN{OFS="\t"}; NR==1{print "name","bin_id",\
"collection_id", "profile_db_path","contigs_db_path"}; \
NR>1 {print}' - > $WORKDIR/08_BINS/Jelly_control_bins.txt

#SJ-SC
#generate bins table 
grep "SJ\\|SC" $WORKDIR/08_BINS/combined_bins_table.txt | \
awk  'BEGIN{OFS="\t"}; NR==1{print "name","bin_id",\
"collection_id", "profile_db_path","contigs_db_path"}; \
NR>1 {print}' - > $WORKDIR/08_BINS/SJ_SC_bins.txt

#estimate functional enrichment
sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR --job-name "comp_fun_enr" \
$REPO_DIR/anvio_scripts/comp_fun_enr.sh