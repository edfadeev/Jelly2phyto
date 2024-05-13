# Anvio workflow for metagenomes assembly, annotation and binning
#login with port forwarding: ssh -L 5678:localhost:5678 slurm

################################
# Define variables
################################
export PATH=/home/user/fadeev/anvio-resources/DAS_Tool/:/home/user/fadeev/anvio-resources/pullseq/src/:/lisc/user/fadeev/anvio-resources/InterProScanParser:/home/user/fadeev/anvio-resources/KrakenTools:/lisc/user/fadeev/.local/bin:$PATH

#Define working directory
WORKDIR=/scratch/oceanography/efadeev/Microcosm_Piran/MicPir_coasem/ && cd $WORKDIR
PROJECT=MicPir_coasm

#generate manually table samples_list.txt that includes the sampled ID and the paths to the fastq files

mkdir $WORKDIR/00_LOGS
################################
# Run metagenome assembly and annotation within Anvio
################################
#start screen
screen -S MicPir

#load modules
module load conda

#Activate conda environment:
conda activate anvio-dev

#Define working directory
WORKDIR=/scratch/oceanography/efadeev/Microcosm_Piran/MicPir_coasem/ && cd $WORKDIR
PROJECT=coassembly

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
# Classify raw reads using Kraken2
################################
mkdir 07_TAXONOMY

#adjust array and run sbatch
sbatch --array=0-26 -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "fastq_tax_classification" ./scripts/kraken2_class_fastq.sh

################################
#add gene taxonomy using kaiju
################################
sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "gene_tax_kaiju" scripts/gene_taxonomy_kaiju.sh

################################
# Metagenome binning using metawrap
################################
#the binning is conducted using the metawrap pipeline (https://github.com/bxlab/metaWRAP)

#create output directory
mkdir $WORKDIR/08_BINS/

#TO start the pipeline we need unzipped fastq files with 1/2 ending (instead of 1/R2) 
#and the contigs of the coassembly

#rename fastq files
gzip -d $WORKDIR/01_QC/*.fastq.gz
  
rename -- _R _ $WORKDIR/01_QC/*.fastq


#start the metawrap pipeline
sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "metawrap_binning" $REPO_DIR/scripts/metawrap_binning.sh

#import bins into anvio
sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "metawrap2anvio" $REPO_DIR/scripts/metawrap2anvio.sh

#estimate metabolism for each bin
sbatch -D `pwd` --export=ALL,WORKDIR=$WORKDIR,PROJECT=$PROJECT --job-name "est_metabolism_bins" $REPO_DIR/anvio_scripts/est_metabolism_bins.sh


anvi-interactive -p $WORKDIR/06_MERGED/$PROJECT/PROFILE.db \
-c $WORKDIR/03_CONTIGS/$PROJECT-contigs.db \
--additional-layers $WORKDIR/09_SUMMARY/bins_TAXONOMY.txt -C metawrap \
--server-only -P 6785

