#!/bin/bash
#SBATCH --job-name=nexflow_achro
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=1G

module load Nextflow/23.10.0
module load Biopython/1.79-foss-2022a
module load SciPy-bundle/2022.05-foss-2022a
module load matplotlib/3.5.2-foss-2022a
module load CD-HIT/4.8.1-GCC-11.3.0
module load BLAST+/2.13.0-gompi-2022a

#screen -dmS achro_clean bash -c ./achro_clean.sh
now=`date +"%s"`
foi=$VSC_SCRATCH_VO_USER/achromatium/full_achromatium

cd $foi
#echo logfile at $foi/run_report_achro_clean_$now.nflog

nextflow run $VSC_SCRATCH_VO_USER/achromatium/achromatium_pipeline/main.nf -resume\
    -profile hydra \
    --data $foi/genomes_clean \
    --outFolder $foi/og_clean \
    --orthoGroupSeqs $foi/og_clean/orthofinder/Results_Nov29 \
    --structures  $foi/og_clean/structures --predictRemaining true \
    --minCoverage 0.8 \
    --outGroup allochromatiumVinosum \
    --foldseekClusters $foi/og_clean/structures 
#    |& tee  $foi/run_report_achro_clean_$now.nflog
#sessionName=$(sed -n '2s/.*\[\(.*\)\].*/\1/p' $foi/run_report_achro_clean_$now.nflog)
#nextflow log | grep $sessionName >> $foi/run_report_achro_clean_$now.nflog
# --predictRemaining true