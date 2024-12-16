#!/bin/bash
#SBATCH --job-name=nexflow_achro
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2G


#screen -dmS nextflowalign bash -c ./magic_hydra.sh $sdir

#for sdir in $VSC_SCRATCH_VO_USER/achromatium/ogs_for_simsapiper/* ; do
    #echo $sdir
    #cd $sdir
    #echo $pwd
    #screen -dmS nextflowalign bash -c ./$VSC_SCRATCH_VO_USER/achromatium/bin/simsapiper.sh 
    #cd $VSC_SCRATCH_VO_USER/achromatium/
#done

#sbatch ../../../achromatium_pipeline/utilscripts/simsapiper.sh ../seqs/OG0000591.toalign.fasta 

module load Nextflow/23.10.0
house=$VSC_SCRATCH_VO_USER/achromatium/full_achromatium/og_clean

#for seqsf in $house/structures/seqs/*.fasta ;do 
    #seqsf=$1
    #ogn=$(basename $seqsf .toalign.fasta | sed 's/\..*$//')
ogf="$1"
echo $ogf
while IFS= read -r line; do
    echo "Processing line: $line"
    ogn=$line
    seqsf=$house/structures/seqs/${ogn}.toalign.fasta
    sdir=$house/simsa/$ogn
    mkdir -p $sdir
    cd $sdir
    nextflow run $VSC_SCRATCH_VO_USER/simsapiper/simsapiper.nf -resume\
        -profile hydralocal \
        --data $sdir \
        --seqs $house/structures/seqs/$ogn.toalign \
        --structures $house/structures/$ogn/pdb \
        --outName $ogn \
        --dssp \
        --squeeze "H,E"\
        --createSubsets 30 --minSubsetID min \
        --outFolder $sdir \
        --tree true
done < $ogf