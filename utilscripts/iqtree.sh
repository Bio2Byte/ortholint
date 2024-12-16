#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --time=1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G

module purge 
module load IQ-TREE/2.2.2.7-gompi-2023a

house=$VSC_SCRATCH_VO_USER/achromatium

for data in OG0000554 ;do
    sdir=$house/clean_achromatium/structures/$data
    mkdir -p $sdir/iqtree
    cd $sdir/iqtree

    #cp $sdir/simsa/msas/squeezed_merged_${data}.fasta ${data}_simsa.fasta

    iqtree2 -s ${data}_simsa.fasta -nt AUTO -B 10000

done



#in each iqtree logfile: 
#Total tree length (sum of branch lengths): 14.7697
#Sum of internal branch lengths: 6.3421 (42.9400% of tree length)
#if % of tree length < 50 OR ttl >10: flag for long branch?

#else:
#simsa ali mapped onto nuc seqs from table
#launch bioenvada?