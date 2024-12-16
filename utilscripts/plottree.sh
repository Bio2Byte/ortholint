#!/bin/bash
#SBATCH --job-name=treeplot
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1G


APPTAINERCACHE=$VSC_SCRATCH_VO_USER/.apptainer

GHOME=$VSC_SCRATCH_VO_USER/achromatium/full_achromatium/og_clean/bea

for tree in $GHOME/OG*/tree/*_rooted.treefile ; do
    echo $tree
    apptainer run --bind $VSC_SCRATCH_VO_USER:/data $APPTAINERCACHE/slheidig-ad_ete3.img python $VSC_SCRATCH_VO_USER/bioenvada/bioenvada/bin/treePlot.py $tree 1
done