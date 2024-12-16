#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2G


APPTAINERCACHE=$VSC_SCRATCH_VO_USER/.apptainer
#apptainer build $APPTAINERCACHE/orthofinder.sif docker://davidemms/orthofinder:3.0.1b1

GHOME=$VSC_SCRATCH_VO/vsc10579/achromatium


#apptainer run --bind $GHOME:/data  $APPTAINERCACHE/orthofinder.sif orthofinder -h 
apptainer run --bind $GHOME:/data  $APPTAINERCACHE/orthofinder.sif orthofinder -f $VSC_SCRATCH_VO_USER/achromatium/only_genomes_cellbin -o $VSC_SCRATCH_VO_USER/achromatium/orthologs/sel_of2