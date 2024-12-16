#!/bin/bash
#SBATCH --job-name=cdhit
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

#run 20% cdhit

module load CD-HIT/4.8.1-GCC-11.3.0
module load BLAST+/2.13.0-gompi-2022a
module load SciPy-bundle/2022.05-foss-2022a


GHOME=$VSC_SCRATCH_VO/vsc10579/achromatium

seqsf=$1


data=$(basename $seqsf .fasta | sed 's/\..*$//')
echo $data

WD=$GHOME/nftester/cdhit/cdhit/work/
mkdir -p $WD
cp $seqsf $WD
cd $WD

$GHOME/achromatium_pipeline/bin/psi-cd-hit.pl -i ${data}.fasta -o  ${data}  -c 0.2 

cp ${data}.clstr $GHOME/nftester/cdhit/cdhit