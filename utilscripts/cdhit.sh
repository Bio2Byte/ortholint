#!/bin/bash
#SBATCH --job-name=cdhit
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

#run 20% cdhit

module load CD-HIT/4.8.1-GCC-11.3.0
module load BLAST+/2.13.0-gompi-2022a
module load SciPy-bundle/2022.05-foss-2022a


GHOME=$VSC_SCRATCH_VO/vsc10579/achromatium
mkdir $GHOME/cdhit/work
for seqsf in $GHOME/cdhit/*.fasta ; do
    sbatch $GHOME/bin/batch_cdhit.sh $seqsf
done

<<Block 
    data=$(basename $seqsf .fasta | sed 's/\..*$//')
    echo $data
    cp $seqsf $GHOME/cdhit/work/
    cd $GHOME/cdhit/work

    $GHOME/bin/psi-cd-hit.pl -i ${data}.fasta -o  ${data}  -c 0.2 

    cp ${data}.clstr $GHOME/cdhit/
    cd ../
done
Block
    

rm -r $GHOME/cdhit/work/
python3 $GHOME/bin/cdhit_parser.py 18 $GHOME/cdhit/

