#!/bin/bash
#SBATCH --job-name=processsfoldseek
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

module load SciPy-bundle/2022.05-foss-2022a
module load matplotlib/3.5.2-foss-2022a

GHOME=$VSC_SCRATCH_VO/vsc10579/achromatium


python3 $GHOME/bin/process_foldseek.py $GHOME/esmfold $GHOME/ogs_for_simsapiper

<<Block
for folder in $GHOME/esmfold/* ; do
    echo $folder
    cd $folder
    for seqf in $folder/*.fasta ; do

        noseqs=$(grep -c ">" $seqf)
        if (( noseqs > 200 )); then
            echo "to the resubmit list!"
            echo $noseqs $seqf >> $GHOME/esmfold/ogs_to_resubmit.txt
        fi

        sbatch $GHOME/bin/esmfold.sh $seqf
    done
    cd ../
done
Block

<<Block
GHOME=$VSC_SCRATCH_VO/vsc10579/achromatium
for folder in ${GHOME}/esmfold/* ; do
    echo $folder
    cd $folder
    sbatch $GHOME/bin/foldseek.sh $folder
    cd ../
done
Block