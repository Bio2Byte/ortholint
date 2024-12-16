#!/bin/bash
#SBATCH --job-name=foldseek
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G


APPTAINERCACHE=$VSC_SCRATCH_VO_USER/.apptainer
header="query	target	seq_len_q	seq_len_t	aln_len	alntmscore	lddt	rmsd	prob_homolog	no_mismatch	coverage_q	coverage_t	e_value	perc_seq_id	no_res_id"

#apptainer build $APPTAINERCACHE/foldseek.sif docker://mcpalumbo/foldseek:1
#folder can only contain pdb files!

GHOME=$1
rm $GHOME/general_structure_similarity.tsv
rm -r $GHOME/tmp

mkdir -p $GHOME/foldseek
echo $PATH
for exa in $GHOME/pdb/*.pdb ;do
    echo $exa
    data=$(basename $exa .pdb | sed 's/\..*$//')
    echo $data
    apptainer run --bind $GHOME:/data  $APPTAINERCACHE/foldseek.sif easy-search $exa $GHOME/pdb $GHOME/foldseek/${data}_strucsim.tsv $GHOME/tmp --max-seqs 2000 --exhaustive-search true --format-output "query,target,qlen,tlen,alnlen,alntmscore,lddt,rmsd,prob,mismatch,qcov,tcov,evalue,pident,nident"
    
    cat $GHOME/foldseek/${data}_strucsim.tsv >> $GHOME/general_structure_similarity.tsv
    sed -i "1s/.*/$header/" $GHOME/foldseek/${data}_strucsim.tsv
    rm -r $GHOME/tmp
    echo '-----------------------------------------------------------'
    echo '-----------------------------------------------------------'
    echo '-----------------------------------------------------------' 
    echo '-----------------------------------------------------------'
    echo '-----------------------------------------------------------' 
    echo '-----------------------------------------------------------' 
done

sed -i "1s/.*/$header/" $GHOME/general_structure_similarity.tsv
echo 'completed'