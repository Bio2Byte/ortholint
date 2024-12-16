#!/bin/bash
#SBATCH --partition=ampere_gpu 
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=8
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=3G

#pascal_gpu has not enough mem for larger models
module load legacy-software
module load ESM-2/2.0.0-foss-2021a-CUDA-11.3.1
module load OpenFold/1.0.1-foss-2021a-CUDA-11.3.1


export TORCH_HOME=/databases/bio/ESM-2-2.0.0/torch

GHOME=$VSC_SCRATCH_VO/vsc10579/achromatium

cd $GHOME/clean_achromatium/structures/

seqstofold=$1
python3 $GHOME/achromatium_pipeline/bin/esmfold_inference.py --chunk-size 32  --max-tokens-per-batch 0 -i $seqstofold -o .
#python3 $GHOME/achromatium_pipeline/bin/extract_plddt.py pdb esm_fold_statistics.csv .

#for file in *.pae *.pdb; do
#    id=$(echo "$file" | cut -d'_' -f1)
#    extension="${file##*.}"
#    mkdir -p "$id/$extension"
#    mv "$file" "$id/$extension/"
#done