#!/bin/bash
#screen -S htbea bash -c ./bioenvada.sh

module load Nextflow/23.10.0

house=$VSC_SCRATCH_VO_USER/achromatium


#first: map id simsa -> nuc ! split ogtag_ ; redo nucmap or do externaly??
for data in OG0000554 OG0001136 OG0000614 ;do
    sdir=$house/clean_achromatium/structures/$data

    nextflow run $VSC_SCRATCH_VO_USER/bioenvada/bioenvada/main.nf -resume\
        -profile hydra \
        --type both \
        --targetSequences "$sdir/iqtree/${data}_simsa.fasta" \
        --nucToMap \
        --preprocessing "protein" \
        --alignSequences false\
        --buildTreeEvo \
        --efoldmine \
        --disomine \
        --fetchStructures false \
        --buildLogo \
        --cladePlots \
        --csubst \
        --branchIds 'all' \
        --foregroundFile $house/foreground.txt \
        --eteEvol 'M7,M8' \