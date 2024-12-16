# ortholint
A nextflow pipeline that creates (or takes) orthologous groups from gene-called genomes and runs sequence- and structure based QC to remove outliers

## vision

![Vision of future ortholint workflow!](ortholint_vision.png "Vision of future ortholint workflow")

## setup
You have to define your own dependency profile in *nextflow.config* unless you have access to the VUB-HPC. Use docker/apptainer to make your life easier. 

## parameters
--data "$launchDir/data" \
--outFolder "$launchDir/" \
--structures false [path/to/structure/dir] \
--predictRemaining false \
--orthoGroupSeqs  [path/to/ogseqs/dir] \
--minCoverage 0.8 \
--outGroup  false [str_with_og_name] \
--foldseekClusters false 

