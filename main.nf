
log.info """\
module load Nextflow/23.10.0
nextflow run achromatium_pipeline/main.nf -resume -profile hydra --data nftester/assemblies --outFolder nftester
================================================================================
                                LIST OF PARAMETERS
================================================================================
                                GENERAL

Launch dir      : $launchDir
Project dir     : $projectDir
Execution time  : $params.executionTimestamp
================================================================================
                                INPUT FILES

Input file folder (--data): $params.data
================================================================================
                                OUTPUT FILES

Output folder (--outFolder): $params.outFolder
================================================================================
================================================================================
"""

include {
    mapAssemblies;
    runOrthofinder;
    runCDHIT;
    predictStructures;
    //runFoldseek;
} from "$projectDir/modules/io"

include {
    analyseOrthofinder;
    analyseCDHIT;
    analyseFoldseek;
} from "$projectDir/modules/data_cleaning"




workflow {

genomes_ch= Channel.fromPath(params.data)
mapped_genomes_ch=mapAssemblies(genomes_ch)

genomes_aa_ch =Channel.fromPath("${params.data}/*/*.faa")
//genomes_aa_ch.view()
no_assemblies=genomes_aa_ch.count()
no_assemblies.view{'No. of assemblies: ' + it}
//problem:hypothetical file

if (params.orthoGroupSeqs){
    //cdhit_seqs=Channel.fromPath("${params.orthoGroupSeqs}/*.fasta")
    geneCount_ch=Channel.fromPath("${params.orthoGroupSeqs}/Orthogroups/Orthogroups.GeneCount.tsv")
    ogseqs_ch=Channel.fromPath("${params.orthoGroupSeqs}/Orthogroup_Sequences")
    
    of=geneCount_ch.combine(ogseqs_ch).combine(no_assemblies)
}else{
    runOrthofinder(genomes_ch)
        of_out_t=runOrthofinder.out.of_og_info_t
    of=of_out_t.combine(no_assemblies)
}

analyseOrthofinder(of)
    cdhit_seqs=analyseOrthofinder.out.cd_seqs.flatten()

runCDHIT(cdhit_seqs)
    cd_cluster=runCDHIT.out.cdhit_cluster

cdhit_seqs_t=cdhit_seqs.map(it -> tuple( it.baseName , it))
cd_cluster_t=cd_cluster.map(it -> tuple( it.baseName , it))


cdhit_seqs_t.join(cd_cluster_t, remainder:true)
    .branch {
        seqNotFound:        it[1] == null
        clusterNotFound:       it[2] == null
        mapped:             true
    }.set{matchCDHIT}
    matchCDHIT.mapped.count().view{"mapped for cdhit " + it}

    matchCDHIT.mapped.combine(no_assemblies) | analyseCDHIT
    seqs_to_model=analyseCDHIT.out.esm_seq
        

all_seqs_to_model = seqs_to_model
        .splitFasta( record: [header: true,  sequence: true ])
        .map{tuple(it.header, it.sequence)}
all_seqs_to_model.count().view{"total number of sequences to model "+it}
    
if (params.predictRemaining){
    foundStrucs=Channel.fromPath("$params.structures/*/pdb/*.pdb").map(it -> tuple(it.baseName, it))
    foundStrucs.count().view{"total number of found structures "+it}
    //foundStrucs.buffer(size: 5).first().view()
    
    all_seqs_to_model.join(foundStrucs, remainder:true)
        .branch {
            seqNotFound:        it[1] == null
            strucNotFound:       it[2] == null
            mapped:             true
        }.set{matchSeqStru}

        matchSeqStru.mapped.count().view{"mapped seqs to structures " + it}
        matchSeqStru.strucNotFound.count().view{"structures to model " + it}
    onefile_seqs_to_model=matchSeqStru.strucNotFound
} else {
    onefile_seqs_to_model=all_seqs_to_model
}


onefile_seqs_to_model.buffer(size: 10).first().view()
onefile_seqs_to_model
    .collectFile( name: "seqs_to_model.fasta" ,storeDir:"${params.outFolder}") {
            item -> '>' + item[0] + '\n' + item[1]}
    .splitFasta(by: 5000 , file: true )
    .set{batch_seqs_to_model_ch}
batch_seqs_to_model_ch.view()

predictStructures(batch_seqs_to_model_ch)
    struc_dirs=predictStructures.out.struc_dir
    gate=struc_dirs.count()

if (params.structures){
    struc_dirs_ex=Channel.fromPath("${params.structures}/OG*/pdb", type: 'dir' )
}else{
    struc_dirs_ex=Channel.empty()
}
struc_dirs_ex.count().view{"mapped struc dirs  " + it}
//struc_dirs_ex.buffer(size: 5).first().view()

struc_dir_t=struc_dirs.collect()
    .mix(struc_dirs_ex)
    .flatten()
    .unique()
struc_dir_t.count().view{"struc dirs for foldseek " + it}

if (params.foldseekClusters){
    sims=Channel.fromPath("${params.foldseekClusters}/OG*/general_structure_similarity.tsv")
    sims.buffer(size: 3).first().view()
    sims.count().view{"found foldseek " + it}
    sims_t=sims.map(it -> tuple(it.parent.baseName ,it)) //give it ogid
}else{

    //struc_dir_t.map(it -> tuple(it.parent.baseName ,it)).combine(gate) | runFoldseek
    //    sims_t=runFoldseek.out.struc_sims
}

seqs_to_model_t=seqs_to_model.map(it -> tuple(it.baseName[0..8], it))
seqs_to_model_t.join(sims_t, remainder:true)
    .branch {
        seqNotFound:        it[1] == null
        strucSimNotFound:       it[2] == null
        mapped:             true
    }.set{matchSeqStruSim}
    //matchSeqStruSim.mapped.count().view{"mapped for foldseek analysis " + it}

matchSeqStruSim.mapped.combine(no_assemblies) | analyseFoldseek
    seqs_to_align=analyseFoldseek.out.foldseek_seq
    seqs_to_align.count().view{"Final number of OGs for alignment " + it}
}

workflow.onComplete {
    println "Pipeline completed at               : $workflow.complete"
    println "Time to complete workflow execution : $workflow.duration"
    println "Execution status                    : ${workflow.success ? 'Success' : 'Failed' }"
    println "Output folder                       : $params.outFolder"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: \n ${workflow.errorReport}"    
}

