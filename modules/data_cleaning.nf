process analyseOrthofinder {
    publishDir "$params.outFolder/", mode: "copy" 
    publishDir "$params.outFolder/", mode: "copy" , pattern: "cdhit/*.fasta"
    //publishDir "$params.outFolder/ogs/", mode: "copy" , pattern: "*.fasta"

    input:
    tuple path(of_ogs_info), path(of_og_seqs), val(genome_occ)

    output:
    //path "*seq_len.png" //per og plot of seq selection based seq len dist
    //path "*og.png"  //general seq len dist pre filtering 
    //path "*info.tsv"  //general seq len dist pre filtering 

    path "ogs_info/" 
    path "cdhit/*.fasta" , emit : cd_seqs //fastas for cdhit  //

    //TODO:
    //path "*.fasta" //excluded bc too many unchar AA (X) or seq len
    //upper/lower boundary parameters

    script:
    """
    python3 $projectDir/bin/seq_qc_len.py ${of_ogs_info} ${of_og_seqs} ogs_info cdhit 100 1000 ${genome_occ} ${params.minCoverage} ${params.outGroup}
    """
    // python3 $projectDir/bin/seq_qc_len.py ${of_ogs_info} ${of_og_seqs} .  . "total min len" "total max len" 
}

process analyseCDHIT {//per entry cdhit_seqs
    //publishDir "$params.outFolder/ogs/${og_name}/esmfold", mode: "copy" 
    publishDir "$params.outFolder/structures/seqs_to_fold", mode: "copy" 
    tag "${og_name}"

    input:
    
    tuple val(og_name), path(cdhit_seq), path(cdhit_cluster), val(genome_occ)

    output:
    path "*.tofold.fasta" , emit: esm_seq , optional: true

    //TODO:
    //path "*.fasta" //excluded bc not enough coverage of input assemblies

    script:
    """
    python3 $projectDir/bin/cdhit_parser.py ${genome_occ} ${cdhit_seq} ${cdhit_cluster} ${params.minCoverage} ${params.outGroup}
    """
}

process analyseFoldseek {//per esm_seq
    publishDir "$params.outFolder/structures/${og_name}/foldseek_summary", mode: "copy"
    publishDir "$params.outFolder/simsa/${og_name}/seq", mode: "copy" , pattern: "*.fasta"

    tag "${esm_seq.baseName}"

    input:
    tuple val(og_name),path(esm_seq), path(struc_sims), val(genome_occ)

    output:
    path "*.toalign.fasta" , emit: foldseek_seq , optional: true
    path "*.tsv"
    path "*.png"

    //TODO:
    //path "*.fasta"//excluded bc not enough coverage of input assemblies

    script:
    """
    python3 $projectDir/bin/process_foldseek.py ${genome_occ} ${esm_seq} ${struc_sims} ${og_name} ${params.minCoverage} ${params.outGroup}
    """
}
