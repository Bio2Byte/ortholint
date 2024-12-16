
process mapAssemblies {
    publishDir "$params.outFolder/assembly_info", mode: "copy" 

    input:
    path genomes_dir

    output:
    path "all_assemblies_info.tsv", emit : assembly_summary
    path '*_info.tsv'

    script:
    """
    python3 $projectDir/bin/map_assemblies.py ${genomes_dir} . all_assemblies_info
    """
}

process runOrthofinder {
    publishDir "$params.outFolder", mode: "copy" 

    input:
    path genomes_dir

    output:
    path "orthofinder/Results_*/*"
    tuple path("orthofinder/Results_*/Orthogroups/Orthogroups.GeneCount.tsv"), path("orthofinder/Results_*/Orthogroup_Sequences") , emit : of_og_info_t

    script:
    """
    mkdir genomes_aseqs
    cp ${genomes_dir}/*/*.faa genomes_aseqs

    orthofinder -f genomes_aseqs -t ${task.cpus} -o orthofinder -os 

    rm -r orthofinder/Results_*/WorkingDirectory
    """
    //orthofinder -f genomes_aseqs -t ${task.cpus} -o ofres -os  -M msa
    //rm genomes_aseqs/*hypotheticals.faa
}

process runCDHIT { //per entry cdhit_seqs
    publishDir "$params.outFolder/cdhit", mode: "copy" 
    tag "${cdhit_seq.baseName}"
    errorStrategy 'ignore'

    input:
    //tuple val(og_name), path(cdhit_seq) 
    path cdhit_seq

    output:
    path "*.clstr" , emit: cdhit_cluster

    script:
    """
    $projectDir/bin/psi-cd-hit.pl -i ${cdhit_seq} -o ${cdhit_seq.baseName}  -c 0.2 
    """
}

process predictStructures {
    publishDir "${params.outFolder}", mode: "copy" 
    
    errorStrategy 'finish'
    //errorStrategy 'retry'
    //maxRetries 3

    input:
    path esm_seq

    output:
    path 'structures/*' 
    path  'structures/*/pdb', emit: struc_dir //cahnge to path  'structures/*', emit: struc_dir
         
    script:
    """
    export TORCH_HOME=/databases/bio/ESM-2-2.0.0/torch
    
    python3 $projectDir/bin/esmfold_inference.py --chunk-size 32  --max-tokens-per-batch 0 -i ${esm_seq} -o .

    mkdir -p structures
    for file in *.pae *.pdb; do
        id=\$(echo "\$file" | cut -d'_' -f1)
        extension="\${file##*.}"
        mkdir -p "structures/\$id/\$extension"
        mv "\$file" "structures/\$id/\$extension/"
    done
    """
}
//python3 $projectDir/bin/extract_plddt.py ${struc_dir}/pdb esm_fold_statistics.tsv .


process runFoldseek {//per og_strucs folder
    publishDir "$params.outFolder/structures/${og_name}/foldseek", mode: "copy" 

    input:
    tuple val(og_name), path (struc_dir)
    val gate

    output:
    path "*.tsv" 
    tuple val(og_name), path("*general_structure_similarity.tsv") , emit: struc_sims
    
    script:
    //struc_dir folder can only contain pdb files!
    """
    header="query	target	seq_len_q	seq_len_t	aln_len	alntmscore	lddt	rmsd	prob_homolog	no_mismatch	coverage_q	coverage_t	e_value	perc_seq_id	no_res_id"
    
    rm *similarity.tsv
    mac=/scratch/brussel/vo/000/bvo00023/vsc10579/.apptainer

    for exa in ${struc_dir}/*.pdb ; do

        exn=\$(\$exa .pdb | sed 's/\..*\$//')

        easy-search \$exa ${struc_dir} \${exn}_strucsim.tsv tmp --max-seqs 2000 --exhaustive-search true --format-output "query,target,qlen,tlen,alnlen,alntmscore,lddt,rmsd,prob,mismatch,qcov,tcov,evalue,pident,nident"
        
        cat \${exn}_strucsim.tsv >> ${og_name}_general_structure_similarity.tsv
        sed -i "1s/.*/\$header/" \${exn}_strucsim.tsv
        rm -r tmp
        echo '-----------------------------------------------------------'
        echo '-----------------------------------------------------------' 
        echo '-----------------------------------------------------------' 
    done

    sed -i "1s/.*/\$header/" ${og_name}_general_structure_similarity.tsv
    echo 'completed'
    """
}//apptainer run --bind \$PWD:/data  \$mac/foldseek.sif easy-search \$exa ${struc_dir} \${exn}_strucsim.tsv tmp --max-seqs 2000 --exhaustive-search true --format-output "query,target,qlen,tlen,alnlen,alntmscore,lddt,rmsd,prob,mismatch,qcov,tcov,evalue,pident,nident"
        
        
