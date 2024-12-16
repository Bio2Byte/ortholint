#this script takes as input folders of annotated bins/assemblies/genomes and map aa/nuc seqs per gene and provide general information

from pathlib import Path
import glob
import sys
import pandas as pd
from Bio import SeqIO

genomes_folder=sys.argv[1]
genomes=glob.glob(genomes_folder+'/*')
#genomes=glob.glob(genomes_folder)
outdir=sys.argv[2] #"genomes_info"
#Path(outdir).mkdir(parents=True, exist_ok=True)

outname=sys.argv[3]


def convert_to_fasta(sequenceFiles,stype,genome_name):
    outfile=genome_name+'2l.'+stype
    count = SeqIO.convert(sequenceFiles, "fasta", outfile, "fasta-2line")
    print(outfile)
    print("Converted %i records" % count)
    return(outfile)

def fasta_to_df(seq_file,stype,genome_name):

    seq_file_2l = convert_to_fasta(seq_file,stype,genome_name)
    df = pd.read_csv(seq_file_2l, header=None, delimiter= '>')
    seq_df = pd.DataFrame({'Locus Tag':df[1].iloc[::2].values, 'seq_'+stype:df[0].iloc[1::2].values})
    seq_df['Locus Tag']= seq_df['Locus Tag'].str.lstrip(' ')
    seq_df['Locus Tag']= seq_df['Locus Tag'].str.lstrip('>')
    seq_df['Locus Tag']= seq_df['Locus Tag'].str.split(' ').str[0]
    
    print(seq_df)
    return (seq_df)


def whole_genome_info(tsv_df):

    try:
        no_contigs= tsv_df['#Sequence Id'].nunique()
    except:
        tsv_df['contigs']=tsv_df['Locus Tag'].str.split('_').str[2]
        no_contigs= tsv_df['contigs'].nunique()
    no_genes,cols=tsv_df.shape
    locustag=tsv_df.loc[0,'Locus Tag'].split('_')[0]

    return ([no_contigs,no_genes,locustag])


wholegenomeinfo_dic={}
for genome in genomes:

    genome_name = genome.split('/')[-1]
    print(genome_name)

    try:
        fna=glob.glob(genome+'/'+genome_name+'*.ffn')[0]
        print(fna)
    except:
        print("!!!!!!!!!!!!!!!!!no genome here!", genome_name)
        continue
    fna_df = fasta_to_df(fna,'nuc',genome_name)

    faa=glob.glob(genome+'/'+genome_name+'*.faa')
    for f in faa:
        if "hypothetical" not in f:
            faa_df = fasta_to_df(f,"aa",genome_name)
    
    try:
        tsv=glob.glob(genome+'/*.tsv')
        for t in tsv:
            if "hypothetical" not in t:
                tsv_df=pd.read_csv(t, skiprows=5 , sep='\t')
        
        tsv_df=tsv_df.merge(faa_df, on='Locus Tag')
        tsv_df=tsv_df.merge(fna_df, on='Locus Tag')
        tsv_df.to_csv(outdir+'/'+genome_name+'_info.tsv', sep='\t')
    except:
        print(fna_df)
        print(faa_df)
        tsv_df=faa_df.merge(fna_df, on='Locus Tag')
    print(tsv_df)
    wholegenomeinfo_dic[genome_name]= whole_genome_info(tsv_df)
    faa_df=fna_df=tsv_df=0

wholegenomeinfo_df=pd.DataFrame.from_dict(wholegenomeinfo_dic, orient='index',columns=['no_contigs','no_genes','anno_tag'])
wholegenomeinfo_df.sort_index(inplace=True)
wholegenomeinfo_df.to_csv(outdir+'/'+outname+'.tsv', sep='\t')





            