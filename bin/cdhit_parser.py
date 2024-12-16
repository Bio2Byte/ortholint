import pandas as pd
import sys
import math
import glob
#from pathlib import Path
from Bio import SeqIO


genome_occ=  int(sys.argv[1])
cdhit_seq= sys.argv[2]
cdhit_cluster= sys.argv[3]
treshold =float(sys.argv[4])
outgroup = sys.argv[5]


#cdhitfolder = sys.argv[2]
outdir =  '.' #'esmfold/'
clusters_files = [cdhit_cluster] #glob.glob(cdhitfolder+'*.clstr')


def convert_to_fasta(sequenceFiles,stype,genome_name):
    outfile=genome_name+'2l.'+stype
    count = SeqIO.convert(sequenceFiles, "fasta", outfile, "fasta-2line")
    print(outfile)
    print("Converted %i records" % count)
    return(outfile)
    
def cd_hit_to_df(clusters_file, seq_file):
    #read cdhit
    with open(clusters_file, 'r') as reps:
        match_dic = {}
        for line in reps:
            if line.startswith('>Cluster'):
                rep='Subset'+line.split(' ')[1].replace('\n','')
            else:
                id=line.split('>')[1].split('..')[0]
                match_dic[id] = rep
    cluster_df=pd.DataFrame.from_dict(match_dic, orient='index')

    cluster_df = cluster_df[0].rename('subset')
    cluster_df= cluster_df.reset_index()
    
    #read fasta
    og_name=seq_file.split('/')[-1].split('.')[0]
    ext=seq_file.split('.')[-1]
    seq_file_2l = convert_to_fasta(seq_file,ext,og_name)
    df = pd.read_csv(seq_file_2l, header=None, delimiter= '>')
    seq_df = pd.DataFrame({'index':df[1].iloc[::2].values, 'seq':df[0].iloc[1::2].values})
    df = pd.merge(seq_df, cluster_df, on='index')

    return(df)

def write_fasta_from_df(df,outname):
    out_name = outname + '.tofold.fasta'

    out_df = df[['index', 'seq']]
    out_df.iloc[:,0]= out_df.iloc[:,0] + '\n'

    rows = out_df.to_string(header=False,index=False,index_names=False).split('\n')

    with open (out_name, 'w') as m:
        for row in rows:
            row = row.replace('\\n','\n').replace(' ','')
            m.write( ">"+ outname+"_"+row + '\n')



for clusters_file in clusters_files:
    og = clusters_file.split('.')[0]
    og_name=og.split('/')[-1]
    print(og_name)
    
    seqs_file=og+'.fasta'
    df = cd_hit_to_df(clusters_file, seqs_file)

    subsets=df.subset.tolist()
    uniq_subset= list(set(subsets))

    noseqs=df.seq.count()

    if len(uniq_subset) == 1:
        outdf=df
        print("only 1 cluster")
    else:
        occurence=df[["seq", 'subset']].groupby('subset').count()
        ols =  occurence.index[occurence['seq']>= 0.05*noseqs].tolist()

        if len(uniq_subset) == len(ols):
            outdf=df
            print("keep all clusters")
        else:     
            #sel big clusters
            outdf=df[df.subset.isin(ols)]
            print ("remove singelton sequences:",df[~df['subset'].isin(ols)].seq.count())
    
    #check if they still contain at least 90% different 5 letter codes :D
    outdf['assembly']=outdf['index'].str.split('_').str[0]
    no_ass=outdf['assembly'].nunique()
    print(outdf)
    print(no_ass)

    if no_ass >= treshold*genome_occ:
        if outgroup in outdf['assembly'].tolist():
        #Path(outdir+og_name).mkdir(parents=True, exist_ok=True)
            write_fasta_from_df(outdf,og_name)
        else:
            print('missing outgroup')
    else:
        print('excluded')

