import pandas as pd
import sys
import math
import glob
from pathlib import Path
import matplotlib.pyplot as plt 

pd.options.mode.chained_assignment = None #hide warning
#module load SciPy-bundle/2022.05-foss-2022a
#module load matplotlib/3.5.2-foss-2022a




genome_occ= int(sys.argv[1])# 
outdir = '.' 
seqs_file =sys.argv[2]
strucfolder =  sys.argv[3] #'esmfold'
og_name=sys.argv[4]

treshold =float(sys.argv[5])
outgroup = sys.argv[6]

def write_fasta_from_df(df,outname):
    out_name = outname + '.toalign.fasta'

    out_df = df[['query', 'seq']]
    out_df.iloc[:,0]= out_df.iloc[:,0] + '\n'

    rows = out_df.to_string(header=False,index=False,index_names=False).split('\n')

    with open (out_name, 'w') as m:
        for row in rows:
            row = row.replace('\\n','\n').replace(' ','')
            m.write( ">"+row + '\n')

def fasta_to_df(seq_file):
    df = pd.read_csv(seq_file, header=None)
    seq_df = pd.DataFrame({'query':df[0].iloc[::2].values, 'seq':df[0].iloc[1::2].values})
    seq_df['query']= seq_df['query'].str.lstrip('>')

    seq_df['seq']= seq_df['seq'].str.upper()
    #seq_df['seq']= seq_df['seq'].str.replace('*','X').str.replace('?','X').str.replace('!','X').str.replace('+','X')

    return (seq_df)

def pairwise_to_distance_matrix(df,og_name,prop):
    # make unique, sorted, common index
    idx = sorted(set(df['query']).union(df['target']))
    # reshape
    dm_df=(df.pivot(index='query', columns='target', values=prop)
    .reindex(index=idx, columns=idx)
    .fillna(0, downcast='infer')
    .pipe(lambda x: x+x.values.T)
    )

    outpath=strucfolder+"/foldseek/"+og_name+'distance_matrix_'+prop
    Path(strucfolder+"/foldseek").mkdir(parents=True, exist_ok=True)
    #outpath=og_name+'distance_matrix_'+prop
    dm_df.to_csv(outpath+'.tsv', sep='\t')

    plt.imshow(dm_df, cmap='viridis', interpolation='nearest', vmin=dm_df.min().min(), vmax=dm_df.max().max())

    plt.title("Distance matrix of " + og_name +" structures")
    plt.colorbar(label=prop.upper())


    plt.gca().invert_yaxis()
    num_ticks = len(dm_df.columns)
    fontsize = max(2, 12 - int(num_ticks / 10)) 
    plt.xticks(ticks=range(len(dm_df.columns)), labels=dm_df.columns, fontsize=fontsize, rotation=90)
    plt.yticks(ticks=range(len(dm_df.index)), labels=dm_df.index, fontsize=fontsize)
    
    plt.savefig(outpath+'.pdf',bbox_inches='tight')
    plt.close()


seqs = fasta_to_df(seqs_file)
strucclus = pd.read_csv(strucfolder+"/general_structure_similarity.tsv",sep='\t')

try:
    pairwise_to_distance_matrix(strucclus[['query','target','lddt']],og_name,'lddt')#or rmsd lddt
    pairwise_to_distance_matrix(strucclus[['query','target','rmsd']],og_name,'rmsd')
except:
    print('missing foldseek values, rerun:', og_name)

#strucclus[['prob_homolog']]=strucclus[['prob_homolog']].apply(pd.to_numeric)

if strucclus['query'].nunique() != seqs['query'].nunique():
    print ('warning, number of seqs dont match')
sel_strucs=strucclus[strucclus['prob_homolog']>=0.5]

sel_strucs_l=sel_strucs['query'].unique().tolist()
outdf=seqs[seqs['query'].isin(sel_strucs_l)]

#check if they still contain at least 90% different 5 letter codes :D
outdf['assembly']=outdf['query'].str.split('_').str[1]
no_ass=outdf['assembly'].nunique()

print(outdf)
print(no_ass)

if no_ass >= treshold*genome_occ:
    if outgroup in outdf['assembly'].tolist():
        mkoutdir=outdir+'/seqs'
        Path(mkoutdir).mkdir(parents=True, exist_ok=True)
        write_fasta_from_df(outdf,mkoutdir+'/'+og_name)
    else:
        print('missing outgroup')
else:
    print('excluded')