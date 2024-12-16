#module load matplotlib/3.5.2-foss-2022a
#module load SciPy-bundle/2022.05-foss-2022a

import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt 
from pathlib import Path
import os
from Bio import SeqIO

#this script takes folder of proposed orthologs 
og_candi_tab        = sys.argv[1] #'orthologs/sel_of/Results_Nov19/Orthogroups/Orthogroups.GeneCount.tsv'  #sys.argv[1] #
og_candi_folder     = sys.argv[2] # 'orthologs/sel_of/Results_Nov19/Orthogroup_Sequences' #sys.argv[2] #

outdir              = sys.argv[3] #"orthogroups_info"
Path(outdir).mkdir(parents=True, exist_ok=True)

outdircdhit         = sys.argv[4] #"cdhit"
Path(outdircdhit).mkdir(parents=True, exist_ok=True)

min_boundary        = int(sys.argv[5]) #100 True
max_boundary        = int(sys.argv[6]) #1000 True #set false to not remove extra long seqs
genome_occ          = int(sys.argv[7])
treshold            =float(sys.argv[8])
minocc              =treshold*genome_occ
outgroup            = sys.argv[9]

def parse_og_candi_tab(og_candi_tab):
    og_candi_df = pd.read_csv(og_candi_tab, sep='\t')
    
    # Generate histogram plot for total occurrence of sequences in orthogroups
    binsize=max(1,round(og_candi_df["Total"].max() * 0.05))
    plt.hist(og_candi_df["Total"], bins=range(0, og_candi_df["Total"].max()+binsize, binsize), edgecolor='black')
    plt.yscale('log')
    plt.xlabel("Number of sequences in orthogroup")
    plt.ylabel("Frequency")
    plt.title("Distribution of total occurrence of sequences in orthogroups")
    plt.savefig(outdir+'/total_occurrence_og.png')
    plt.close()
    
    og_candi_df['unique_occurrence'] = genome_occ-(og_candi_df == 0).astype(int).sum(axis=1)

    # Generate histogram plot for unique occurrence of sequences in orthogroups
    
    binsize=max(1,round(og_candi_df["unique_occurrence"].max() * 0.05))
    plt.hist(og_candi_df["unique_occurrence"], bins=range(0,og_candi_df["unique_occurrence"].max()+binsize, binsize), edgecolor='black')
    plt.xlabel("Number of assemblies in orthogroup")
    plt.ylabel("Frequency")
    plt.title("Distribution of unique occurrence of assemblies in orthogroups")
    plt.savefig(outdir+'/unique_occurrence_og.png')
    plt.close()

    # filter can only be 0 less then 90% of all cols;
    og_candi_df=og_candi_df[og_candi_df['unique_occurrence']>minocc]
    og_occ_sel=og_candi_df['Orthogroup'].tolist()
    return (og_occ_sel)




def convert_to_fasta(sequenceFiles,stype,genome_name):
    outfile=genome_name+'2l.'+stype
    count = SeqIO.convert(sequenceFiles, "fasta", outfile, "fasta-2line")
    print(outfile)
    print("Converted %i records" % count)
    return(outfile)

def fasta_to_df(seq_file):

    og_name=seq_file.split('/')[-1].split('.')[0]
    ext=seq_file.split('.')[-1]
    seq_file_2l = convert_to_fasta(seq_file,ext,og_name)
    
    df = pd.read_csv(seq_file_2l, header=None, delimiter= '>')
    
    seq_df = pd.DataFrame({'Locus Tag':df[1].iloc[::2].values, 'seq':df[0].iloc[1::2].values})
    seq_df['Locus Tag']= seq_df['Locus Tag'].str.lstrip(' ')
    seq_df['Locus Tag']= seq_df['Locus Tag'].str.lstrip('>')

    seq_df['seq']= seq_df['seq'].str.upper()
    #seq_df['seq']= seq_df['seq'].str.replace('*','X').str.replace('?','X').str.replace('!','X').str.replace('+','X')

    return (seq_df)

def write_fasta_from_df(df,labelcol,seqcol,outname):
    out_name = outname + '.fasta'
    out_df = df.filter(items=[labelcol, seqcol])
    out_df[labelcol] =  out_df[labelcol]+"\n"
    rows = out_df.to_string(header=False,index=False,index_names=False).split('\n')
    with open (out_name, 'w') as m:
        for row in rows:
            row=row.replace('*','').replace('#','').replace('+','')
            row=row.replace('\\n','\n').replace(' ',"")
            m.write( '>'+row + '\n')

def seq_len_boundaries(seq_df,og_name):
    #self? min_boundary,max_boundary

    seq_df.sort_values(by=['seq_len'], inplace=True, ascending=False)

    stats_df=seq_df['seq_len'].describe()
    stats_df['mean_median_skew']=stats_df['mean']-stats_df['50%']

    sl_min=best_min_l=stats_df['min']
    sl_max=best_max_l=stats_df['max']
    stats_df['full_spread']=sl_max-sl_min
    
    stats_df['quartile_count']=stats_df['count']/2
    stats_df['quartile_spread']=stats_df['75%']-stats_df['25%'] 
    if stats_df['full_spread'] >0 :
        evaluation_q=(stats_df['count']*stats_df['quartile_spread'])/(stats_df['count']*0.5*stats_df['full_spread'])
    else:
        evaluation_q =0

    stats_df['quartile_factor']=evaluation_q
    final_score=evaluation_q

    sl_SD=stats_df['std']
    sl_mean=stats_df['mean']

    seq_lens=seq_df['seq_len'].tolist()
    for factor in [1,2,3,'0.5x2']:

        if factor =='0.5x2':
            min_lim=sl_mean/2
            max_lim=sl_mean*2
        else:
            min_lim=sl_mean-(factor*sl_SD)
            max_lim=sl_mean+(factor*sl_SD)

        #sequ len between 100-1000 aa
        if min_lim <=sl_min:
            min_lim=sl_min
        if max_lim >=sl_max:
            max_lim=sl_max 
        
        #do not set a min or max len boundary; min max len boundary are parameters now
        #if not min_boundary:
        #    min_lim = sl_min
        #if not max_boundary:
        #    max_lim = sl_max
        
        sd_len_df=seq_df['seq_len'][seq_df['seq_len'].between(min_lim,max_lim,"both")]

        stats_sd_len_df=sd_len_df.describe()
        stats_sd_len_df['full_spread']=stats_sd_len_df['max']-stats_sd_len_df['min']
        sd_count=stats_sd_len_df['count']
        sd_spred=stats_sd_len_df['full_spread']

        if stats_df['full_spread'] >0:
            evaluation=(stats_df['count']*sd_spred)/(sd_count*stats_df['full_spread'])
        else:
            evaluation =0

        stats_df['sd'+str(factor)+'_count']=sd_count
        stats_df['sd'+str(factor)+'_spread']=sd_spred
        stats_df['sd'+str(factor)+'_factor']=evaluation

        if evaluation_q > evaluation:
            best_min_l =min_lim
            best_max_l =max_lim
            final_score = evaluation_q
        else:
            evaluation_q = evaluation

    #ideal seqlen boundary: no_seqs max and seqs_spread min. 
    #hence: max of this should be ideal? (no_full_seqs*sel_spread)/(no_sel_seqs*full_spread) > is max or min better?
    #problem: no reason for this distribution to be on bell curve
    #calc no of seqs and spread if seq_len_boundaries=1,2,3SD
    #if not max_boundary: only calc for lower_boundary

    #alternative: minlim =2/mean, maxlim 2xmean



    # add no_seqs and spread for each seqlenboundary to stats_df

    #ideal seqlen boundary: no_seqs max and seqs_spread min. 
    #hence: max of this should be ideal? (no_full_seqs*sel_spread)/(no_sel_seqs*full_spread)
    #filter df for ideal boundary?
    
    bounded_seq_df=seq_df[seq_df['seq_len'].between(best_min_l,best_max_l,"both")]   
    bounded_seq_df['best_quart_score']=final_score

    if bounded_seq_df['assembly_id'].count()>minocc:
        if outgroup in bounded_seq_df['assembly_id'].tolist():
            write_fasta_from_df(bounded_seq_df,'Locus Tag','seq',outdircdhit+'/'+og_name)
            # Generate histogram plot for sequences length in orthogroups
            binsize=max(1,round(seq_df["seq_len"].max() * 0.05))
            plt.hist(seq_df["seq_len"], bins=range(seq_df["seq_len"].min(),seq_df["seq_len"].max()+binsize, binsize), edgecolor='black' , color='gray')
            plt.hist(bounded_seq_df["seq_len"], bins=range(seq_df["seq_len"].min(),seq_df["seq_len"].max()+binsize, binsize), color='skyblue', edgecolor='white',alpha=0.3)
            plt.xlabel("Length of sequence in orthogroup")
            plt.ylabel("Frequency")
            plt.title("Distribution of sequence length in orthogroups")
            if stats_df['count'] >100:
                plt.yscale('log')
            plt.savefig(outdir+'/'+og_name+'_seq_len.png')
            plt.close()

            stats_df['status']='sel'
        else:
            stats_df['status']='ex_no-outgroup'
    else:
        stats_df['status']='ex_unique-occ'
    print(stats_df.index.tolist())
    cols=stats_df.index.tolist()

    return (stats_df.tolist(),cols)




og_occ_sels = parse_og_candi_tab(og_candi_tab)
all_stats_dic={}
for og_occ_sel in og_occ_sels:
    print(og_candi_folder + '/' +og_occ_sel + ".fa")
    seq_df = fasta_to_df(og_candi_folder + '/' +og_occ_sel + ".fa")

    seq_df['seq_len']=seq_df.seq.str.len()

    #remove all seqs with more then 5% unresolved aminoacids
    seq_df['count_X'] = seq_df.seq.str.count('X')
    seq_df=seq_df[seq_df['count_X']<seq_df['seq_len']*0.05]

    #find ideal remove extremely different lengthed seqs
    #100 AA cutoff

    print("all seqs",seq_df.shape)
    seq_df=seq_df[seq_df['seq_len']>=min_boundary]
    print('rm seqs shorter then AA',min_boundary,seq_df.shape)
    seq_df=seq_df[seq_df['seq_len']<=max_boundary]
    print('rm seqs longer then then AA',max_boundary,seq_df.shape)

    seq_df['assembly_id']=seq_df['Locus Tag'].str.split("_").str[0]

    boundedlist, cols =seq_len_boundaries(seq_df,og_occ_sel) 
    all_stats_dic[og_occ_sel]=boundedlist 

all_stats_df=pd.DataFrame.from_dict(all_stats_dic, orient='index', columns=cols )#['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max', 'mean_median_skew', 'full_spread', 'quartile_count', 'quartile_spread', 'quartile_factor', 'sd1_count', 'sd1_spread', 'sd1_factor', 'sd2_count', 'sd2_spread', 'sd2_factor', 'sd3_count', 'sd3_spread', 'sd3_factor'])
all_stats_df.sort_index(inplace=True)
all_stats_df.to_csv(outdir+'/all_og_info.tsv', sep='\t')
