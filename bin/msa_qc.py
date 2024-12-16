
import pandas as pd
import sys

"""
for elem in ../simsa/OG*/msas/squeezed_merged_OG*.fasta ; do
identifier=$(basename "$elem" | grep -o 'OG[0-9]*')
mkdir $identifier
cp -n $elem $identifier/simsa_$identifier.fasta
done
"""
"""
for elem in OG*/simsa*.fasta ;do 
echo  $elem  
python3 ../../../achromatium_pipeline/bin/msa_qc.py $elem
done
"""

def write_fasta(out_name, rows):
    with open (out_name, 'w') as m:
        for row in rows:
            row = row.replace('\\n','\n').replace(' ','')
            row = row.replace('>'+ogn+'_contig','>'+ogn+'_AGuadaBins3_contig')
            m.write( row + '-\n')  #trailing gap to allow b2btools

msa = sys.argv[1] 
GAP_CHAR = '-'

min_occ=0.1
max_locc_gap=10
coverage= 0.6
assembly_occ= 38*0.8

ogn=msa.split('/')[-1].split('_')[-1].split('.')[0]

df = pd.read_csv(msa, header=None)
msa_df = pd.DataFrame({'label':df[0].iloc[::2].values, 'seq':df[0].iloc[1::2].values})
seq_df = pd.DataFrame(msa_df.seq.apply(list).tolist(), index=msa_df.label)

#check if seqs have all the same length
if seq_df.isnull().values.any() == True:
    raise ValueError("Sequences do not have all the same length, is this an MSA?")

# minimal occupancy 
SEQUENCES_COUNT, RESIDUES_COUNT = seq_df.shape
OCCUPANCY = 1 - (seq_df == GAP_CHAR).sum() / SEQUENCES_COUNT

locc = OCCUPANCY[OCCUPANCY < min_occ].index

#ignore initial/trailing gaps
numbers=locc.tolist()
indices = [i for i in range(len(numbers) - 1) if numbers[i + 1] - numbers[i] > 1]
start_index = max(0,indices[0] - max_locc_gap)
end_index =  min(indices[-1] + max_locc_gap, len(numbers))
low_occupancy_positions = numbers[:start_index] +numbers[indices[0]+1:indices[-1]+1] +numbers[end_index:]

#print("init",seq_df.shape) 
# Identify sequences causing gaps in low occupated positions
seq_df['low_occ_gaps'] = (seq_df.iloc[:, low_occupancy_positions] != GAP_CHAR).sum(axis=1)
seq_df['keep'] = seq_df['low_occ_gaps'] <= max_locc_gap
seq_df = seq_df.drop(columns=['low_occ_gaps'])

keep_true = seq_df[seq_df['keep']].drop(columns=['keep'])
keep_false = seq_df[~seq_df['keep']].drop(columns=['keep'])

#print('low occ filter',min_occ,max_locc_gap,keep_true.shape)
keep_true = keep_true.loc[:, (keep_true != GAP_CHAR).any()]
#print('dropempty cols',keep_true.shape)

# Coverage filter
keep_true['gap_ratios'] = (keep_true == GAP_CHAR).mean(axis=1)
keep_true['keep']=keep_true['gap_ratios'] < coverage
keep_true = keep_true.drop(columns=['gap_ratios'])

keep_true_true = keep_true[keep_true['keep']].drop(columns=['keep'])
keep_false_false = keep_true[~keep_true['keep']].drop(columns=['keep'])

#print('coverage filter',1-coverage,keep_true_true.shape)
keep_true_true = keep_true_true.loc[:, (keep_true_true != GAP_CHAR).any()]



#prep and write outfile
keep_true_true.index= keep_true_true.index + '\n'
keep_true_true.reset_index(drop=False, inplace=True)
#check assembly occ
keep_true_true['ass']=keep_true_true['label'].str.split('_').str[1]

noass=keep_true_true['ass'].nunique()
if keep_true_true['ass'].nunique() > assembly_occ:

    ####remove in reruN!!
    keep_true_true=keep_true_true[keep_true_true['ass'] != 'AspWMS1']

    keep_true_true = keep_true_true.drop(columns=['ass'])
    rows = keep_true_true.to_string(header=False,index=False,index_names=False).split('\n')
    mpath,mname = msa.rsplit('/',1)
    out_name=mpath+'/'+mname.split('_')[-1]
    write_fasta(out_name, rows)


    keep_false_false.index= keep_false_false.index + '\n'
    keep_false_false.reset_index(drop=False,inplace=True)
    ffrows = keep_false_false.to_string(header=False,index=False,index_names=False).split('\n')

    keep_false.index= keep_false.index + '\n'
    keep_false.reset_index(drop=False,inplace=True)
    frows = keep_false.to_string(header=False,index=False,index_names=False).split('\n')

    falserows=frows+ffrows
    false_name=mpath+'/removed_from_simsa_'+mname.split('_')[-1]
    write_fasta(false_name, falserows)
    sel='selected'
else:
    #print('not enough coverage of assemblies')
    #print(keep_true_true['ass'].nunique() ,assembly_occ)
    sel='excluded'


with open('msa_qc_log.tsv', 'a+') as log:
    line='\t'.join(str(x) for x in [ogn ,seq_df.shape ,min_occ,max_locc_gap,keep_true.shape, str(1-coverage),keep_true_true.shape,noass, sel])
    log.write(line +'\n')