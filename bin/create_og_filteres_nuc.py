import pandas as pd
import sys
import math
import glob
#from pathlib import Path
from Bio import SeqIO
from collections import defaultdict

#module load Biopython/1.79-foss-2022a
"""
for elem in OG*/OG*.fasta ;do 
echo  $elem  
python3 ../../../achromatium_pipeline/bin/create_og_filteres_nuc.py $elem ../../genomes_clean/ 
done
"""

aa_file=sys.argv[1]
nuc_folderp=sys.argv[2]
nuc_folders=glob.glob(nuc_folderp+'/*/*.ffn')

nuc_files_dic={}
for nf in nuc_folders:
    name = nf.split("/")[-1].split(".")[0].replace('_',"")
    nuc_files_dic[name]=nf


#outfile= 'ognuc/'+ aa_file.split("/")[-1]
outfile = aa_file.split(".")[0] + ".fna"


with open(aa_file, 'r') as af:
    aa_dic = SeqIO.to_dict(SeqIO.parse(af, "fasta"))

aa_keys = sorted(aa_dic.keys())
aa_key_dic = defaultdict(list)
for elem in aa_keys:
    ass =elem.split('_')[1]
    match_id =elem.split('_',1)[1]
    aa_key_dic[ass].append((match_id,elem))

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

outdic={} #aa_id, nucseq
for key in aa_key_dic.keys():
    if key =="contig":
        print('wtf is contig')
        continue

    with open(nuc_files_dic[key], 'r') as nf:
        nuc_dic = to_dict_remove_dups(SeqIO.parse(nf, "fasta"))
    aa_oi = aa_key_dic[key]
    for elem in aa_oi:
        match_id, full_id =elem
        try:
            mseq = nuc_dic[match_id]
            outdic[full_id]=mseq.seq
        except:
            print("nomatch" ,match_id)

#print(outdic)
with open(outfile, 'w') as of:
    for elem in outdic.keys():
        line = '>'+elem+'\n'+ str(outdic[elem])
        of.write(line+'\n')
#print(outfile)



