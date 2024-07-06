def reverse_complement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'T', 't':'A', 'c':'G', 'g':'C','N':'N','n':'n','-':'-','.':'.', '?':'?','_':'_'}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    return ''.join(bases[::-1])

import pandas as pd
import numpy as np
import csv
import sys

arriba_loc=str(sys.argv[1]) #tsv loc
filename=str(sys.argv[2])
csv_file = arriba_loc + '/tmp/' + filename + '.csv'
txt_file = arriba_loc + '/tmp/' + filename + '.txt'
data = pd.read_csv(arriba_loc+'/fusiontsv/'+filename+'_fusions.tsv',sep='\t')

### 15+15 brk point
data['seq1'],data['seq2'] = data['fusion_transcript'].str.split("|",1).str
data['seq2'] = data['seq2'].str.replace('|', '', regex=True)
seq1=data['seq1'].str[-15:]
seq2 = data['seq2'].str[0:15]

### check name: gene1_gene2_brk
data['brkpoint'] = pd.DataFrame(list(seq1.str.cat(seq2)),columns=['brkpoint'])


###4 cols: gene1, gene2, check name, brk
out = pd.concat([data['#gene1'],data['gene2']],axis=1)
out = pd.concat([out,data['brkpoint']],axis=1)

# FGFR/NTRKrelated fusions left only
o1 = out[(out["#gene1"].str.contains('NTRK1|NTRK2|NTRK3|FGFR2|FGFR3').astype("int")+out["gene2"].str.contains('NTRK1|NTRK2|NTRK3|FGFR2|FGFR3').astype("int"))== 1]
o1.dropna(axis=0,how="any",inplace=True)

for i in range(o1['#gene1'].shape[0]):
    if o1.iloc[i,1] == 'NTRK1' or o1.iloc[i,1] == 'NTRK2' or o1.iloc[i,1] == 'NTRK3' or o1.iloc[i,1] == 'NTRK4' or o1.iloc[i,1] == 'FGFR1' or o1.iloc[i,1] == 'FGFR2' or o1.iloc[i,1] == 'FGFR3':
        print(o1.iloc[i,2])
        o1.iloc[i,2] =reverse_complement(o1.iloc[i,2])
        tmp = o1.iloc[i,0]
        o1.iloc[i,0] = o1.iloc[i,1]
        o1.iloc[i,1] = tmp
o1['fusions'] = o1['#gene1'].str.cat(o1['gene2'],sep='@')
o1['fusions'] = o1['fusions'].str.cat(o1['brkpoint'],sep='@')

o1['brkpoint'] = o1['brkpoint'].str.upper()
o1['brkpoint'] = o1['brkpoint'].fillna('-')
o1['brkpoint'] = o1['brkpoint'].str.replace('?', '@', regex=True)
#o1.to_csv(csv_file,index=0)

file = []
for p,q in o1.iterrows():
    tag = ">" + str(q[3])
    file.append(tag)
    file.append(str(q[2]))

file = pd.DataFrame(file)
file.to_csv(txt_file,sep=' ', index=False, header=False, quoting=csv.QUOTE_NONE, escapechar='|')

