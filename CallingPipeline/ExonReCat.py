import numpy as np
import pandas as pd 
import csv
import sys

output = str(sys.argv[1])
inputname = str(sys.argv[2])
First_Exon_File = output + '/tmp/' + inputname + '_fusions_1stExonPass.csv'
ExonPassed_File = output + '/tmp/' + inputname + '_ExonPass.csv'
Final_File = output + '/result/' + inputname + '_fusions_final.csv'

try:
    First_Exon = pd.read_csv(First_Exon_File,sep=',')
except FileNotFoundError:
    First_Exon = pd.DataFrame()

try:
    ExonPassed = pd.read_csv(ExonPassed_File,sep=',')
    ## clear duplicates
    sorted = ExonPassed.sort_values("3'Evalue").sort_values('Sequence',key=lambda col:col.str.upper())
    duplicated = sorted[sorted.duplicated('Sequence',keep=False)]
    dups = duplicated.drop_duplicates('Sequence',keep='first')
    num = 0
    for dup in dups['Sequence']:
        fusions = duplicated[duplicated['Sequence']==dup]
        evalue = fusions["3'Evalue"].min()
        evalue_dup = fusions[fusions["3'Evalue"] == evalue].shape[0]
        if evalue_dup == 1:
            sorted.drop(fusions.index[1:],inplace=True)
        else:
            sorted.drop(fusions.index,inplace=True)
    Data_simplified_2 = sorted[sorted["5'End"]>=sorted["3'Start"]-1]
    cols_to_keep_2 = ["gene1_gene2_brk", "reads", "exon", "direction", "Sequence", "5'gene","5'Exon", "5'ChrPosition", "3'gene", "3'Exon", "3'ChrPosition",
                    "5'Trans", "5'Start", "5'End",  "5'Type",
                    "3'Trans", "3'Start", "3'End",  "3'Type"]

    Data_simplified_2 = Data_simplified_2[cols_to_keep_2]
    Data_simplified_2['gene1_gene2_brk'] = Data_simplified_2['gene1_gene2_brk'].str.replace('@', '-')
    Data_simplified_2.columns.values[0] = 'Fusions'
except FileNotFoundError:
    Data_simplified_2 = pd.DataFrame()


Final_Concat= pd.concat([First_Exon, Data_simplified_2], ignore_index=True)

if not Final_Concat.empty:
    Final_Concat.to_csv(Final_File,header=True,index=False)


