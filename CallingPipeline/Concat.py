def Chromosome_Position_Get(Tmp,Start,End):
    Tmp.index=(range(len(Tmp)))
    Tmp['SeqEnd'] = Tmp['ExonLength'].cumsum()
    Tmp['SeqStart'] = Tmp['SeqEnd']- Tmp['ExonLength'] + 1

    if Tmp.iloc[0,5] == '+':
        ChrStart = Start- \
                (Tmp[(Tmp['SeqStart']<= Start) & (Tmp['SeqEnd']>= Start)]['SeqStart'].astype(int)).values[0]+ \
                (Tmp[(Tmp['SeqStart']<= Start) & (Tmp['SeqEnd']>= Start)]['Start'].astype(int)).values[0]
        ExonStart = (Tmp[(Tmp['SeqStart']<= Start) & (Tmp['SeqEnd']>= Start)]['Exon'].astype(int)).values[0]
        
        ChrEnd = End - \
                (Tmp[(Tmp['SeqStart'] <= End) & (Tmp['SeqEnd'] >= End)]['SeqStart'].astype(int)).values[0] + \
                (Tmp[(Tmp['SeqStart'] <= End) & (Tmp['SeqEnd'] >= End)]['Start'].astype(int)).values[0]
        ExonEnd = (Tmp[(Tmp['SeqStart'] <= End) & (Tmp['SeqEnd'] >= End)]['Exon'].astype(int)).values[0]


    else:
        ChrStart = (Tmp[(Tmp['SeqStart']<= Start) & (Tmp['SeqEnd']>= Start)]['End'].astype(int)).values[0]- \
                (Start- \
                (Tmp[(Tmp['SeqStart']<= Start) & (Tmp['SeqEnd']>= Start)]['SeqStart'].astype(int)).values[0])
        ExonStart = (Tmp[(Tmp['SeqStart']<= Start) & (Tmp['SeqEnd']>= Start)]['Exon'].astype(int)).values[0]

                
        ChrEnd = (Tmp[(Tmp['SeqStart'] <= End) & (Tmp['SeqEnd'] >= End)]['End'].astype(int)).values[0] - \
                (End - \
                (Tmp[(Tmp['SeqStart'] <= End) & (Tmp['SeqEnd'] >= End)]['SeqStart'].astype(int)).values[0])
        ExonEnd = (Tmp[(Tmp['SeqStart'] <= End) & (Tmp['SeqEnd'] >= End)]['Exon'].astype(int)).values[0]

    return ChrStart, ChrEnd, ExonStart, ExonEnd

def Exon_filter(Primer, Primer_dir,Gene, ChrStart, ChrEnd):
    data = Primer[(Primer['Gene'] == Gene) & (Primer['Direction'] == Primer_dir)]
    primer = data[(data[['End', 'Start']].max(axis=1) <= max(ChrStart, ChrEnd)) & (data[['End', 'Start']].min(axis=1) >= min(ChrStart, ChrEnd))]
    return primer


import numpy as np 
import pandas as pd 
import csv
import sys

output = str(sys.argv[1])
inputname = str(sys.argv[2])
Primer_loc = str(sys.argv[3])
Raw_data_loc = output + '/result/' + inputname + '_fusion_Raw.csv'
exon_loc = output + '/tmp/' + inputname + '_exon.csv'
FinalFile = output + '/tmp/' + inputname + '_fusions_1stExonPass.csv'


exon = pd.read_csv(exon_loc)
exon['ExonLength']= exon['End'] - exon['Start']+1
Primers = pd.read_csv(Primer_loc)
Raw_data = pd.read_csv(Raw_data_loc)
cols_to_keep = ["gene1_gene2_brk", "reads", "exon", "direction", "Sequence", 
                "5'Trans", "5'Start", "5'End", "5'gene", "5'Type", 
                "3'Trans", "3'Start", "3'End", "3'gene", "3'Type","3'Evalue"]
Data_simplified = Raw_data[cols_to_keep]
ExonCheck = []

Data_simplified["5'ChrPosition"] = 0
Data_simplified["3'ChrPosition"] = 0
Data_simplified["ExonFlag"] = 0
Data_simplified["5'Exon"] = 0
Data_simplified["3'Exon"] = 0


for i in range(Raw_data["5'Trans"].shape[0]):
    Primer_dir = Raw_data.iloc[i,3]
    Chr_5 = Raw_data.iloc[i,15]
    Dir_5 = Raw_data.iloc[i,18]
    Trans_5 = Raw_data.iloc[i,5]
    Gene_5 = Raw_data.iloc[i,14]
    Tmp_5 = exon[exon['Transcript'] == Trans_5]
    Start_5 = Raw_data.iloc[i,11]
    End_5 = Raw_data.iloc[i,12]

    ChrStart_5, ChrEnd_5, ExonStart_5, ExonEnd_5 = Chromosome_Position_Get(Tmp_5,Start_5,End_5)
    Data_simplified.iloc[i,16] = 'Chr'+str(Chr_5)+':'+str(ChrStart_5)+':'+str(ChrEnd_5)+':'+str(Dir_5)
    Data_simplified.iloc[i,19] = ExonEnd_5
    PrimerCheck = Exon_filter(Primers, Primer_dir,Gene_5, ChrStart_5, ChrEnd_5)
    ExonFlag = PrimerCheck.shape[0]
    if ExonFlag == 1:
        Data_simplified.iloc[i,18] = 1
    elif ExonFlag > 1:
        PrimerCheck['Gene'] = Raw_data.iloc[i,0]
        PrimerCheck = PrimerCheck[PrimerCheck['Exon'] != Raw_data.iloc[i,2]]
        if len(ExonCheck) == 0:
            ExonCheck = PrimerCheck
        else:
            ExonCheck = pd.concat([ExonCheck, PrimerCheck], ignore_index=True)

    Chr_3 = Raw_data.iloc[i,31]
    Dir_3 = Raw_data.iloc[i,34]
    Trans_3 = Raw_data.iloc[i,21]
    Tmp_3 = exon[exon['Transcript'] == Trans_3]
    Start_3 = Raw_data.iloc[i,27]
    End_3 = Raw_data.iloc[i,28]

    ChrStart_3, ChrEnd_3, ExonStart_3, ExonEnd_3 = Chromosome_Position_Get(Tmp_3,Start_3,End_3)
    Data_simplified.iloc[i,17] = 'Chr'+str(Chr_3)+':'+str(ChrStart_3)+':'+str(ChrEnd_3)+':'+str(Dir_3)
    Data_simplified.iloc[i,20] = ExonEnd_3


## clear exon number
Exondup = Data_simplified[Data_simplified['ExonFlag'] == 0]
elements_in_Exondup = set(Exondup['gene1_gene2_brk'])
counts = Data_simplified['gene1_gene2_brk'].value_counts() ## counts for each fusion
Exon_filtered = Data_simplified
for element, count in counts.items():
    if element in elements_in_Exondup and count < 2:
        # delete rows once in dataframe2 only
        Data_simplified = Data_simplified[Data_simplified['gene1_gene2_brk'] != element]
Data_Exonfiltered = Exon_filtered[~Exon_filtered.isin(Data_simplified)].dropna()


## clear duplicates
sorted = Data_simplified.sort_values("3'Evalue").sort_values('Sequence',key=lambda col:col.str.upper())
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

if not Data_simplified_2.empty:
    Data_simplified_2.to_csv(FinalFile,header=True,index=False)
    
if not Data_Exonfiltered.empty:
    File_exonFiltered = output + '/tmp/' + inputname + '_ExonFiltered.csv'
    File_ExonCheck = output + '/tmp/' + inputname + '_ExonReCheck.csv'
    Data_Exonfiltered.to_csv(File_exonFiltered,header=True,index=False)
    ExonCheck.to_csv(File_ExonCheck,header=True,index=False)
