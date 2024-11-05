#!/bin/sh

func() {
    echo "Usage:"
    echo "Pipeline.sh [-i input] [-o outputDir] [-e]"
    echo "input, the path to the input file (support .fastq.gz only)"
    echo "optional -e, remain exon transcript"
    echo "outputDir, the path of destination"
    exit -1
}

module="all"

while getopts 'i:o:e' OPT; do
    case $OPT in
        i) input="$OPTARG";;
        o) output="$OPTARG";;
        e) module="exon";;
        ?) func;;
    esac
done

inputname=`basename $input .fastq.gz`
TargetsFile="/Path/To/targets.csv" # primer design files for reference
pseudo="/Path/To/pseudo.csv" # black list for genes such as pseudogene

cd ${output}


Code=$(dirname "$PWD")
source ${Code}/code/fastp.sh
source ${Code}/code/RunArriba.sh
source ${Code}/code/SpacerFilter.sh
source ${Code}/code/blast.sh
source ${Code}/code/Partners.sh
source ${Code}/code/ExonReCheck.sh

if [ ! -d ${output}/result ];then
    mkdir result
fi

if [ ! -d ${output}/qc ];then
    mkdir qc
fi
if [ ! -d ${output}/fusiontsv ];then
    mkdir fusiontsv
fi
if [ ! -d ${output}/tmp ];then
    mkdir tmp
fi

if [ -f ${output}/result/${inputname}_fusion_Raw.csv ];then
    rm -f ${output}/result/${inputname}_fusion_Raw.csv
fi

if [ -f ${output}/tmp/${inputname}.fasta ];then
    rm -f ${output}/tmp/${inputname}.fasta
fi

if [ -f "${output}/tmp/${inputname}_exon.csv" ];then
    rm -f ${output}/tmp/${inputname}_exon.csv
fi

# QC and outputs stored in output/qc/inputname
QC ${input} ${output} ${inputname}
### Raw arriba and outputs stored in output/fusiontsv/inputname_fusions.tsv
RunArriba ${output} ${inputname}
###
python ${Code}/code/RawFilter.py ${output} ${inputname}

# IPF+Spacer list for targeted genes

for reads in `cat ${output}/tmp/${inputname}.txt`; do
    if [ ${reads:0:1} = ">" -a ${reads:1} = "nan" ]; then
        echo "nan"  ## fusion name don't exist
        
    elif [ ${reads:0:1} = ">" -a ${reads:1} != "nan" ];then
        name=`echo ${reads:1} | sed -e 's/,/|/g'` ## FGFR1@part2-1,part2-2 -> FGFR1@part2-1|part2-2
        echo "Query" > ${output}/tmp/${inputname}_blasttmp.fasta  #blast require
        current_query=$reads
        echo ${reads} >> ${output}/tmp/${inputname}.fasta
        cat ${output}/tmp/${inputname}.fasta
        
    elif [ ${reads:0:1} != ">" -a ${reads:0:1} != "-" ]; then
        seq=`zcat ${input} | grep ${reads} | sort -n | uniq -c | sort -n | tail -n 1 | awk '{print $2}'`
        num=`zcat ${input} | grep ${reads} | sort -n | uniq -c | sort -n | tail -n 1 | awk '{print $1}'`
        count=`zcat ${input} | grep ${reads} | sort -n | uniq -c | sort -n | tail -n 1 | awk 'BEGIN{FS=" ";OFS=","}{print $1,$2}'`
        target=`echo ${name} | cut -d '@' -f 1`
        if [ ! -z "${seq}" ] && [ ${num} -gt "5" ];then   # filter:reads>5
            result=$(FindFPAndSpacer ${TargetsFile} ${target} ${seq})
            read -r direction exon FP flag <<< "$result" # flag =1 : start with IFP+Spacer
            if [ "${flag}" -eq 1 ];then
                echo ${seq} >> ${output}/tmp/${inputname}_blasttmp.fasta
                blastInput=${output}/tmp/${inputname}_blasttmp.fasta
                blastOutput=${output}/tmp/${inputname}_blasttmp.blast
                database="/Path/to/reference/blastn/h38cDNA"
                Ref="/Path/to/reference/ensemble/Homo_sapiens.GRCh38.cdna.all.fa"
                blastFlag=$(blast ${blastInput} ${blastOutput} ${database} ${Ref})
                blastModified=${output}/tmp/${inputname}_modified.csv
                
                if [ ${blastFlag} -eq 1 ];then
                ## check pseudos
                    awk -F ',' '
                    FNR==NR {
                        gene_names[$0] = 1;  # store names in array
                        next;
                      }
                      {
                        gene_symbol = $10;  # gene_symbol--col 10
                        if (!(gene_symbol in gene_names)) {
                          print $0;
                        }
                      }
                    ' ${pseudo} ${blastModified} > ${output}/tmp/${inputname}_modified1.csv
                    Partners ${name} ${output}/tmp/${inputname}_modified1.csv ${inputname}
                    partner1=`echo ${name} | cut -d '@' -f 2 | cut -d '|' -f 1 | cut -d '(' -f 1`
                    partner2=`echo ${name} | cut -d '@' -f 2 | cut -d '|' -f 2 | cut -d '(' -f 1`
                    rm -f ${output}/tmp/${inputname}_modified1.csv
                    rm -f ${output}/tmp/${inputname}_modified.csv
                    if [ -f ${output}/tmp/${name}_single.csv ];then
                        Rows=`cat ${output}/tmp/${name}_single.csv | wc -l`
                        echo "Rows: ${Rows}"
                        if [ $Rows -eq 2 ]; then
                            echo -n "${name},${num},${exon},${direction},${seq}," >> ${output}/result/${inputname}_fusion_Raw.csv
                            head -n 1 "${output}/tmp/${name}_single.csv" | tr -d '\n' >> ${output}/result/${inputname}_fusion_Raw.csv
                            echo -n "," >> ${output}/result/${inputname}_fusion_Raw.csv
                            sed -n '2p' "${output}/tmp/${name}_single.csv" >> ${output}/result/${inputname}_fusion_Raw.csv
                            echo ${seq} >> ${output}/tmp/${inputname}.fasta
                        elif [ $Rows -eq 3 ]; then
                            breakseq=`echo ${name} | awk -F "@" '{print $3}'`
                            echo -n "${target}@${partner1}@${breakseq},${num},${exon},${direction},${seq}," >> ${output}/result/${inputname}_fusion_Raw.csv
                            head -n 1 "${output}/tmp/${name}_single.csv" | tr -d '\n' >> ${output}/result/${inputname}_fusion_Raw.csv
                            echo -n "," >> ${output}/result/${inputname}_fusion_Raw.csv
                            sed -n '2p' "${output}/tmp/${name}_single.csv" >> ${output}/result/${inputname}_fusion_Raw.csv
                            sed -i '$d' ${output}/tmp/${inputname}.fasta
                            echo “>${target}@${partner1}@${breakseq} >> ${output}/tmp/${inputname}.fasta
                            echo "${target}@${partner2}@${breakseq},${num},${exon},${direction},${seq}," >> ${output}/result/${inputname}_fusion_Raw.csv
                            head -n 1 "${output}/tmp/${name}_single.csv" | tr -d '\n' >> ${output}/result/${inputname}_fusion_Raw.csv
                            echo -n "," >> ${output}/result/${inputname}_fusion_Raw.csv
                            tail -n 1 "${output}/tmp/${name}_single.csv" >> ${output}/result/${inputname}_fusion_Raw.csv
                            echo "" >> ${output}/result/${inputname}_fusion_Raw.csv
                            echo “>${target}@${partner2}@${breakseq} >> ${output}/tmp/${inputname}.fasta
                            echo ${seq} >> ${output}/tmp/${inputname}.fasta
                        else
                            sed -i '$d' ${output}/tmp/${inputname}.fasta
                        fi
                    else
                        sed -i '$d' ${output}/tmp/${inputname}.fasta
                    fi
                else
                    sed -i '$d' ${output}/tmp/${inputname}.fasta
                fi
            else
                sed -i '$d' ${output}/tmp/${inputname}.fasta
            fi
        else
            sed -i '$d' ${output}/tmp/${inputname}.fasta
        fi
        rm -f ${output}/tmp/${name}_single.csv
    fi
    #clear tem files
    rm -f ${output}/tmp/${inputname}_blasttmp.fasta
    rm -f ${output}/tmp/${inputname}_blasttmp.blast
done

## table header
if [ -s "${output}/result/${inputname}_fusion_Raw.csv" ]; then
    sed -i "1igene1_gene2_brk,reads,exon,direction,Sequence,5'Trans,5'Percent,5'Length,5'Mismatch,5'Start,5'End,5'ChrStart,5'ChrEnd,5'Evalue,5'gene, 5'Chr,5'PosStart,5'PosEnd,5'Dir,5'Version,5'Type,3'Trans,3'Percent,3'Length,3'Mismatch,3'Start,3'End,3'ChrStart,3'ChrEnd,3'Evalue,3'gene,3'Chr,3'PosStart,3'PosEnd,3'Dir,3'Version,3'Type" ${output}/result/${inputname}_fusion_Raw.csv
fi
if [ -s "${output}/tmp/${inputname}_exon.csv" ]; then
    sed -i "1iTranscript,Gene,Chromosome,Start,End,Dir,Exon" ${output}/tmp/${inputname}_exon.csv
    python ${Code}/code/Concat.py ${output} ${inputname} ${TargetsFile}
fi

if [ -s "${output}/tmp/${inputname}_ExonFiltered.csv" ]; then
    ExonReCheck ${input} ${output}
fi

if [ -s "${output}/tmp/${inputname}_fusions_1stExonPass.csv" -o -s "${output}/tmp/${inputname}_ExonPass.csv" ];then
    python ${Code}/code/ExonReCat.py  ${output} ${inputname}
fi

