#!/bin/sh
Code=$(dirname "$PWD")

function ExonReCheck() {
    num_args=$#
    if [ $num_args -eq 0 ]; then
        echo "No arguments provided."
        exit 1
    fi
    input=$1
    output=$2
    inputname=`basename $input .fastq.gz`
    ExonFiltered=${output}/tmp/${inputname}_ExonFiltered.csv
    ExonReCheck=${output}/tmp/${inputname}_ExonReCheck.csv
    ExonPass=${output}/tmp/${inputname}_ExonPass.csv
    while IFS=',' read -r fusion reads exon direction _; do
        if [ "$fusion" != "gene1_gene2_brk" ]; then
            record=`cat ${ExonFiltered} | grep ${fusion}`
            extracted_sequence=$(echo "$fusion" | rev | cut -d '@' -f 1 | rev | tr '[:lower:]' '[:upper:]')
            exons=`cat ${ExonReCheck} | grep ${extracted_sequence}`
            ExonReCheck_Flag=0
            while IFS=',' read -r Gene Direction Exon FP _; do
                ExonCheckReads=$(zcat ${input} | grep ${extracted_sequence} | grep ${FP} | sort -n | uniq -c | sort -n | tail -n 1)
                if [ -n "$ExonCheckReads" ]; then
                    ExonReCheck_Flag=$((ExonReCheck_Flag + 1))
                fi
            done < ${ExonReCheck}
            rows=$(cat ${ExonReCheck} | wc -l)
            if [ "${rows}" -eq "$((ExonReCheck_Flag + 1))" ]; then
               printf "%s\n"  "${record}"
            fi
        fi
    done < ${ExonFiltered} > ${ExonPass}
    
    if [ -s "${ExonPass}" ]; then
        sed -i "1igene1_gene2_brk,reads,exon,direction,Sequence,5'Trans,5'Start,5'End,5'gene,5'Type,3'Trans,3'Start,3'End,3'gene,3'Type,3'Evalue,5'ChrPosition,3'ChrPosition,ExonFlag,5'Exon,3'Exon" ${ExonPass}
    else
        rm -f ${ExonPass}
    fi
}
