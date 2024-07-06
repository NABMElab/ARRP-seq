#!/bin/sh
Code=$(dirname "$PWD")
source ${Code}/code/ExonExtraction.sh

function Partners() {
    num_args=$#
    if [ $num_args -eq 0 ]; then
        echo "No arguments provided."
        exit 1
    fi
    name=$1
    blastModified=$2
    inputname=$3
    loc=`dirname ${blastModified}`
    partner1=`echo ${name} | cut -d '@' -f 2 | cut -d '|' -f 1 | cut -d '(' -f 1`
    partner2=`echo ${name} | cut -d '@' -f 2 | cut -d '|' -f 2 | cut -d '(' -f 1`
    panelgene=`echo ${name} | cut -d '@' -f 1`
    gene1=`cat ${blastModified} | grep ${panelgene}`
    Gtf="/data/wuy/reference/GENECODE/gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz"
    ## get gene1(FGFR/NTRK position
    if [ ! -z "${gene1}" ]; then
        pc_head=`cat ${blastModified} | grep -w ${panelgene} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 | awk 'BEGIN{FS=",";OFS="-"}{print $9}'` # evalue Ascending order --9  transcript descending order--15
        min_head=`cat ${blastModified} | grep -w ${panelgene} | sort -k 9g -k 16r -t , | head -n 1 | awk 'BEGIN{FS=",";OFS="-"}{print $9}'`
        if [ ${pc_head} = ${min_head} ]; then
            gene1_pos=`cat ${blastModified} | grep -w ${panelgene} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 | awk 'BEGIN{FS=",";OFS="-"}{print $5,$6}'`
            Transcript=$(cat ${blastModified} | grep -w ${panelgene} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 | awk -F "," '{print $1}')
            ExonExtraction ${Gtf} ${Transcript} ${panelgene} ${loc} ${inputname}
            cat ${blastModified} | grep -w ${panelgene} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 > ${loc}/${name}_single.csv
        else
            gene1_pos=`cat ${blastModified} | grep -w ${panelgene} | sort -k 9g -k 16r -t , | head -n 1 | awk 'BEGIN{FS=",";OFS="-"}{print $5,$6}'`
            Transcript=$(cat ${blastModified} | grep -w ${panelgene} | sort -k 9g -k 16r -t , | head -n 1 | awk -F "," '{print $1}')
            ExonExtraction ${Gtf} ${Transcript} ${panelgene} ${loc} ${inputname}
            cat ${blastModified} | grep -w ${panelgene} | sort -k 9g -k 16r -t , | head -n 1 > ${loc}/${name}_single.csv
        fi
    else
        gene1_pos=' '
        echo "NA,NA,NA,NA,NA,NA,NA,NA,NA,${gene1},NA,NA,NA,NA,NA" > ${loc}/${name}_single.csv
    fi
    
    
    ## most genefusion partner1=partner2
    if [ ${partner1} = ${partner2} ];then
        gene2=`cat ${blastModified} | grep -w ${partner1}`
            if [ ! -z "${gene2}" ]; then
                pc_head=`cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 | awk 'BEGIN{FS=",";OFS="-"}{print $9}'`
                min_head=`cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , | head -n 1 | awk 'BEGIN{FS=",";OFS="-"}{print $9}'`
                if [ "${pc_head}" = "${min_head}" ]; then
                    Transcript=$(cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 | awk -F "," '{print $1}')
                    ExonExtraction ${Gtf} ${Transcript} ${partner1} ${loc} ${inputname}
                    cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 >> ${loc}/${name}_single.csv
                else
                    Transcript=$(cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , |  head -n 1 | awk -F "," '{print $1}')
                    ExonExtraction ${Gtf} ${Transcript} ${partner1} ${loc} ${inputname}
                    cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , | head -n 1 >> ${loc}/${name}_single.csv
                fi
            fi
        
        
    ## rare genefusion partner1!=partner2
    else
        gene2_1=`cat ${blastModified} | grep -w ${partner1}`
        gene2_2=`cat ${blastModified} | grep -w ${partner1}`
        if [ ! -z "${gene2_1}" ] && [ ! -z "${gene2_2}" ];then
            pc_head1=`cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 | awk 'BEGIN{FS=",";OFS="-"}{print $9}'`
            min_head1=`cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , | head -n 1 | awk 'BEGIN{FS=",";OFS="-"}{print $9}'`
            pc_head2=`cat cat ${blastModified} | grep -w ${partner2} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 | awk 'BEGIN{FS=",";OFS="-"}{print $9}'`
            min_head2=`cat ${blastModified} | grep -w ${partner2} | sort -k 9g -k 16r -t , | head -n 1 | awk 'BEGIN{FS=",";OFS="-"}{print $9}'`
        
            #### partner1
            if [ ${pc_head1} = ${min_head1} ]; then
                Transcript=$(cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 | awk -F "," '{print $1}')
                ExonExtraction ${Gtf} ${Transcript} ${partner1} ${loc} ${inputname}
                cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 >> ${loc}/${name}_single.csv
            else
                Transcript=$(cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , |  head -n 1 | awk -F "," '{print $1}')
                ExonExtraction ${Gtf} ${Transcript} ${partner1} ${loc} ${inputname}
                cat ${blastModified} | grep -w ${partner1} | sort -k 9g -k 16r -t , | head -n 1 >> ${loc}/${name}_single.csv
            fi
            #### partner2
            if [ ${pc_head2} = ${min_head2} ]; then
                Transcript=$(cat ${blastModified} | grep -w ${partner2} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 | awk -F "," '{print $1}')
                ExonExtraction ${Gtf} ${Transcript} ${partner2} ${loc} ${inputname}
                cat ${blastModified} | grep -w ${partner2} | sort -k 9g -k 16r -t , | grep -w "protein_coding" | head -n 1 >> ${loc}/${name}_single.csv
            else
                Transcript=$(cat ${blastModified} | grep -w ${partner2} | sort -k 9g -k 16r -t , | head -n 1 | awk -F "," '{print $1}')
                ExonExtraction ${Gtf} ${Transcript} ${partner2} ${loc} ${inputname}
                cat ${blastModified} | grep -w ${partner2} | sort -k 9g -k 16r -t , | head -n 1 >> ${loc}/${name}_single.csv
            fi
        fi
    fi

    
}
