#!/bin/sh

function FindFPAndSpacer() {
    flag=0
    direction=0
    exon=0
    FP=0
    record_validation=0
    Targets=$1
    Targetgene=$2
    seq=$3

    for record in `cat ${Targets} | grep ${Targetgene}`;do
        FP=`echo ${record} | cut -d ',' -f 4 | tr -cd 'AGCT'`
        if [[ "${seq}" =~ ^"${FP}".* ]]; then
            direction=`echo ${record} | cut -d ',' -f 2`
            exon=`echo ${record} | cut -d ',' -f 3`
            flag=1
        fi
    done
    echo "$direction $exon $FP $flag"
    
# $1 Targets file
# $2 target gene
# $3 seq
}

