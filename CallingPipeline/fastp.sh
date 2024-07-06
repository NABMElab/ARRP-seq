#!/bin/sh

function QC() {
    num_args=$#
    if [ $num_args -eq 0 ]; then
        echo "No arguments provided."
        exit 1
    fi
    input=$1
    output_loc=$2
    inputname=$3
    echo ${input}
    echo ${output_loc}
    echo ${inputname}
    fastp \
    -i ${input} \
    -o ${output_loc}/qc/${inputname}.fastq.gz \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 20 \
    --length_required 30 \
    --length_limit 151 \
    --adapter_sequence=auto
# $1 input file
# $2 output dir
# $3 inputname
}



