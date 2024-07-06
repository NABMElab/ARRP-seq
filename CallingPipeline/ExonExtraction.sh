#!/bin/sh

function ExonExtraction() {
    num_args=$#
    if [ $num_args -eq 0 ]; then
        echo "No arguments provided."
        exit 1
    fi
    Gtf=$1
    Transcript=$2
    gene=$3
    loc=$4
    inputname=$5
    echo "${loc}/${inputname}_exon.csv"
    if [ -f "${loc}/${inputname}_exon.csv" ]; then
        if [ -z "$(cat ${loc}/${inputname}_exon.csv | grep "${Transcript}")" ];then
            zcat ${Gtf} | grep "${Transcript}" | awk -F '\t' -v OFS=',' '$3 == "exon" && $9 ~ /transcript_id[[:space:]]+"'${Transcript}'"/ { match($9, /exon_number ([0-9]+)/, num); chr = gensub(/^chr([0-9]+).*$/, "\\1", "g", $1); print "'${Transcript}'", "'${gene}'", chr, $4, $5, $7, num[1] }' >> ${loc}/${inputname}_exon.csv
            
        fi
    else
        zcat ${Gtf} | grep "${Transcript}" | awk -F '\t' -v OFS=',' '$3 == "exon" && $9 ~ /transcript_id[[:space:]]+"'${Transcript}'"/ { match($9, /exon_number ([0-9]+)/, num); chr = gensub(/^chr([0-9]+).*$/, "\\1", "g", $1); print "'${Transcript}'", "'${gene}'", chr, $4, $5, $7, num[1] }' > ${loc}/${inputname}_exon.csv
    fi


}
