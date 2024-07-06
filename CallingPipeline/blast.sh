#!/bin/sh

function blast() {
    num_args=$#
    if [ $num_args -eq 0 ]; then
        echo "No arguments provided."
        exit 1
    fi
    blastflag=0
    blastInput=$1
    blastOutput=$2
    database=$3
    Ref=$4
    loc=`dirname ${blastOutput}`
    Filename=`basename "${blastInput}" | sed 's/_blasttmp\.fasta$//'`
    modified_File=${loc}/${Filename}_modified.csv
    blastn \
        -query ${blastInput} \
        -out ${blastOutput} \
        -db ${database} \
        -outfmt 6 \
        -evalue 0.02 \
        -num_threads 20 \
        -task blastn-short
    if [ -s ${blastOutput} ];then
        awk '
        FNR==NR{
            if ($0 ~ /^>/) {
                split($0, parts, " ");
                split($1, id,">");
                if (split(parts[6], transcript_biotype_parts, ":") && transcript_biotype_parts[1] == "transcript_biotype") {
                    transcript_biotypes[id[2]] = transcript_biotype_parts[2];}
                else {
                    transcript_biotypes[id[2]] = "NA"
                }
                split(parts[3],positions,":");
                if (split(parts[7], gene_symbol_parts, ":") && gene_symbol_parts[1] == "gene_symbol") {
                    gene_symbols[id[2]] = gene_symbol_parts[2];
                    chrs[id[2]] = positions[3];
                    starts[id[2]] = positions[4];
                    ends[id[2]] = positions[5];
                    dirs[id[2]] = positions[6];
                    versions[id[2]] = positions[2];
                }
            }
        next;
        }
        {
            enst_id = $2;
            gene_symbol = gene_symbols[enst_id];
            start = starts[enst_id];
            end = ends[enst_id];
            dir = dirs[enst_id];
            chr = chrs[enst_id];
            version = versions[enst_id];
            transcript = transcript_biotypes[enst_id];
            print $2 "," $3 "," $4 "," $5 "," $7 "," $8 "," $9 "," $10 "," $11 "," gene_symbol "," chr "," start "," end "," dir "," version "," transcript
        }
        ' ${Ref} ${blastOutput} > ${modified_File}
        blastflag=1
        echo ${blastflag}
    else
        echo ${blastflag}
    fi
}
