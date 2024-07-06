#!/bin/sh
## set reference for STAR and arriba respectively
function RunArriba() {
    num_args=$#
    if [ $num_args -eq 0 ]; then
        echo "No arguments provided."
        exit 1
    fi
    STAR \
        --runThreadN 40 \
        --genomeDir /Path/to/reference/ref_Arriba/index \
        --genomeLoad NoSharedMemory \
        --readFilesIn $1/qc/$2.fastq.gz \
        --readFilesCommand zcat \
        --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
        --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |
    arriba \
        -x /dev/stdin \
        -o $1/fusiontsv/$2_fusions.tsv \
        -O $1/fusiontsv/$2_fusions.discarded.tsv \
        -a /Path/to/reference/ref_Arriba/refgenome/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        -g /Path/to/reference/ref_Arriba/refgenome/Homo_sapiens.GRCh38.106.gtf \
        -b /Path/to/tools/arriba/arriba_v2.3.0/database/blacklist_hg38_GRCh38_v2.3.0.tsv.gz  \
        -k /Path/to/tools/arriba/arriba_v2.3.0/database/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz \
        -t /Path/to/tools/arriba/arriba_v2.3.0/database/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz \
        -p /Path/to/tools/arriba/arriba_v2.3.0/database/protein_domains_hg38_GRCh38_v2.3.0.gff3 \
        -f duplicates,no_coverage,merge_adjacent

# $1 output dir
# $2 inputname
}




