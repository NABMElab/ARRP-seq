# ARRP-seq

This repository contains code for processing fusion genes with ARRP-seq. The code is designed to streamline the analysis of fusion gene data, facilitating research and development in this area.

## Usage

The main script for running the analysis is `pipeline.sh`. To use the code, please follow these steps:

1. **Modify the Configuration**:
   Open `pipeline.sh` and update the following variables:
   - `TargetsFile`: The location of your primer design whose detailed demand could be found on section tagets file.
   - `database`: The desired location for your database built for local blast.
   - `Ref`: The desired location for you reference file, e.g. Homo_sapiens.GRCh38.cdna.all.fa.
   - `pseudo`: The location for pseudo genes or other excluded genes.

  Open `runArriba.sh` and update the variables for refernce including `--genomeDir` for STAR and `-a`,`-g` for Arriba. some of these files are included in files installed.
  
  Open `Partners.sh` and update the following variables:
  - `Gtf`: The location of your gtf file, `.gz` is demanded.

open `RawFilter.py` and update the target gene designed. Please replace the NTRK1|NTRK2|NTRK3|FGFR2|FGFR3 with any other target genes covered by your p   anel in the code `o1 = out[(out["#gene1"].str.contains('NTRK1|NTRK2|NTRK3|FGFR2|FGFR3').astype("int")+out["gene2"].str.contains('NTRK1|NTRK2|NTRK3|FGFR2|FGFR3').astype("int"))== 1]`  

2. **Run the Pipeline**:
   Execute the script using the command:
   ```bash
   sh pipeline.sh -i input.fastq.gz -o output_location


## Targets file demand
### Gene Primer Design Table

| Gene        | Direction | Exon | FP                 | Start | End   | Primer             |
|-------------|-----------|------|--------------------|-------|-------|--------------------|
| TargetGene1 | +         | 1    | ATCGTAGCxxxxx   | 1000  | 1100  | ATCGTAGC            |
| TargetGene2 | -         | 2    | GCTAGCTAxxxxx   | 2000  | 2100  | GCTAGCTA            |
| TargetGene3 | +         | 3    | TTAGGCCAxxxxx   | 3000  | 3100  | TTAGGCCA            |

This file is in CSV format and includes the following content:

- **Gene**: Name of the target gene
- **Direction**: the direction for primer  designed based on the `+` or `-` strand of the chromosome
- **FP**: Contains the designed primer sequence and spacer sequence
- **Start**: Corresponds to the starting position on the chromosome
- **End**: Corresponds to the ending position on the chromosome
- **Primer**: The sequence of the designed primer sequence

## Example
test fastq file is uploaded on figshare that can be downloaded from https://doi.org/10.6084/m9.figshare.27093502.v1
run the `run.sh`, the reference version for all referecnce such as fastq, gtf and database for blast used in example is hg38. Target file is included in `Example` file.


