# ARRP-seq

This repository contains code for processing fusion genes with ARRP-seq. The code is designed to streamline the analysis of fusion gene data, facilitating research and development in this area.

## Usage

The main script for running the analysis is `pipeline.sh`. To use the code, please follow these steps:

1. **Modify the Configuration**:
   Open `pipeline.sh` and update the following variables:
   - `TargetsFile`: The location of your primer design whose detailed demand could be found on section tagets file.
   - `database`: The desired location for your database built for local blast.
   - `Ref`: The desired location for you reference file, e.g. Homo_sapiens.GRCh38.cdna.all.fa.

  Open `runArriba.sh` and update the variables for refernce including `--genomeDir` for STAR and `-a`,`-g` for Arriba. some of these files are included in files installed.
  
  Open `Partners.sh` and update the following variables:
  - `Gtf`: The location of your gtf file, `.gz` is demanded.

2. **Run the Pipeline**:
   Execute the script using the command:
   ```bash
   sh pipeline.sh -i input -o output
