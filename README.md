```
 ▄▄▄▄▄▄▄▄▄▄▄  ▄▄        ▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄ 
▐░░░░░░░░░░░▌▐░░▌      ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
▐░█▀▀▀▀▀▀▀▀▀ ▐░▌░▌     ▐░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌
▐░▌          ▐░▌▐░▌    ▐░▌▐░▌          ▐░▌       ▐░▌▐░▌       ▐░▌
▐░▌          ▐░▌ ▐░▌   ▐░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌
▐░▌          ▐░▌  ▐░▌  ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
▐░▌          ▐░▌   ▐░▌ ▐░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ 
▐░▌          ▐░▌    ▐░▌▐░▌▐░▌          ▐░▌          ▐░▌          
▐░█▄▄▄▄▄▄▄▄▄ ▐░▌     ▐░▐░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░▌          ▐░▌          
▐░░░░░░░░░░░▌▐░▌      ▐░░▌▐░░░░░░░░░░░▌▐░▌          ▐░▌          
 ▀▀▀▀▀▀▀▀▀▀▀  ▀        ▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀            ▀           
```                                                                 


cNePP  
Overview

This repository contains a comprehensive pipeline for immuno-oncology neoepitope multipe feature analysis. The pipeline automates various tasks, from preparing data folders to processing and analyzing genomic data, including SNV (Single Nucleotide Variants), indels (insertions and deletions), fusions, HLA (Human Leukocyte Antigen) typing, neoantigen prediction, and more. It supports different types of genomic data and various steps in the analysis process.
Features

    HLA typing: Process and analyze HLA typing data.
    SNV/Indel/Fusion: Analyze single nucleotide variants, indels, and fusion genes.
    Neoantigen Prediction: Predict potential neoantigens from genomic data.
    Excel Output: Convert results into XLSX format for easy reporting and sharing.
    Flexible Configuration: Customize various options like run IDs, steps, data type, and more via command-line arguments.

Prerequisites

    Bash shell
    Access to the necessary genomic data files
    Required dependencies for each step of the pipeline (see individual steps below)
    GitHub repository and script access

Usage

To run the pipeline, use the following command:

bash <script_name>.sh -r <runID> -s <steps> [additional options]

Options

    -r runID: Unique identifier for the run.
    -e netMHCpanID: NetMHCpan version (e.g., netMHCpan4_1, netMHCstabpan).
    -p patientID: Patient identifier, defaults to the value of runID if not provided.
    -t tcga: TCGA data type (e.g., RNAseq, TCGA-READ).
    -s steps: Specify the pipeline steps to execute. The available steps include:
        s0: Prepare folder
        s1a_HLA: HLA typing step 1a
        s1b_HLA: HLA typing step 1b
        s2_snv: SNV analysis
        s3_add_expression: Add expression data
        s8a_filter: Filter data
        s8b_xlsx_to_public: Convert data to public format
        i4a_indel_predict: Indel prediction
        f1a: Run STAR and Arriba (fusion detection)
        f2: Prepare HLA data
        f3: Neoantigen prediction
        f4: MER21 prediction
        f5: Convert to XLSX format
        More steps listed in the script.
    -d debug: Enable debug mode for troubleshooting.
    -h hlaID: Specify the HLA ID (e.g., NCT_IP for internal patient data).
    -a dataType: Specify the data type (e.g., RNAbamOpt).
    -c cgiTumorType: Specify the tumor type (e.g., CANCER).
    -b vcfOnly: Process only VCF files (e.g., origin, pathology, promise).
    -w wg: Data type (e.g., wes for exome sequencing, wgs for whole-genome sequencing).
    -o logDir: Directory to store log files.
    -l merlength: Set the MER length (e.g., 21, 27).

Example

bash all_in_one_pipeline.sh -r $runID -s snv-indel-fusion -e netMHCpan4_1 -p patientID_123 -t RNAseq -d debug

This will run the pipeline with the specified steps: preparing the folder, processing HLA data (step s1b), SNV analysis, and converting the results to public XLSX format. Debug mode is enabled, and the output will be logged.

Main Pipeline Steps

    Prepare folder (s0)
    Run HLA typing (s1b)
    SNV, indel, and fusion based neoepitpe analysis
    LOH/CGI feature
    Convert results to XLSX format

License

This project is licensed under the GNU3.0 License - see the LICENSE file for details.
Contact

For questions or issues, please contact the project maintainers or open an issue in this repository.

This README provides a comprehensive guide on how to use the pipeline, the available steps, and how to configure and run the analysis with different parameters.