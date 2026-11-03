# holbaek_olink_pqtls

- The preprint is on medRxiv: https://www.medrxiv.org/content/10.1101/2025.03.18.25324169v2.
- GWAS summary statistics will be uploaded to GWAS catalogue upon publication.
- This repository provides several custom scripts used in this project.

## Contents
Folder                      | Description
------------------------- | --------------------------------------
0_GEMMA | Simple example script of association testing with GEMMA.
1_GWAS_processing_subcohorts |  Snakemake pipeline to process summary statistics in population-based cohort (n = 2,095) and obesity clinic cohort (n = 1,758).
2_GWAS_meta_analysis | Snakemake pipeline to meta-analyse summary statistics from the population-based and obesity clinic cohorts (n=3,853). Annotation of meta-analysed summary statistics.
3_Stratified_pQTL_analysis | R markdown script to perform obesity-stratified comparison of pQTLs.
4_Mendelian_randomization | R scripts to perform proxy-identification, two-sample Mendelian randomization, and colocalisation between pQTLs and cardiometabolic outcomes.
other_scripts | Other relevant scripts (Olink proteomics QC, TSS extraction from Biomart)

## Other resources
- The pipeline for colocalisation analysis with GTEx eQTLs can be found here: https://github.com/Hansen-Group/GTEX_eqtl_colocalization

## Study design
<img width="2869" height="3328" alt="pQTL study overview" src="https://github.com/user-attachments/assets/44f7b03f-a1f8-4415-a3eb-089b30b92183" />


