# holbaek_olink_pqtls

- The preprint is on medRxiv: https://www.medrxiv.org/content/10.1101/2025.03.18.25324169v2.
- GWAS summary statistics were deposited on GWAS catalogue under the accessions GCST90589888–GCST90590241.
- This repository provides several custom scripts used in this project.

## Contents
Folder                      | Description
------------------------- | --------------------------------------
0_GEMMA | Simple example script of association testing with GEMMA.
1_GWAS_processing_subcohorts |  Snakemake pipeline to process summary statistics in population-based cohort (n = 2,095) and obesity clinic cohort (n = 1,758).
2_GWAS_meta_analysis | Snakemake pipeline to meta-analyse summary statistics from the population-based and obesity clinic cohorts (n=3,853). Annotation of meta-analysed summary statistics.
3_Stratified_pQTL_analysis | R markdown script to perform obesity-stratified comparison of lead cis-pQTLs.
4_Mendelian_randomization | R markdown script to perform one-sample Mendelian randomization in the population cohort.
other_scripts | Other relevant scripts (Olink proteomics QC, TSS extraction from Biomart)

## Study design
![pQTL study overview](https://github.com/user-attachments/assets/f1b20bf4-e0bc-4cc0-8c69-9084d53adb2f)
