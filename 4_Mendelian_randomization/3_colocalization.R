args <- commandArgs(trailingOnly = TRUE)
trait_name <- args[1]
sample_size_trait <- as.numeric(args[2])
study_type <- args[3]
cc_ratio <- ifelse(args[4] == "NA", NA, as.numeric(args[4]))

# Log standard output
log_file <- paste0("log/", trait_name, ".coloc.txt")
log_conn <- file(log_file, open = "wt")
sink(log_conn)
sink(log_conn, type = "message")

library(tidyverse)
library(coloc)
library(TwoSampleMR)

## Read in cardio sumstat
out_path <- "1_harmonized_sumstats/"
cardio_full_sumstat <- read.table(paste0(out_path, trait_name, ".tsv.gz"), sep = "\t", header = TRUE)

## Read in pqtl information
lead_cis_pqtls <- read_tsv("gw_cis_for_MR.tsv")
lead_cis_n <- lead_cis_pqtls %>%
  dplyr::rename(protein = protein_name)

coloc_results_df <- data.frame(
  protein = character(),
  trait = character(),
  nsnps = numeric(),
  pp_h4_abf = numeric(),
  stringsAsFactors = FALSE
)

for (protein_identifier in lead_cis_n$protein) {
  chr_lead = lead_cis_n %>% dplyr::filter(protein == protein_identifier) %>% pull(chr)
  pos_lead = lead_cis_n %>% dplyr::filter(protein == protein_identifier) %>% pull(pos)
  gwas_target = lead_cis_n %>% dplyr::filter(protein == protein_identifier) %>% pull(gwas_target)
  pqtl_sumstat <- read_tsv(paste0("1_pQTL_regions/", gwas_target, ".meta.assoc.extracted.txt"), show_col_types = FALSE)
  print(protein_identifier)

  ## select correct region from cardio sumstat
  cardio_sumstat <- cardio_full_sumstat %>%
  dplyr::filter(chromosome == chr_lead, abs(base_pair_location - pos_lead) < 500000) %>%
  dplyr::mutate(chromosome = as.numeric(chromosome))

  ## Match cardiometabolic trait alleles with pQTL alleles using TwoSampleMR.
  ## This approach handles strand flips, matches on chr:pos, and removes duplicate chr:pos entries.
  ## 1. Format exposure
   exp_dat_fmt <- format_data(pqtl_sumstat %>% dplyr::mutate(protein = protein_identifier,
                                                            chr_pos = paste0(chromosome, "_", base_pair_location)) %>%
                               dplyr::group_by(chr_pos) %>%
                               dplyr::filter(dplyr::n() == 1) %>%
                               dplyr::ungroup(),
                           phenotype_col = trait_name,
                           type = "exposure",
                           snp_col = "chr_pos",
                           beta_col = "beta",
                           se_col = "standard_error",
                           effect_allele_col = "effect_allele",
                           other_allele_col = "other_allele",
                           eaf_col = "effect_allele_frequency",
                           pval_col = "p_value"
)
  ## 2. Format outcome
  outcome_fmt <- format_data(cardio_sumstat %>% dplyr::mutate(trait = trait_name,
                                                              chr_pos = paste0(chromosome, "_", base_pair_location)) %>%
                               dplyr::group_by(chr_pos) %>%
                               dplyr::filter(dplyr::n() == 1) %>%
                               dplyr::ungroup(),
                  type = "outcome",
                  phenotype_col = trait_name,
                  snp_col = "chr_pos",
                  beta_col = "beta",
                  se_col = "standard_error",
                  eaf_col = "effect_allele_frequency",
                  effect_allele_col = "effect_allele",
                  other_allele_col = "other_allele",
                  pval_col = "p_value",
                  chr_col = "chr",
                  pos_col = "position"
      )

  harmonized <- TwoSampleMR::harmonise_data(exp_dat_fmt, outcome_fmt, action = 2) %>%
  dplyr::filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) &
           (eaf.outcome > 0.001 & eaf.outcome < 0.999))

  # Only apply filter when there are >0 harmonized variants
  if (nrow(harmonized) > 0) {harmonized <- harmonized %>% dplyr::filter(mr_keep == TRUE)}

  if (nrow(harmonized) == 0) {
    coloc_results_df <- dplyr::bind_rows(coloc_results_df, tibble::tibble(
      protein = protein_identifier,
      trait = trait_name,
      nsnps = 0,
      pp_h4_abf = NA_real_))
  next
  }

  ## pQTLs: Create columns necessary for coloc
  pqtl_coloc <- list(
    beta = harmonized$beta.exposure,
    varbeta = harmonized$se.exposure^2,
    MAF = ifelse(harmonized$eaf.exposure > 0.5, 1-harmonized$eaf.exposure, harmonized$eaf.exposure),         # optional, can be omitted if unavailable
    snp = paste0(harmonized$SNP,"_", harmonized$other_allele.exposure, "_", harmonized$effect_allele.exposure),
    pvalues = harmonized$pval.exposure,
    position = sub(".*_", "", harmonized$SNP),
    N = 3853,
    type = "quant"
  )

  cardio_coloc <- list(
    beta = harmonized$beta.outcome,
    varbeta = harmonized$se.outcome^2,
    MAF = ifelse(harmonized$eaf.outcome > 0.5, 1-harmonized$eaf.outcome, harmonized$eaf.outcome),         # optional, can be omitted if unavailable
    snp = paste0(harmonized$SNP,"_", harmonized$other_allele.outcome, "_", harmonized$effect_allele.outcome),
    pvalues = harmonized$pval.outcome,
    position = sub(".*_", "", harmonized$SNP),
    N = sample_size_trait,
    type = study_type
  )

  if (!is.na(cc_ratio) && study_type == "cc") {
    cardio_coloc$s = cc_ratio
  }

  coloc_results <- coloc::coloc.abf(pqtl_coloc, cardio_coloc)

  coloc_results_df <- dplyr::bind_rows(coloc_results_df, tibble::tibble(
    protein = protein_identifier,
    trait = trait_name,
    nsnps = coloc_results$summary["nsnps"],
    pp_h0_abf = coloc_results$summary["PP.H0.abf"],
    pp_h1_abf = coloc_results$summary["PP.H1.abf"],
    pp_h2_abf = coloc_results$summary["PP.H2.abf"],
    pp_h3_abf = coloc_results$summary["PP.H3.abf"],
    pp_h4_abf = coloc_results$summary["PP.H4.abf"]
  ))
}
write_tsv(coloc_results_df, paste0("3_coloc_output_twosamplemr/", trait_name, ".coloc.txt"))

sink(type = "message")
sink()
close(log_conn)