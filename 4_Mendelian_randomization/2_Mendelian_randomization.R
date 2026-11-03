library(TwoSampleMR)
library(magrittr)
library(tidyverse)
library(vroom)
library(data.table)
library(dplyr)
library(purrr)
library(readr)
library(stringr)

wd <- "./"

output_dir <- paste0(wd, "2_output/2_twosampleMR/2_meta_analysis/1_ivw/")

dirs <- c("harmonized", "or", "proxy")
lapply(dirs, function(dir) system(paste0('mkdir -p ', output_dir, dir, '/')))

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
outcome_path <- args[1]
outcomename <- gsub("\\.tsv.gz$", "", basename(outcome_path))

full_exp_gwas <- read_delim("gw_cis_for_ivw_meta.tsv", delim = "\t")

snps_of_interest <- readRDS(paste0(output_dir, "snps_of_interest.rda"))

snps_of_interest_rsID <- snps_of_interest %>%
  dplyr::select(rsID, chr, pos)

# ---- load proxy ----
proxy_df <- read.csv(paste0(output_dir, "proxy/combined_find_proxies.csv"))

proxy_outcomename <- proxy_df %>%
  filter(outcome == outcomename)

# ---- 2SMR ----
# Load list of IV SNPs
iv_list <- snps_of_interest$rsID

# Filter proxies for this outcome
outcome_proxies <- proxy_outcomename %>%
  dplyr::select(rsID, RS_Number)

# Full list of SNPs to extract from exposure: all original IVs, plus any proxies used for this outcome
all_snps_for_exp <- unique(c(iv_list, outcome_proxies$RS_Number))

proxy_outcomename <- proxy_outcomename %>%
  dplyr::mutate(chr_proxy = chr, pos_proxy = position)

# Extract proxy info from protein gwas, get gwas_target for RS_Number
# rsID = pQTL, RS_Number = proxy (need chr_proxy, pos_proxy to extract info from protein gwas)
proxy_outcomename_slim <- proxy_outcomename %>%
  dplyr::select(rsID, chr_proxy, pos_proxy, proxy, R2, RS_Number) %>%
  distinct()

proxy_rsID <- snps_of_interest %>%
  left_join(proxy_outcomename_slim, by = "rsID") %>%
  filter(!is.na(RS_Number)) %>%
  mutate(MarkerName_proxy = paste(chr_proxy, pos_proxy, sep = "_")) %>%
  distinct()

full_exp_gwas_path <- "3_meta_analysis/"

get_proxy_exposure_info <- function(markername, gwas_target, rsid) {
  gwas_protein_file <- list.files(
    full_exp_gwas_path,
    pattern = paste0("^", gwas_target, ".*meta\\.stderr\\.assoc\\.info\\.tbl$"), # start of the filename (gwas_target)
    full.names = TRUE
  )

  if (length(gwas_protein_file) == 0) {
    warning(paste("No assoc found for", gwas_target))
    return(NULL)
  }

  gwas_protein_data <- fread(gwas_protein_file[1])

  proxy_info <- gwas_protein_data %>%
    mutate(
      chr = as.character(chr),
      pos = as.numeric(pos),
      chr_pos = paste(chr, pos, sep = "_"),
      `P-value` = as.numeric(`P-value`)
    ) %>%
    filter(chr_pos %in% markername)

  if (nrow(proxy_info) == 0) {
    message(paste("Proxy", markername, "not found in", gwas_target))
    return(NULL)
  }

  proxy_info$gwas_target <- gwas_target
  proxy_info$rsID <- rsid
  proxy_info$proxy <- TRUE

  return(proxy_info)
}

if (nrow(proxy_rsID) > 0) {
  proxy_exp_list <- purrr::pmap(
    list(
      markername = proxy_rsID$MarkerName_proxy,
      gwas_target = proxy_rsID$gwas_target,
      rsid = proxy_rsID$RS_Number
    ),
    get_proxy_exposure_info
  )

  if (length(proxy_exp_list) > 0 && any(!sapply(proxy_exp_list, is.null))) {
    proxy_exp_data <- proxy_exp_list %>%
      purrr::compact() %>%
      bind_rows() %>%
      mutate(
        RS_Number = rsID,
        chr = as.character(chr),
        pos = as.numeric(pos)
      )
  } else {
    proxy_exp_data <- snps_of_interest[0, ] %>%
      mutate(proxy = logical(), RS_Number = character())
  }
} else {
  proxy_exp_data <- snps_of_interest[0, ] %>%
    mutate(proxy = logical(), RS_Number = character())
}

# Combine original and proxy IVs for 2SMR
snps_of_interest_proxy <- bind_rows(snps_of_interest, proxy_exp_data)

# Extract exposure data for required SNPs
exp_dat <- snps_of_interest_proxy %>%
  dplyr::filter(rsID %in% all_snps_for_exp) %>%
  distinct(rsID, gwas_target, .keep_all = TRUE)

# Format exposure data
exp_dat_fmt <- format_data(exp_dat,
                           phenotype_col = "gwas_target",
                           type = "exposure",
                           snp_col = "rsID",
                           beta_col = "Effect",
                           se_col = "StdErr",
                           effect_allele_col = "Allele1",
                           other_allele_col = "Allele2",
                           eaf_col = "Freq1",
                           pval_col = "P-value"
)

exp_dat_fmt <- exp_dat_fmt %>%
  mutate(snp_type = case_when(
    SNP %in% iv_list ~ "IV",
    SNP %in% snps_of_interest_proxy$RS_Number ~ "proxy",
    TRUE ~ "unknown"
  ))

proxy_outcomename <- proxy_outcomename %>%
  mutate(IV = case_when(
    proxy == TRUE & !is.na(RS_Number) ~ RS_Number,
    proxy == FALSE ~ rsID))

outcome_gwas_filtered <- proxy_outcomename %>%
  filter(IV %in% exp_dat_fmt$SNP) %>%
  mutate(SNP = IV,
         phenotype = outcome)

# Format outcome
outcome_fmt <- format_data(
  as.data.frame(outcome_gwas_filtered),
  type = "outcome",
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  chr_col = "chr",
  pos_col = "position")

# Harmonize
harmonized <- harmonise_data(exp_dat_fmt, outcome_fmt, action = 2) %>%
  filter((eaf.exposure > 0.001 & eaf.exposure < 0.999) &
           (eaf.outcome > 0.001 & eaf.outcome < 0.999))

harmonized <- harmonized %>%
  left_join(dplyr::select(exp_dat_fmt, SNP, id.exposure, snp_type) %>% distinct(),
            by = c("SNP", "id.exposure"))

if (nrow(harmonized) == 0) {
  message("No harmonized SNPs for ", outcomename)
  quit(save = "no")
}

# Save harmonized data
dir.create(paste0(output_dir, "harmonized/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_dir, "or/"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(output_dir, "sensitivity/"), showWarnings = FALSE, recursive = TRUE)

write_tsv(harmonized, paste0(output_dir, "harmonized/", outcomename, ".harmonized.infer_forward_strand.txt"))

# ---- Outcome metadata for Steiger ----
# Columns (no header): outcome name, total N, study type (cc or quant), proportion of cases
outcome_meta <- read_tsv("info_for_coloc.tsv",
                         col_names = c("outcome", "n", "study_type", "prop_cases"),
                         col_types = "cccc") %>%
  mutate(n          = as.numeric(n),
         prop_cases = as.numeric(prop_cases),
         type       = if_else(study_type == "cc", "binary", "continuous"),
         ncase      = if_else(type == "binary", round(n * prop_cases), NA_real_),
         ncontrol   = n - ncase)

meta <- outcome_meta %>%
  filter(outcome == outcomename)
stopifnot(nrow(meta) == 1)
if (meta$type == "binary") stopifnot(!is.na(meta$ncase))

N_EXPOSURE <- 3853   # pQTL meta-analysis sample size

# Use only SNPs retained after harmonisation (the set used by mr())
dat <- harmonized %>%
  filter(mr_keep) %>%
  mutate(samplesize.exposure = N_EXPOSURE,
         samplesize.outcome  = meta$n)

# Variance explained, on the liability scale for binary outcomes
dat$r.exposure <- get_r_from_bsen(dat$beta.exposure, dat$se.exposure, dat$samplesize.exposure)

if (meta$type == "binary") {
  dat$r.outcome <- get_r_from_lor(lor = dat$beta.outcome,
                                  af  = coalesce(dat$eaf.outcome, dat$eaf.exposure),
                                  ncase = meta$ncase, ncontrol = meta$ncontrol,
                                  prevalence = 0.1)
} else {
  dat$r.outcome <- get_r_from_bsen(dat$beta.outcome, dat$se.outcome, dat$samplesize.outcome)
}

# directionality_test() uses r.*, steiger_filtering() uses rsq.* and effective_n.*
dat <- dat %>%
  mutate(rsq.exposure = r.exposure^2,
         rsq.outcome  = r.outcome^2,
         effective_n.exposure = samplesize.exposure,
         effective_n.outcome  = if (meta$type == "binary") {
           2 / (1 / meta$ncase + 1 / meta$ncontrol)
         } else {
           samplesize.outcome
         })

if (any(is.na(dat$r.outcome))) {
  warning(sum(is.na(dat$r.outcome)), " SNPs with missing r.outcome for ", outcomename)
}

# ---- 1. MR estimates (Wald ratio, IVW, MR-Egger, weighted median, weighted mode) ----
mr_results <- mr(dat, method_list = c("mr_wald_ratio", "mr_ivw",
                                      "mr_egger_regression",
                                      "mr_weighted_median", "mr_weighted_mode"))
mr_results_single_snp <- mr_singlesnp(dat) %>%
  filter(str_starts(SNP, "rs"))

# ---- 2. Cochran's Q (2+ SNPs) ----
mr_het <- mr_heterogeneity(dat, method_list = c("mr_ivw", "mr_egger_regression"))

# ---- 3. MR-Egger intercept (3+ SNPs) ----
mr_pleio <- mr_pleiotropy_test(dat)

# ---- 4. MR Steiger (per exposure-outcome pair and per SNP) ----
mr_steiger     <- directionality_test(dat)
mr_steiger_snp <- steiger_filtering(dat) %>%
  dplyr::select(SNP, exposure, id.exposure, rsq.exposure, rsq.outcome, steiger_dir, steiger_pval)

# ---- 5. F-statistics and I2_GX ----
f_stats <- dat %>%
  mutate(F_snp = (beta.exposure / se.exposure)^2) %>%
  group_by(id.exposure, id.outcome, exposure) %>%
  summarise(nsnp     = n(),
            F_mean   = mean(F_snp),
            F_min    = min(F_snp),
            R2_total = sum(r.exposure^2),
            F_overall = (R2_total * (first(samplesize.exposure) - nsnp - 1)) /
                        ((1 - R2_total) * nsnp),
            I2_GX    = if (n() >= 3) Isq(abs(beta.exposure), se.exposure) else NA_real_,
            .groups  = "drop")

# ---- 6. One summary row per exposure-outcome pair ----
res_wide <- mr_results %>%
  mutate(method = recode(method,
                         "Wald ratio" = "Wald",
                         "Inverse variance weighted" = "IVW",
                         "MR Egger" = "Egger",
                         "Weighted median" = "WMedian",
                         "Weighted mode" = "WMode")) %>%
  dplyr::select(id.exposure, id.outcome, method, b, se, pval) %>%
  pivot_wider(names_from = method, values_from = c(b, se, pval),
              names_glue = "{method}_{.value}")

het_wide <- mr_het %>%
  mutate(method = if_else(method == "Inverse variance weighted", "IVW", "Egger")) %>%
  dplyr::select(id.exposure, id.outcome, method, Q, Q_df, Q_pval) %>%
  pivot_wider(names_from = method, values_from = c(Q, Q_df, Q_pval),
              names_glue = "{method}_{.value}")

pleio_slim <- mr_pleio %>%
  dplyr::select(id.exposure, id.outcome,
                Egger_intercept = egger_intercept,
                Egger_intercept_se = se,
                Egger_intercept_pval = pval)

steiger_slim <- mr_steiger %>%
  dplyr::select(id.exposure, id.outcome,
                snp_r2.exposure, snp_r2.outcome,
                correct_causal_direction, steiger_pval)

summary_tbl <- list(f_stats, res_wide, het_wide, pleio_slim, steiger_slim) %>%
  reduce(left_join, by = c("id.exposure", "id.outcome")) %>%
  mutate(outcome = outcomename)

# ---- Write ----
mr_results$outcome <- outcomename
mr_results_single_snp$outcome <- outcomename
write_tsv(mr_results, paste0(output_dir, "or/", outcomename, ".or.infer_forward_strand.txt"))
write_tsv(mr_results_single_snp, paste0(output_dir, "or/", outcomename, ".or.single_snp.infer_forward_strand.txt"))
write_tsv(summary_tbl, paste0(output_dir, "sensitivity/", outcomename, ".sensitivity_summary.txt"))
write_tsv(mr_steiger_snp, paste0(output_dir, "sensitivity/", outcomename, ".steiger_per_snp.txt"))
