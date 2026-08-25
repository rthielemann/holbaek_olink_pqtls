library(dplyr)
library(readr)
library(data.table)
library(stringr)
library(LDlinkR)

wd <- "./"

output_dir <- paste0(wd, "2_output/2_twosampleMR/2_meta_analysis/1_ivw/")

dirs <- c("proxy")
lapply(dirs, function(dir) system(paste0('mkdir -p ', output_dir, dir, '/')))

# Read arguments
args <- commandArgs(trailingOnly = TRUE)
outcome_path <- args[1]
gwas_files <- list.files(outcome_path, pattern = "\\.tsv\\.gz$", full.names = TRUE)

full_exp_gwas <- read_delim("gw_cis_for_ivw_meta.tsv", delim = "\t")

# ---- find proxy per outcome ----

# Load the .bim file to get chr and pos for rsID
bim <- fread("holbaek.bim", col.names = c("chr", "rsid", "cm", "pos", "allele1", "allele2"))

bim <- bim %>%
  mutate(MarkerName = paste0(chr, ":", pos, ":", allele2, ":", allele1))

# Load SNPs of interest
snps_of_interest <- full_exp_gwas %>%
  left_join(bim %>% dplyr::select(MarkerName, rsid), by = "MarkerName")

snps_of_interest <- snps_of_interest %>%
  mutate(rsID = if_else(
    rsid == "9:90347450:G:T", "rs3128511", rsid)) %>%
  mutate(chr = as.character(chr), pos = as.numeric(pos))

saveRDS(snps_of_interest, paste0(output_dir, "snps_of_interest.rda"))

snps_of_interest_rsID <- snps_of_interest %>%
  dplyr::select(rsID, chr, pos)

## Load example GWAS to test if proxy is in there
holbaek_example_gwas <- read_tsv("9_nobmi_OPG.population.assoc.info.filtered.sample_size") %>%
  dplyr::select('CHROM:POS:REF:ALT', rs)

ld_token <- Sys.getenv("LDLINK_TOKEN")

# ---- find proxy ----
find_proxies <- function(lead_snp, ld_token, r2_threshold = 0.8) {
  tryCatch({
    SNP_LD <- LDproxy(
      snp = lead_snp, pop = "EUR", r2d = "r2",
      token = ld_token, genome_build = "grch37"
    )

    required_cols <- c("RS_Number", "Coord", "R2")
    if (!all(required_cols %in% colnames(SNP_LD))) {
      warning(paste("LDproxy output for", lead_snp, "missing expected columns"))
      return(data.frame(RS_Number = character(), CHROM = character(), POS = numeric(), R2 = numeric()))
    }

    SNP_LD %>%
      dplyr::filter(R2 > r2_threshold) %>%
      dplyr::mutate(
        CHROM = str_extract(Coord, "chr[^:]*"),
        POS = as.numeric(str_extract(Coord, "(?<=:)[0-9]*")),
        CHROM = gsub("chr", "", CHROM),
        CHROM = as.character(CHROM),
        POS = as.numeric(POS)
      ) %>%
      dplyr::select(RS_Number, CHROM, POS, R2)
  }, error = function(e) {
    message(paste("Skipping", lead_snp, ":", e$message))
    return(data.frame(RS_Number = character(), CHROM = character(), POS = numeric(), R2 = numeric()))
  })
}

combined_results <- list()
skipped_snps <- list()

for (gwas_file in gwas_files) {
  outcomename <- gsub("\\.tsv.gz$", "", basename(gwas_file))
  message("Processing file: ", gwas_file)

  gwas_data <- readr::read_tsv(gwas_file, show_col_types = FALSE) %>%
    dplyr::rename(position = base_pair_location, chr = chromosome) %>%
    dplyr::mutate(
      chr = as.character(chr),
      position = as.numeric(position)
    )

  snp_results <- list()

  for (i in 1:nrow(snps_of_interest_rsID)) {
    snp <- snps_of_interest_rsID[i, ]

    snp_found <- gwas_data %>%
      dplyr::filter(
        chr == as.character(snp$chr),
        position == as.numeric(snp$pos)
      )

    if (nrow(snp_found) > 0) {
      snp_results[[i]] <- snp_found %>%
        dplyr::mutate(
          chr = as.character(chr),
          position = as.numeric(position),
          rsID = snp$rsID,
          proxy = FALSE,
          R2 = NA_real_,
          RS_Number = NA_character_
        )
    } else {
      proxies <- find_proxies(snp$rsID, ld_token) %>%
        dplyr::mutate(
          CHROM = as.character(CHROM),
          POS = as.numeric(POS)
        )

      if (nrow(proxies) > 0) {
        proxy_found <- gwas_data %>%
          ## filters if proxy in outcome gwas, filter for proxy in exposure gwas is below
          dplyr::filter(chr %in% proxies$CHROM, position %in% proxies$POS) %>%
          dplyr::mutate(
            chr = as.character(chr),
            position = as.numeric(position)
          )

        if (nrow(proxy_found) > 0) {
          best_proxy <- proxy_found %>%
            dplyr::left_join(proxies, by = c("chr" = "CHROM", "position" = "POS")) %>%
            ## Only keep rsids also present in the reference cohort
            dplyr::filter(RS_Number %in% holbaek_example_gwas$rs)

          if (nrow(best_proxy) > 0) {
            best_proxy <- best_proxy %>%
              dplyr::arrange(desc(R2)) %>%
              dplyr::slice(1) %>%
              dplyr::mutate(
                chr = as.character(chr),
                position = as.numeric(position),
                rsID = snp$rsID,
                proxy = TRUE
              )

            snp_results[[i]] <- best_proxy
          } else {
            snp_results[[i]] <- data.frame(
              chr = as.character(snp$chr),
              position = as.numeric(snp$pos),
              rsID = snp$rsID,
              proxy = TRUE,
              R2 = NA_real_,
              note = "Proxy not in exposure GWAS"
            )
          }
        } else {
          snp_results[[i]] <- data.frame(
            chr = as.character(snp$chr),
            position = as.numeric(snp$pos),
            rsID = snp$rsID,
            proxy = TRUE,
            R2 = NA_real_,
            note = "Proxy not in outcome GWAS"
          )
        }
      } else {
        snp_results[[i]] <- data.frame(
          chr = as.character(snp$chr),
          position = as.numeric(snp$pos),
          rsID = snp$rsID,
          proxy = TRUE,
          R2 = NA_real_,
          note = "No proxy returned"
        )
      }
    }
  }

  outcome_results <- dplyr::bind_rows(snp_results) %>%
    dplyr::mutate(
      chr = as.character(chr),
      position = as.numeric(position),
      outcome = outcomename
    )

  combined_results[[gwas_file]] <- outcome_results
}

proxy_df <- dplyr::bind_rows(combined_results)

write.csv(proxy_df, paste0(output_dir, "proxy/", "combined_find_proxies.csv"), row.names = FALSE)