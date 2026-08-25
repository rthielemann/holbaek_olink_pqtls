library(dplyr)
library(readr)
library(data.table)
library(bignum)    

# Read in dataframe containing genome-wide significant, unclumped SNPs
path <- snakemake@input[[1]]
df <- fread(
  input      = path,
  sep        = "\t",
  header     = TRUE,
  stringsAsFactors = FALSE,
  fill       = TRUE,
  colClasses = list(character = "P-value")   # ← read P-values as text
)

# calculate the lowest p-value in the mhc region
print(df)
mhc_min_p_df <- df %>%
  dplyr::filter(chr == 6 & ps >= 29691116 & ps <= 33054976) %>%
  dplyr::rename(p_val_str = p_wald) %>%            # original string
  dplyr::mutate(p_val = bigfloat(p_val_str)) %>%       # high-precision number
  dplyr::filter(p_val == min(p_val))
mhc_min_p_ <- mhc_min_p_df$p_val
print(mhc_min_p_)

# clamp to 1 if somehow empty/NA
mhc_min_p <- min(mhc_min_p_, bigfloat("1"))
print(mhc_min_p)

# Filter out all MHC SNPs whose p_val > mhc_min_p, then restore original P-value col
df <- df %>%
  dplyr::rename(p_val_str = p_wald) %>%
  dplyr::mutate(p_val = bigfloat(p_val_str),
    neg_log10_p  = -log10(p_val)) %>%
  dplyr::filter(
    !(chr == 6 & ps >= 29691116 & ps <= 33054976 & p_val > mhc_min_p)
  ) %>%
  # put the original text back into a P-value column
  dplyr::mutate(p_wald = p_val_str) %>%
  dplyr::select(-p_val_str, -p_val)

write_tsv(df, snakemake@output[[1]])
