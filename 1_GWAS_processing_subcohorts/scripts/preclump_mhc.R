library(dplyr)
library(readr)
library(data.table)

# Read iun dataframe containing genome-wide significant, unclumped SNPs
path <- snakemake@input[[1]]
df <- fread(input = path, sep = "\t", header = T, stringsAsFactors = F, fill=TRUE)

# calculate th elowest p-value in the mhc region
print(df)
mhc_min_p_df <- df %>%
	dplyr::filter((chr = 6 & ps <= 33448354 & ps >= 28477797)) %>%
	dplyr::filter(p_wald == min(p_wald))
mhc_min_p_ <- mhc_min_p_df$p_wald
# This keeps mhc_min_p from being empty, which causes problem in the comparison
mhc_min_p <- min(mhc_min_p_, 1)

# Filter the dataframe to exclude all values in the MHC region which are not the highest p-values
df <- df %>% 
    dplyr::filter(!(chr == 6 & ps <= 33448354 & ps >= 28477797 & p_wald > mhc_min_p)) 

write_tsv(df, snakemake@output[[1]])
