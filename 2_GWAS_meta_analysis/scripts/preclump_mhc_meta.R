library(dplyr)
library(readr)
library(data.table)

# Read iun dataframe containing genome-wide significant, unclumped SNPs
path <- snakemake@input[[1]]
df <- fread(input = path, sep = "\t", header = T, stringsAsFactors = F, fill=TRUE)

# calculate the lowest p-value in the mhc region
print(df)
mhc_min_p_df <- df %>%
	dplyr::filter((chr = 6 & pos <= 33448354 & pos >= 28477797)) %>%
	# rename to skip problem with - in name
	dplyr::rename(p_val = 'P-value') %>% 
	dplyr::filter(p_val == min(p_val))
mhc_min_p_ <- mhc_min_p_df$p_val
print(mhc_min_p_)
# This keeps mhc_min_p from being empty, which causes problem in the comparison
mhc_min_p <- min(mhc_min_p_, 1)
print(mhc_min_p)
# Filter the dataframe to exclude all values in the MHC region which are not the highest p-values
df <- df %>% 
	dplyr::rename(p_val = 'P-value') %>%
	dplyr::filter(!(chr == 6 & pos <= 33448354 & pos >= 28477797 & p_val > mhc_min_p)) %>%
	dplyr::rename('P-value' = p_val)
	#reverse renaming to not break later steps in pipeline
print(df)
write_tsv(df, snakemake@output[[1]])