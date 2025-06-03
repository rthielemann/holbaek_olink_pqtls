library(dplyr)
library(readr)
library(data.table)

path <- snakemake@input[[1]]
df <- fread(input = path, sep = "\t", header = T, stringsAsFactors = F, fill=TRUE)

df <- df %>% 
    dplyr::filter(tss_chr != 'NA') %>%
    dplyr::mutate(across(c(pos, tss_position), as.numeric)) %>%
    dplyr::mutate(across(c(chr, tss_chr), as.character)) %>%
    dplyr::mutate(cis_trans = case_when(
        (chr == tss_chr & pos > tss_position - 1000000 & pos < tss_position + 1000000) ~ 'cis',
         TRUE ~ 'trans'))


    
write_tsv(df, snakemake@output[[1]])
