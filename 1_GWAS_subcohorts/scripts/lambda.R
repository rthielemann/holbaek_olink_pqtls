library(dplyr)
library(data.table)
library(qqman)

sum_stats <- snakemake@input[[1]]
print(sum_stats)

gwasResults <- fread(input = sum_stats, sep = "\t", header = T, stringsAsFactors = F)

print(head(gwasResults))

gwasResults <- gwasResults %>% 
    dplyr::rename(SNP = 'CHROM:POS:REF:ALT', BP = ps, P = p_wald, CHR = 'chr', ) %>% # rename for manhattan package
    dplyr::select(SNP, CHR, BP, P) %>%
    dplyr::mutate(across(c(CHR, BP, P), as.numeric))
    
## QQ plot calculations
p <- gwasResults$P[!is.na(gwasResults$P)]
n <- length(p)
x2obs <- qchisq(p,1,lower.tail=FALSE)
x2exp <- qchisq((1:n - 0.5)/n,1,lower.tail=FALSE)
lambda <- median(x2obs)/median(x2exp) # calculates your lambda value to check for inflation


print(snakemake@wildcards[["sample"]])
print(snakemake@output[["lambDa"]])
write(paste(c(snakemake@wildcards[["sample"]], lambda), collapse = "\t"), file=snakemake@output[["lambDa"]])
