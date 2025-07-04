#!/usr/bin/env bash

echo 'SCHEME   STDERR' > $3'.metal.txt'
echo 'GENOMICCONTROL OFF' >> $3'.metal.txt'
echo '' >> $3'.metal.txt'
echo 'AVERAGEFREQ ON' >> $3'.metal.txt'
echo 'MINMAXFREQ ON' >> $3'.metal.txt'
echo '' >> $3'.metal.txt'
echo 'CUSTOMVARIABLE    TotalSampleSize' >> $3'.metal.txt'
echo 'LABEL TotalSampleSize as sample_size' >> $3'.metal.txt'
echo 'SEPARATOR TAB' >> $3'.metal.txt'
echo 'MARKER   CHROM:POS:REF:ALT' >> $3'.metal.txt'
echo 'ALLELE   allele1 allele0' >> $3'.metal.txt'
echo 'FREQ     af' >> $3'.metal.txt'
echo 'EFFECT   beta' >> $3'.metal.txt'
echo 'STDERR   se' >> $3'.metal.txt'
echo 'PVAL     p_wald' >> $3'.metal.txt'
echo 'WEIGHTLABEL sample_size' >> $3'.metal.txt'
echo '' >> $3'.metal.txt'
echo '# process POPULATION cohort' >> $3'.metal.txt'
echo 'PROCESS' $1 >> $3'.metal.txt'
echo '' >> $3'.metal.txt'
echo 'SEPARATOR TAB' >> $3'.metal.txt'
echo 'MARKER   CHROM:POS:REF:ALT' >> $3'.metal.txt'
echo 'ALLELE   allele1 allele0' >> $3'.metal.txt'
echo 'FREQ     af' >> $3'.metal.txt'
echo 'EFFECT   beta' >> $3'.metal.txt'
echo 'STDERR   se' >> $3'.metal.txt'
echo 'PVAL     p_wald' >> $3'.metal.txt'
echo 'WEIGHTLABEL sample_size' >> $3'.metal.txt'
echo '' >> $3'.metal.txt'
echo '#process OBERSITY_CLINIC cohort' >> $3'.metal.txt'
echo 'PROCESS' $2 >> $3'.metal.txt'
echo '' >> $3'.metal.txt'
echo 'OUTFILE' $4 .tbl >> $3'.metal.txt'
echo 'ANALYZE HETEROGENEITY' >> $3'.metal.txt'
echo '' >> $3'.metal.txt'
echo 'QUIT' >> $3'.metal.txt'

metal $3'.metal.txt'
