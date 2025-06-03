# Association analysis were run with GEMMA 0.98.3 
# Input files:
# - holbaek.genotypes (Genotype dosages in the BIMBAM mean genotype format)
# - phenotypes.csv (Phenotype file with individuals as rows and one protein per column)
# - holbaek.kinship.cXX.txt ( Relatedness matrix generated with GEMMA 0.98.3) 
# - covariates.csv (Covariate file with individuals as rows and covariates as columns)

gemma \
-g holbaek.genotypes \
-p phenotypes.csv \
-k holbaek.kinship.cXX.txt \
-c covariates.csv \
-lmm 1 \
-hwe 0.000001 \
-o ${SLURM_ARRAY_TASK_ID}'_'${protein}'' \
-n ${SLURM_ARRAY_TASK_ID}
