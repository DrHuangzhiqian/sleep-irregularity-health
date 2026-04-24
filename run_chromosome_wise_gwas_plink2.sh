#!/bin/bash

# Accept chromosome number as input
CHR=$1

echo "Processing chromosome ${CHR}"

# Define required input and output paths
PLINK_PATH="/path/to/plink2"
INPUT_BED="/path/to/genotype/ukb_imp_chr${CHR}_v3"
KEEP_FILE="/path/to/keepIDs.txt"
PHENO_FILE="/path/to/pheno_data.txt"
COVAR_FILE="/path/to/covariates.txt"
OUTPUT_PREFIX="/path/to/output/results/chr${CHR}"

# Run GWAS using PLINK2
$PLINK_PATH \
  --bfile $INPUT_BED \
  --keep $KEEP_FILE \
  --pheno $PHENO_FILE \
  --covar $COVAR_FILE \
  --pheno-name SRI_rint \
  --rm-dup force-first \
  --maf 0.01 \
  --geno 0.05 \
  --hwe 1e-6 0 \
  --glm hide-covar cols=chrom,pos,ref,alt,a1freq,nobs,beta,se,tz,p \
  --variance-standardize \
  --no-input-missing-phenotype \
  --vif 1000 \
  --out $OUTPUT_PREFIX
