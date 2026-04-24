library(data.table)
library(dplyr)
library(TwoSampleMR)
library(ggplot2)

source("/path/to/ld_clump.R")

gwas_all <- fread("/path/to/GWAS_summary_statistics.csv")

gwas_all <- gwas_all %>%
  mutate(
    SNP = ID,
    beta = BETA,
    se = SE,
    pval = P,
    effect_allele = A1,
    eaf = A1_FREQ,
    samplesize = OBS_CT,
    other_allele = ifelse(A1 == ALT, REF, ALT)
  ) %>%
  filter(
    !is.na(SNP),
    !is.na(beta),
    !is.na(se),
    !is.na(pval),
    !is.na(effect_allele),
    !is.na(other_allele),
    !is.na(eaf),
    !is.na(samplesize),
    pval > 0,
    pval <= 1,
    CHR %in% 1:22
  )

expo_file <- gwas_all %>%
  transmute(
    SNP = SNP,
    beta = beta,
    se = se,
    A1 = effect_allele,
    A2 = other_allele,
    p = pval
  )

expo_rt <- format_data(
  expo_file,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "p"
)

expo_rt$samplesize.exposure <- XXX
expo_rt$exposure <- "XXX"

expo_rt <- expo_rt[expo_rt$pval.exposure < XXX, ]

expo_rt$id <- expo_rt$id.exposure
expo_rt$rsid <- expo_rt$SNP
expo_rt$pval <- expo_rt$pval.exposure

expo_rt <- ld_clump_local(
  dat = expo_rt,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p = 1,
  bfile = "/path/to/reference_panel",
  plink_bin = "/path/to/plink"
)

expo_rt$id <- NULL
expo_rt$rsid <- NULL
expo_rt$pval <- NULL

outc_rt <- read_outcome_data(
  snps = expo_rt$SNP,
  filename = "/path/to/outcome_GWAS_summary_statistics.gz",
  sep = "\t",
  snp_col = "XXX",
  beta_col = "XXX",
  se_col = "XXX",
  effect_allele_col = "XXX",
  other_allele_col = "XXX",
  pval_col = "XXX"
)

harm_rt <- harmonise_data(
  exposure_dat = expo_rt,
  outcome_dat = outc_rt,
  action = 2
)

harm_rt <- steiger_filtering(harm_rt)
harm_rt <- subset(harm_rt, mr_keep == TRUE)

harm_rt$R2 <- (2 * harm_rt$beta.exposure^2) /
  (2 * harm_rt$beta.exposure^2 +
     2 * harm_rt$samplesize.exposure * harm_rt$se.exposure^2)

harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
harm_rt$meanf <- mean(harm_rt$f, na.rm = TRUE)
harm_rt <- harm_rt[harm_rt$f > 10, ]

outcome <- "XXX"
out_dir <- file.path("/path/to/output_directory", outcome)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

fwrite(
  harm_rt,
  file = file.path(out_dir, paste0("MR_harmonised_XXX_", outcome, ".txt")),
  sep = "\t"
)

mr_result <- mr(harm_rt)
OR <- generate_odds_ratios(mr_result)

write.csv(
  mr_result,
  file = file.path(out_dir, paste0("MR_result_raw_", outcome, ".csv")),
  row.names = FALSE
)

write.csv(
  OR,
  file = file.path(out_dir, paste0("MR_result_OR_", outcome, ".csv")),
  row.names = FALSE
)

heterogeneity_res <- mr_heterogeneity(harm_rt)
pleiotropy_res <- mr_pleiotropy_test(harm_rt)
singlesnp_res <- mr_singlesnp(harm_rt)
leaveoneout_res <- mr_leaveoneout(harm_rt)

write.csv(
  heterogeneity_res,
  file = file.path(out_dir, paste0("MR_heterogeneity_", outcome, ".csv")),
  row.names = FALSE
)

write.csv(
  pleiotropy_res,
  file = file.path(out_dir, paste0("MR_pleiotropy_", outcome, ".csv")),
  row.names = FALSE
)

write.csv(
  singlesnp_res,
  file = file.path(out_dir, paste0("MR_singlesnp_", outcome, ".csv")),
  row.names = FALSE
)

write.csv(
  leaveoneout_res,
  file = file.path(out_dir, paste0("MR_leaveoneout_", outcome, ".csv")),
  row.names = FALSE
)

p1 <- mr_scatter_plot(mr_result, harm_rt)

ggsave(
  filename = file.path(out_dir, paste0("MR_scatter_plot_", outcome, ".pdf")),
  plot = p1[[1]],
  width = 6,
  height = 6
)

p2 <- mr_forest_plot(singlesnp_res)

ggsave(
  filename = file.path(out_dir, paste0("MR_forest_plot_", outcome, ".pdf")),
  plot = p2[[1]],
  width = 7,
  height = 10
)

p3 <- mr_leaveoneout_plot(leaveoneout_res)

ggsave(
  filename = file.path(out_dir, paste0("MR_leaveoneout_plot_", outcome, ".pdf")),
  plot = p3[[1]],
  width = 7,
  height = 10
)

p4 <- mr_funnel_plot(singlesnp_res)

ggsave(
  filename = file.path(out_dir, paste0("MR_funnel_plot_", outcome, ".pdf")),
  plot = p4[[1]],
  width = 7,
  height = 7
)