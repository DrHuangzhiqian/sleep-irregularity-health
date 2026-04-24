## =============================================================================
## 批量疾病代谢中介分析（Stage 1 正式版 | 防卡死护航版）
## 核心优化：每跑完一个疾病自动重启并行核心，强制释放内存，杜绝跑到一半卡死！
## =============================================================================

BFI_factor <- 313
sims_n <- 1000
num_cores <- 8        # 你可按机器调整
min_n <- 200           # 每个 mediator 完整样本不足则跳过

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(openxlsx)
  library(mediation)
  library(parallel)
  library(foreach)
  library(doParallel)
})

# ==============================================================================
# 1. 路径与目录设定
# ==============================================================================
setwd("D:/ukb data/Sleep regularity/Sleep regularity")

outcome_folder_name <- "16_outcomes"
cox_result_folder   <- "Cox_single_omics_outcomes/Results_Cox_Metabolomics"
result_dir <- "Mediation/single_feature_mediation/Results_Mediation_Metab_260306"
if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

outcome_files_all <- list.files(outcome_folder_name, pattern = "_outcome_.*\\.tsv$")
disease_list <- unique(sapply(strsplit(outcome_files_all, "_outcome_"), `[`, 1))

cat("共检测到", length(disease_list), "个疾病结局需要分析。\n")

# ==============================================================================
# 2. 线性回归筛选
# ==============================================================================
lm_result_file <- "Linear_reg/constant_SRI_linea/SRI_Omics_All_Results_260126.xlsx"
df_lm_res <- read.xlsx(lm_result_file, sheet = "Metabo_Model2")
df_lm_res <- df_lm_res %>% mutate(BFI = P_value * BFI_factor)
lm_sig <- df_lm_res %>% filter(!is.na(BFI), BFI < 0.05) %>% pull(Outcome) %>% unique()

# ==============================================================================
# 5. 读取代谢原始数据并裁剪
# ==============================================================================
Metab_file_path <- "Metabolomic_processed_260113.tsv"
df_Metab <- readr::read_tsv(Metab_file_path, show_col_types = FALSE)
valid_mediators <- intersect(lm_sig, colnames(df_Metab))
df_metab_subset <- df_Metab %>% dplyr::select(eid, all_of(valid_mediators))

df_base <- GGIR_selected %>%
  inner_join(CovariatesImputed, by = "eid") %>%
  inner_join(df_metab_subset, by = "eid") %>%
  as.data.frame()

# （注意：删除了在这里全局启动并行核心的代码，移到了循环内部）

# ==============================================================================
# 7. 开始全疾病大循环
# ==============================================================================
for (disease in disease_list) {
  cat("\n======================================================\n")
  cat(">>> 正在处理疾病:", disease, "\n")
  cat("======================================================\n")
  
  cox_file <- file.path(cox_result_folder, paste0("Result_Cox_Metabolomics_", disease, ".csv"))
  if (!file.exists(cox_file)) next
  
  df_cox <- read.csv(cox_file, stringsAsFactors = FALSE)
  if (!all(c("Feature", "P_Value") %in% colnames(df_cox))) next
  
  df_cox <- df_cox %>% mutate(BFI = P_Value * BFI_factor)
  cox_sig <- df_cox %>% filter(!is.na(BFI), BFI < 0.05) %>% pull(Feature) %>% unique()
  
  target_mediators <- intersect(valid_mediators, cox_sig)
  if (length(target_mediators) == 0) next
  
  outcome_file <- list.files(outcome_folder_name, pattern = paste0("^", disease, "_outcome_.*\\.tsv$"), full.names = TRUE)[1]
  if (is.na(outcome_file)) next
  
  outcome_df <- readr::read_tsv(outcome_file, show_col_types = FALSE) %>%
    dplyr::select(eid, target_status, target_time)
  
  df_analysis0 <- df_base %>%
    inner_join(outcome_df, by = "eid") %>%
    filter(target_time > 0) %>%
    as.data.frame()
  
  exposure_var <- "SRI_scaled"
  covariates_string <- "age_test + MVPA + season + sex + race + tdi + smk + alc + bmi + fastingtime + SleepDurationInSpt_AD_T5A5_mn"
  needed_covs <- c("age_test","MVPA","season","sex","race","tdi","smk","alc","bmi","fastingtime","SleepDurationInSpt_AD_T5A5_mn")
  
  cat("  -> 需跑并行中介代谢物:", length(target_mediators), "个。正在唤醒并行节点...\n")
  
  # 【防御机制 1：在这个疾病开始前，现唤醒并行核心】
  if (num_cores < 1) num_cores <- 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  results_fast <- foreach(
    mediator = target_mediators,
    .packages = c("mediation","stats","dplyr"),
    .combine = rbind,
    .errorhandling = "pass"
  ) %dopar% {
    
    needed_cols <- c("target_status", exposure_var, mediator, needed_covs)
    df_analysis <- df_analysis0[, needed_cols]
    df_analysis <- df_analysis[complete.cases(df_analysis), , drop = FALSE]
    
    if (nrow(df_analysis) < min_n) return(NULL)
    
    out <- tryCatch({
      f_m <- as.formula(paste(mediator, "~", exposure_var, "+", covariates_string))
      model_m <- lm(f_m, data = df_analysis)
      model_m$call$formula <- f_m
      
      f_y <- as.formula(paste("target_status ~", exposure_var, "+", mediator, "+", covariates_string))
      model_y <- glm(f_y, family = binomial(link = "logit"), data = df_analysis)
      model_y$call$formula <- f_y
      
      set.seed(123)
      med_out <- mediate(model_m, model_y, treat = exposure_var, mediator = mediator, boot = TRUE, sims = sims_n)
      s <- summary(med_out)
      
      data.frame(
        Disease = disease, Mediator = mediator, N_Analysis = nrow(df_analysis),
        ACME_Est = s$d0, ACME_P = s$d0.p, ADE_Est = s$z0, ADE_P = s$z0.p,
        Total_Est = s$tau.coef, Total_P = s$tau.p, Prop_Med_Pct = s$n0 * 100, Prop_Med_P = s$n0.p,
        Error = NA_character_
      )
    }, error = function(e) { return(NULL) })
    
    out
  }
  
  # 【防御机制 2：该疾病跑完立刻强制关闭并行节点，释放全部垃圾内存！】
  stopCluster(cl)
  
  # 【防御机制 3：主进程主动召唤一次垃圾回收，清空无用缓存】
  rm(df_analysis0, outcome_df)
  gc() 
  
  if (!is.null(results_fast) && nrow(results_fast) > 0) {
    save_path <- file.path(result_dir, paste0("Mediation_Metab_", disease, "_sims", sims_n, "_", Sys.Date(), ".csv"))
    write.csv(results_fast, save_path, row.names = FALSE)
    cat("  ✅ 保存成功：", save_path, "\n")
  } else {
    cat("  ⚠️ 无有效结果生成。\n")
  }
} 

cat("\n🎯 全部任务圆满结束！内存安全控制发挥作用。\n")