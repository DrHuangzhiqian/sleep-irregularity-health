## INFLA-score的计算与中介分析

# ==============================================================================
# 0. 环境设置与包加载
# ==============================================================================
library(data.table)
library(dplyr)
library(mediation) 
library(foreach)
library(doParallel)
library(stringr)

# --- 设置工作目录 ---
# 请修改为你实际的路径
setwd("D:/ukb data/Sleep regularity/Sleep regularity")

# --- 并行计算设置 ---
num_cores <- parallel::detectCores() - 12
if(num_cores < 1) num_cores <- 1
registerDoParallel(cores = num_cores)
cat("已注册并行核心数:", num_cores, "\n")

# ==============================================================================
# 1. 数据准备：INFLA-score 计算
# ==============================================================================
cat("正在处理炎症数据并计算 INFLA-score...\n")

# 1.1 读取炎症数据
inflam_data <- fread("Blood_inflam_processed_260115.tsv") 

# 1.2 计算 INFLA-score
# 逻辑：将4个指标分为10个等级(Deciles)，根据特定规则赋分求和
bloodcount_index <- inflam_data %>%
  dplyr::select(eid, White_blood_cell, Platelet, NLR, CRP) %>%
  na.omit() 

INFLA_score_df <- bloodcount_index %>%
  mutate(
    # 计算分位数
    Decile_White_blood_cell = ntile(White_blood_cell, 10),
    Decile_Platelet = ntile(Platelet, 10),
    Decile_NLR = ntile(NLR, 10),
    Decile_CRP = ntile(CRP, 10),
    
    # 赋分
    Score_White_blood_cell = case_when(
      Decile_White_blood_cell %in% 1:4 ~ Decile_White_blood_cell - 5,
      Decile_White_blood_cell %in% 5:6 ~ 0,
      Decile_White_blood_cell %in% 7:10 ~ Decile_White_blood_cell - 6,
      TRUE ~ 0
    ),
    Score_Platelet = case_when(
      Decile_Platelet %in% 1:4 ~ Decile_Platelet - 5,
      Decile_Platelet %in% 5:6 ~ 0,
      Decile_Platelet %in% 7:10 ~ Decile_Platelet - 6,
      TRUE ~ 0
    ),
    Score_NLR = case_when(
      Decile_NLR %in% 1:4 ~ Decile_NLR - 5,
      Decile_NLR %in% 5:6 ~ 0,
      Decile_NLR %in% 7:10 ~ Decile_NLR - 6,
      TRUE ~ 0
    ),
    Score_CRP = case_when(
      Decile_CRP %in% 1:4 ~ Decile_CRP - 5,
      Decile_CRP %in% 5:6 ~ 0,
      Decile_CRP %in% 7:10 ~ Decile_CRP - 6,
      TRUE ~ 0
    ),
    # 总分
    INFLA = Score_White_blood_cell + Score_Platelet + Score_NLR + Score_CRP
  ) %>%
  dplyr::select(eid, INFLA)

# --- 保存 INFLA_score 数据 ---
if (!dir.exists("Blood_Inflam")) dir.create("Blood_Inflam")
fwrite(INFLA_score_df, "Blood_Inflam/INFLA_score_calculated.csv")
cat(">>> INFLA score 数据已保存至: Blood_Inflam/INFLA_score_calculated.csv\n")
# ------------------------------------

cat("INFLA-score 计算完成。有效样本量:", nrow(INFLA_score_df), "\n")

# ==============================================================================
# 2. 数据准备：整合协变量与暴露变量
# ==============================================================================
cat("正在导入协变量 (含睡眠时长) 和 暴露变量 (SRI)...\n")

# 2.1 导入 GGIR
GGIR_selected <- fread("GGIR_selected.csv") %>%
  dplyr::select(eid, SleepRegularityIndex_AD_T5A5_mn, age_test, MVPA, season, SleepDurationInSpt_AD_T5A5_mn) %>%
  mutate(season = as.factor(season)) %>%
  # 标准化 SRI
  mutate(SRI_scaled = as.numeric(scale(SleepRegularityIndex_AD_T5A5_mn))) %>%
  dplyr::select(-SleepRegularityIndex_AD_T5A5_mn) 

# 2.2 导入其他协变量
CovariatesImputed <- read.csv("CovariatesImputed.csv", stringsAsFactors = FALSE) %>%
  dplyr::select(eid, sex, race, tdi, smk, alc, bmi, fastingtime) %>%
  mutate(across(c(sex, race, smk, alc), as.factor)) %>%
  mutate(across(c(tdi, bmi, fastingtime), as.numeric))

# 2.3 合并所有基线数据 (X + M + Covariates)
Base_Pheno <- inner_join(GGIR_selected, CovariatesImputed, by = "eid")
Mediation_Base_Data <- inner_join(Base_Pheno, INFLA_score_df, by = "eid")

cat("基线数据合并完成 (X+M+Cov)。总样本量:", nrow(Mediation_Base_Data), "\n")

# ==============================================================================
# 3. 定义中介分析函数
# ==============================================================================
run_mediation_analysis <- function(outcome_name, outcome_data, base_data) {
  
  # 合并结局 (只保留 target_status)
  analysis_df <- inner_join(base_data, outcome_data, by = "eid") %>%
    na.omit() 
  
  if (nrow(analysis_df) < 500) return(NULL)
  
  # --- Model M: INFLA ~ SRI + Covariates ---
  f_m <- as.formula("INFLA ~ SRI_scaled + age_test + MVPA + season + sex + race + tdi + smk + alc + bmi + fastingtime + SleepDurationInSpt_AD_T5A5_mn")
  model_m <- lm(f_m, data = analysis_df)
  
  # 【关键修正 1】强制将公式对象植入 call 中，防止 "object 'f_m' not found"
  model_m$call$formula <- f_m
  
  # --- Model Y: Outcome ~ SRI + INFLA + Covariates ---
  f_y <- as.formula("target_status ~ SRI_scaled + INFLA + age_test + MVPA + season + sex + race + tdi + smk + alc + bmi + fastingtime + SleepDurationInSpt_AD_T5A5_mn")
  model_y <- glm(f_y, family = binomial(link = "logit"), data = analysis_df)
  
  # 【关键修正 2】同上，强制植入公式
  model_y$call$formula <- f_y
  
  # --- 运行中介 (Quasi-Bayesian) ---
  set.seed(123)
  med_out <- tryCatch({
    mediate(model.m = model_m, 
            model.y = model_y, 
            treat = "SRI_scaled", 
            mediator = "INFLA", 
            boot = TRUE,
            sims = 1000) 
  }, error = function(e) return(NULL))
  
  if (is.null(med_out)) return(NULL)
  
  # 保存 RDS
  rds_file_path <- file.path("Blood_Inflam/INFLA_score_mediation", paste0("med_INFLA_score_", outcome_name, ".rds"))
  saveRDS(med_out, file = rds_file_path)
  
  res_summary <- summary(med_out)
  
  # 提取结果
  res <- data.frame(
    Outcome = outcome_name,
    N = nrow(analysis_df),
    
    ACME_Est = res_summary$d0,
    ACME_P   = res_summary$d0.p,
    ACME_CI_Low = res_summary$d0.ci[1],
    ACME_CI_Up  = res_summary$d0.ci[2],
    
    ADE_Est  = res_summary$z0,
    ADE_P    = res_summary$z0.p,
    ADE_CI_Low  = res_summary$z0.ci[1],
    ADE_CI_High = res_summary$z0.ci[2],
    
    Total_Est = res_summary$tau.coef,
    Total_P   = res_summary$tau.p,
    Total_Low = res_summary$tau.ci[1],
    Total_High = res_summary$tau.ci[2],
    
    Prop_Med_Pct = res_summary$n0 * 100,
    Prop_Med_P   = res_summary$n0.p,
    Prop_Med_Low = res_summary$n0.ci[1],
    Prop_Med_High = res_summary$n0.ci[2]
  )
  return(res)
}

# ==============================================================================
# 4. 主循环
# ==============================================================================

# 创建多级输出目录 (recursive = TRUE 确保父文件夹不存在时也会创建)
output_folder <- "Blood_Inflam/INFLA_score_mediation"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

cat("输出目录已创建:", output_folder, "\n")
cat("待分析结局:", length(outcome_files), "\n")

# --- 并行执行 ---
start_time_all <- Sys.time()

results_list <- foreach(file_path = outcome_files, .combine = rbind, .packages = c("dplyr", "mediation", "data.table")) %dopar% {
  
  # 获取结局名称
  clean_name <- tools::file_path_sans_ext(basename(file_path))
  clean_name <- gsub("_outcome_260116", "", clean_name) 
  
  # 读取结局
  outcome_df <- fread(file_path, select = c("eid", "target_status"))
  
  # 运行分析
  res <- run_mediation_analysis(clean_name, outcome_df, Mediation_Base_Data)
  
  return(res)
}

end_time_all <- Sys.time()
cat("分析完成！总耗时:", round(difftime(end_time_all, start_time_all, units="mins"), 2), "分钟。\n")

# ==============================================================================
# 5. 汇总保存
# ==============================================================================
if (!is.null(results_list)) {
  final_file <- file.path(output_folder, "Total_Mediation_INFLA_Summary.csv")
  write.csv(results_list, final_file, row.names = FALSE)
  cat("汇总结果已保存至:", final_file, "\n")
} else {
  cat("无结果生成。\n")
}

stopImplicitCluster()