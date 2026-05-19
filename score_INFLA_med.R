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

# --- 并行计算设置 ---
num_cores <- parallel::detectCores() - 12
if(num_cores < 1) num_cores <- 1
registerDoParallel(cores = num_cores)
cat("已注册并行核心数:", num_cores, "\n")

# 定义中介分析函数
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
