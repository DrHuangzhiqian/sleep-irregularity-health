## signature评分与重大慢病及死亡的cox回归关联分析
# 获取所有结局文件名
outcome_files <- list.files(outcome_folder_name, pattern = "\\.tsv$", full.names = TRUE)

# 检查是否找到了文件
if(length(outcome_files) == 0) {
  stop("未在指定文件夹中找到 .tsv 结局文件，请检查 outcome_folder_name 路径或文件后缀。")
} else {
  print(paste("成功找到", length(outcome_files), "个结局文件。"))
}

# 定义两套协变量
# Model 1: 基础协变量
covariates_model1 <- c("age_test", "MVPA", "season", "sex", "race", "tdi", "smk", "alc", "bmi")

# Model 2: 基础协变量 + 睡眠时长
covariates_model2 <- c(covariates_model1, "SleepDurationInSpt_AD_T5A5_mn")

# 定义一个提取 Cox 结果的辅助函数
# 只需要提取 Q1 (相对于 Q4) 的结果
extract_cox_res <- function(model, outcome_name) {
  sum_cox <- summary(model)
  coefs <- sum_cox$coefficients
  conf_int <- sum_cox$conf.int
  
  # 找到 Q1 的行名 (因为 ref 是 Q4，所以变量名应该是 met_quartileQ1)
  target_row <- grep("met_quartileQ1", rownames(coefs))
  
  if (length(target_row) > 0) {
    res <- data.frame(
      Outcome = outcome_name,
      HR = coefs[target_row, "exp(coef)"],
      CI_low = conf_int[target_row, "lower .95"],
      CI_high = conf_int[target_row, "upper .95"],
      P_value = coefs[target_row, "Pr(>|z|)"]
    )
    return(res)
  } else {
    return(NULL)
  }
}

# 初始化结果列表
results_list_model1 <- list()
results_list_model2 <- list()

# ==============================================================================
# 4. 循环分析结局
# ==============================================================================

for (file_path in outcome_files) {
  
  # 4.1 提取结局名称 (假设文件名即结局名，去掉 .csv)
  raw_name <- tools::file_path_sans_ext(basename(file_path))
  
  # 清理文件名，例如 "anaemia_outcome_260116" -> "anaemia"
  outcome_name <- gsub("_outcome_.*", "", raw_name)
  
  message(paste("正在分析:", outcome_name, "..."))
  
  # 4.2 读取结局数据
  # 假设结局文件里有 eid, target_status, target_time 等列
  outcome_df <- read_tsv(file_path, show_col_types = FALSE) %>% 
    dplyr::select(eid, target_status, target_time) # 根据实际文件列名调整
  
  # 4.3 合并数据
  df_analysis <- df_base_fixed %>%
    inner_join(outcome_df, by = "eid")
  
  # 4.4 过滤 target_time > 0 (图片特别要求) 确保随访时间>0
  df_analysis <- df_analysis %>%
    filter(target_time > 0)
  
  # 4.5 运行 Model 1 (不含睡眠时长)
  # 公式: Surv ~ Q1-Q4 + covariates1
  f1 <- as.formula(paste("Surv(target_time, target_status) ~ met_quartile +", 
                         paste(covariates_model1, collapse = " + ")))
  
  fit1 <- coxph(f1, data = df_analysis)
  res1 <- extract_cox_res(fit1, outcome_name)
  results_list_model1[[outcome_name]] <- res1
  
  # 4.6 运行 Model 2 (含睡眠时长)
  # 公式: Surv ~ Q1-Q4 + covariates2
  f2 <- as.formula(paste("Surv(target_time, target_status) ~ met_quartile +", 
                         paste(covariates_model2, collapse = " + ")))
  
  fit2 <- coxph(f2, data = df_analysis)
  res2 <- extract_cox_res(fit2, outcome_name)
  results_list_model2[[outcome_name]] <- res2
}

# ==============================================================================
# 5. 结果汇总与校正
# ==============================================================================

# --- 处理 Model 1 结果 ---
final_res_model1 <- bind_rows(results_list_model1)
# FDR 校正
final_res_model1$FDR <- p.adjust(final_res_model1$P_value, method = "fdr")
# 整理顺序
final_res_model1 <- final_res_model1 %>%
  dplyr::select(Outcome, HR, CI_low, CI_high, P_value, FDR)

# --- 处理 Model 2 结果 ---
final_res_model2 <- bind_rows(results_list_model2)
# FDR 校正
final_res_model2$FDR <- p.adjust(final_res_model2$P_value, method = "fdr")
# 整理顺序
final_res_model2 <- final_res_model2 %>%
  dplyr::select(Outcome, HR, CI_low, CI_high, P_value, FDR)

# ==============================================================================
# 6. 导出结果到 Excel
# ==============================================================================

wb <- createWorkbook()

# Sheet 1: Model 1
addWorksheet(wb, "Model 1 (Base Covariates)")
writeData(wb, "Model 1 (Base Covariates)", final_res_model1)

# Sheet 2: Model 2
addWorksheet(wb, "Model 2 (+ Sleep Duration)")
writeData(wb, "Model 2 (+ Sleep Duration)", final_res_model2)

# 保存文件
output_filename <- paste0("Results/Cox_Met_Signature_Analysis_Results_", Sys.Date(), ".xlsx")
saveWorkbook(wb, output_filename, overwrite = TRUE)

message(paste("分析完成！结果已保存至:", output_filename))
