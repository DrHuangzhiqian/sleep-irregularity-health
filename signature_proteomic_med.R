############################################准备数据##############################################
## 1. 加载所需 R 包
library(data.table)
library(readr)
library(dplyr)
library(tidyr)
library(survival)

library(mediation)   # 方法1：逻辑回归
library(CMAverse)    # 方法2：Cox 中介分析

## 标准化 SRI 蛋白评分
pro_score$SRI_pro_score_scaled <- as.numeric(
  scale(pro_score$SRI_pro_score)
)

## 仅保留需要的变量
pro_score <- pro_score %>%
  dplyr::select(eid, SRI_pro_score_scaled)


## 3. 导入 GGIR 睡眠数据(暴露)
GGIR_selected <- fread("HZQ-260121/GGIR_selected.csv")  ## 约 91544 人

GGIR_selected <- GGIR_selected %>%
  dplyr::select(
    eid,
    calendar_date,
    SleepDurationInSpt_AD_T5A5_mn,
    SleepRegularityIndex_AD_T5A5_mn,
    age_test,
    MVPA,
    season
  ) %>%
  mutate(season = as.factor(season))

## 4. 导入协变量数据
CovariatesImputed <- read.csv(
  "HZQ-260121/CovariatesImputed.csv",
  stringsAsFactors = FALSE
)  ## 约 502124 人，无 NA

CovariatesImputed <- CovariatesImputed %>%
  mutate(
    across(
      c(sex, race, smk, alc, cl_med, hbp, dm),
      as.factor
    )
  )

## 5. 合并固定底表
df_base_fixed <- pro_score %>%
  inner_join(GGIR_selected, by = "eid") %>%
  inner_join(CovariatesImputed, by = "eid")

print(paste("固定底表人数:", nrow(df_base_fixed)))

## 6. 批量读取结局并分析
# 获取目录下所有符合 "outcome_260116.tsv" 结尾的文件
# full.names = TRUE 保证读取时能找到路径
# 如果你的文件就在 setwd 的目录下，path 可以去掉或写 "."
outcome_files <- list.files(path = outcome_folder_name, pattern = "\\.tsv$", full.names = TRUE)

if(length(outcome_files) == 0) {
  stop("错误：在指定文件夹中未找到 .tsv 文件，请检查路径或文件名！")
} else {
  message(paste("共找到", length(outcome_files), "个结局文件，准备开始循环分析..."))
}

# 定义公式组件
covariates_formula <- "age_test + MVPA + season + sex + race + tdi + smk + alc + bmi" # 注意：请确认 covariates 列名是否都在 df_base_fixed 里
exposure_var <- "SleepRegularityIndex_AD_T5A5_mn"
mediator_var <- "SRI_pro_score_scaled"

# 创建结果存储容器
results_table <- data.frame()
dir.create("results", showWarnings = FALSE) # 创建结果目录

## 7. 开始循环 ---
for (file_path in outcome_files) {
  
  # A. 从文件名提取疾病名称 (例如 "anaemia_outcome_260116.tsv" -> "anaemia")
  file_name <- basename(file_path)
  outcome_name <- gsub("_outcome_260116\\.tsv", "", file_name)
  
  message(paste(">>> 正在处理:", outcome_name, "<<<"))
  
  # B. 读取结局数据
  outcome_data <- read_tsv(file_path, show_col_types = FALSE)
  
  # C. 自动识别 状态列(status) 和 时间列(time)
  # 假设列名包含 "status" 和 "time" (如 diabetes_status, diabetes_time 或 target_status)
  col_names <- names(outcome_data)
  status_col <- col_names[grep("status", col_names, ignore.case = TRUE)]
  time_col   <- col_names[grep("time", col_names, ignore.case = TRUE)]
  
  # 简单检查是否找到列
  if (length(status_col) != 1 || length(time_col) != 1) {
    warning(paste("跳过", outcome_name, ": 无法自动识别唯一的 status 或 time 列名。请检查文件表头。"))
    next
  }
  
  # D. 与底表合并
  # 注意：inner_join 会自动根据 eid 匹配
  df_loop <- df_base_fixed %>%
    inner_join(outcome_data, by = "eid") %>%
    filter(.data[[time_col]] > 0) %>%        # 剔除随访时间 <= 0
    filter(!is.na(.data[[status_col]]))      # 剔除结局状态缺失
  
  # E. 建模
  # (1) Mediator Model
  f_med <- as.formula(paste(mediator_var, "~", exposure_var, "+", covariates_formula))
  med_model <- lm(f_med, data = df_loop)
  
  # (2) Outcome Model (Logistic)
  f_out <- as.formula(paste(status_col, "~", exposure_var, "+", mediator_var, "+", covariates_formula))
  out_model <- glm(f_out, family = binomial(link = "logit"), data = df_loop)
  
  # F. 中介分析 (Bootstrap)
  set.seed(123)
  # 测试时可以将 sims 设为 50 或 100，正式跑建议 500 或 1000
  med_res <- mediate(
    model.m = med_model,
    model.y = out_model,
    treat = exposure_var,
    mediator = mediator_var,
    boot = TRUE,
    sims = 1000 
  )
  
  # G. 提取并汇总结果
  res_summary <- summary(med_res)
  
  temp_res <- data.frame(
    Disease = outcome_name,
    N_Analysis = nrow(df_loop), # 实际进入分析的人数
    
    # ACME (间接效应)
    ACME_Est = res_summary$d0,
    ACME_P   = res_summary$d0.p,
    ACME_CI_Low = res_summary$d0.ci[1],
    ACME_CI_Up  = res_summary$d0.ci[2],
    
    # ADE (直接效应)
    ADE_Est  = res_summary$z0,
    ADE_P    = res_summary$z0.p,
    ADE_CI_Low  = res_summary$z0.ci[1],
    ADE_CI_High = res_summary$z0.ci[2],
    
    # Total Effect
    Total_Est = res_summary$tau.coef,
    Total_P   = res_summary$tau.p,
    Total_Low = res_summary$tau.ci[1],
    Total_High = res_summary$tau.ci[2],
    
    # Proportion Mediated (%)
    Prop_Med_Pct = res_summary$n0 * 100,
    Prop_Med_P   = res_summary$n0.p,
    Prop_Med_Low = res_summary$n0.ci[1],
    Prop_Med_High = res_summary$n0.ci[2]
  )
  
  # 添加到总表
  results_table <- rbind(results_table, temp_res)
  
  # 备份单个 RDS (防止程序崩溃白跑)
  saveRDS(med_res, file = paste0("results/med_res_pro_", outcome_name, ".rds"))
  
  print(paste(outcome_name, "处理完成。"))
}

############################################ 3. 保存最终总表 ##############################################

write_csv(results_table, "results/Final_Mediation_All_16_Diseases_Proteomics.csv")
print("所有分析完成！结果已保存。")

