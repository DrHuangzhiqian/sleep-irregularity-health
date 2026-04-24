## SRI ~ omics(met, pro, infla) 线性回归

library(tidyverse)
library(readxl)
library(data.table)
library(tidyverse)
library(dplyr)

## 1.设置文件路径
setwd("D:/ukb data/Sleep regularity/Sleep regularity/Linear_reg")

## 2. 导入代谢/蛋白/炎症
df_met <- read_tsv("exp_cov_omics/Metabolomic_Processed_260113.tsv",show_col_types = FALSE)
df_pro <- read_tsv("exp_cov_omics/Proteomic_processed_260113.tsv",show_col_types = FALSE)
df_infla <- read_tsv("exp_cov_omics/Blood_inflam_processed_260115.tsv",show_col_types = FALSE)

## 3. 导入GGIR 睡眠数据(暴露)
GGIR_selected <- fread("exp_cov_omics/GGIR_selected.csv")  ## 约 91544 人

## 4. 导入协变量数据
CovariatesImputed <- read.csv("exp_cov_omics/CovariatesImputed.csv",stringsAsFactors = FALSE)  ## 约 502124 人，无 NA

CovariatesImputed <- CovariatesImputed %>%
  mutate(
    across(
      c(sex, race, smk, alc, cl_med, hbp, dm),
      as.factor
    )
  )

# 基本对象
id.var <- "eid"
exposure.var <- "SRI_scaled"
met.vars <- setdiff(colnames(df_met), "eid")
pro.vars <- setdiff(colnames(df_pro), "eid")
infla.vars <- setdiff(colnames(df_infla), "eid")

length(met.vars) #记录数量 met313
length(pro.vars) #pro2920
length(infla.vars) #infla 12+2

# 协变量 
# Model 1
covs_m1 <- c("age_test", "MVPA", "season", "sex", "race", 
             "tdi", "smk", "alc", "bmi", "fastingtime")

# Model 2 (追加睡眠时长)
covs_m2 <- c(covs_m1, "SleepDurationInSpt_AD_T5A5_mn")

## 5. 合并
# 代谢
df_merge_met <- df_met %>%
  inner_join(GGIR_selected, by = "eid") %>%
  inner_join(CovariatesImputed, by = "eid")

print(paste("分析样本量:", nrow(df_merge_met))) #89436
head(df_merge_met)
colnames(df_merge_met)

# 蛋白
df_merge_pro <- df_pro %>%
  inner_join(GGIR_selected, by = "eid") %>%
  inner_join(CovariatesImputed, by = "eid")

print(paste("分析样本量:", nrow(df_merge_pro))) #8938

# 炎症指标
df_merge_infla <- df_infla %>%
  inner_join(GGIR_selected, by = "eid") %>%
  inner_join(CovariatesImputed, by = "eid")

print(paste("分析样本量:", nrow(df_merge_infla))) #87667

## 6. 定义循环函数
run_regression_loop <- function(data, outcome_list, exposure, covariates) {
  
  # 构造公式右半部分
  rhs_formula <- paste(c(exposure, covariates), collapse = " + ")
  
  # 初始化一个空列表来装结果
  results_list <- list()
  
  # 获取总任务数，用于显示进度
  total_vars <- length(outcome_list)
  
  cat(paste0("开始分析: 共 ", total_vars, " 个指标...\n"))
  
  # 开启循环
  for (i in seq_along(outcome_list)) {
    y <- outcome_list[i]
    
    # 打印进度 (每跑 10 个显示一次，防止刷屏)
    if (i %% 10 == 0) cat(paste0("\r进度: ", i, " / ", total_vars, " (", round(i/total_vars*100, 1), "%)"))
    
    # 1. 构造公式
    form <- as.formula(paste(y, "~", rhs_formula))
    
    # 2. 运行模型 (使用 tryCatch 跳过报错的指标)
    fit <- tryCatch({
      lm(form, data = data)
    }, error = function(e) return(NULL))
    
    if (!is.null(fit)) {
      # 3. 提取系数
      coefs <- summary(fit)$coefficients
      
      # 4. 提取暴露变量结果
      if (exposure %in% rownames(coefs)) {
        results_list[[y]] <- data.frame(
          Outcome = y,
          Beta = coefs[exposure, "Estimate"],
          SE = coefs[exposure, "Std. Error"],
          P_value = coefs[exposure, "Pr(>|t|)"],
          N = nobs(fit),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  cat("\n分析完成，正在整理结果...\n")
  
  # 合并结果
  final_df <- do.call(rbind, results_list)
  
  # FDR 校正
  if (!is.null(final_df) && nrow(final_df) > 0) {
    final_df$FDR <- p.adjust(final_df$P_value, method = "BH")
    final_df <- final_df %>% arrange(P_value)
  }
  
  return(final_df)
}

# ==============================================================================
# 7. 执行分析 (Model 1 & Model 2)
# ==============================================================================

## 代谢
# --- Model 1 ---
cat("\n=== 正在运行 Model 1 ===\n")
res_met_model1 <- run_regression_loop(
  data = df_merge_met,
  outcome_list = met.vars,
  exposure = exposure.var,
  covariates = covs_m1
)

# --- Model 2 ---
cat("\n=== 正在运行 Model 2 ===\n")
res_met_model2 <- run_regression_loop(
  data = df_merge_met,
  outcome_list = met.vars,
  exposure = exposure.var,
  covariates = covs_m2
)

## 蛋白
# --- Model 1 ---
cat("\n=== 正在运行 Model 1 ===\n")
res_pro_model1 <- run_regression_loop(
  data = df_merge_pro,
  outcome_list = pro.vars,
  exposure = exposure.var,
  covariates = covs_m1
)

# --- Model 2 ---
cat("\n=== 正在运行 Model 2 ===\n")
res_pro_model2 <- run_regression_loop(
  data = df_merge_pro,
  outcome_list = pro.vars,
  exposure = exposure.var,
  covariates = covs_m2
)

## 炎症指标
# --- Model 1 ---
cat("\n=== 正在运行 Model 1 ===\n")
res_infla_model1 <- run_regression_loop(
  data = df_merge_infla,
  outcome_list = infla.vars,
  exposure = exposure.var,
  covariates = covs_m1
)

# --- Model 2 ---
cat("\n=== 正在运行 Model 2 ===\n")
res_infla_model2 <- run_regression_loop(
  data = df_merge_infla,
  outcome_list = infla.vars,
  exposure = exposure.var,
  covariates = covs_m2
)

# ==============================================================================
# 8. 导出 Excel
# ==============================================================================
library(openxlsx)

cat("\n正在整理最终结果...\n")

# 1. 创建输出列表
# 列表的 "名称" 就是 Excel 下方的 Sheet 名字
# 注意：Excel 的 Sheet 名字不能超过 31 个字符，所以用了简写
output_list <- list(
  "Inflam_Model1"    = res_infla_model1,
  "Inflam_Model2"    = res_infla_model2,
  "Metabo_Model1"    = res_met_model1,
  "Metabo_Model2"    = res_met_model2,
  "Proteo_Model1"    = res_pro_model1,
  "Proteo_Model2"    = res_pro_model2
)

# 2. 定义文件名 (带上今天的日期)
# 格式类似: SRI_Omics_Results_260126.xlsx
file_name <- paste0("SRI_Omics_All_Results_", format(Sys.Date(), "%y%m%d"), ".xlsx")

# 3. 写入 Excel
# write.xlsx 会自动把 list 里的每个 element 存成一个 sheet
write.xlsx(output_list, file = file_name, rowNames = FALSE)

cat(paste0("\n✅ 任务完成！\n"))
cat(paste0("文件已保存为: ", file_name, "\n"))
cat(paste0("保存路径: ", getwd(), "\n"))

# 4. (可选) 简单的结果预览
cat("\n--- 结果概览 (显著 P < 0.05 / FDR < 0.05 的数量) ---\n")
count_sig <- function(df, name) {
  if(is.null(df)) return(paste0(name, ": 结果为空"))
  n_p <- sum(df$P_value < 0.05, na.rm = TRUE)
  n_fdr <- sum(df$FDR < 0.05, na.rm = TRUE)
  return(paste0(name, " -> P<0.05: ", n_p, " 个 | FDR<0.05: ", n_fdr, " 个"))
}

print(count_sig(res_infla_model1, "炎症 Model1"))
print(count_sig(res_met_model1,   "代谢 Model1"))
print(count_sig(res_pro_model1,   "蛋白 Model1"))

print(count_sig(res_infla_model2, "炎症 Model2"))
print(count_sig(res_met_model2,   "代谢 Model2"))
print(count_sig(res_pro_model2,   "蛋白 Model2"))