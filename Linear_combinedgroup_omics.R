## =============================================================================
## 睡眠规律性与时长联合分组 (combined_group) ~ 组学指标线性回归
## 暴露: combined_group (Reference = normal_high)
## 协变量: age_test + MVPA + season + sex + race + tdi + smk + alc + bmi + fastingtime
## 结局: 炎症、代谢组、蛋白组
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(doParallel)
  library(foreach)
  library(openxlsx)
})

# ==============================================================================
# 1. 路径设置与准备底表
# ==============================================================================
# 设置基础工作目录
setwd("D:/ukb data/Sleep regularity/Sleep regularity")

# 创建专门用于存放联合分组线性回归结果的输出文件夹
output_dir <- "Linear_reg/category_combined_linea"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("正在读取组学结局数据...\n")
df_met   <- read_tsv("Metabolomic_Processed_260113.tsv", show_col_types = FALSE)
df_pro   <- read_tsv("Proteomic_processed_260113.tsv", show_col_types = FALSE)
df_infla <- read_tsv("Blood_inflam_processed_260115.tsv", show_col_types = FALSE)

cat("正在读取暴露与协变量数据...\n")
# 1. 基础协变量
GGIR_selected <- fread("GGIR_selected.csv") %>%
  as_tibble() %>%
  dplyr::select(eid, age_test, MVPA, season) %>%
  mutate(season = as.factor(season))

# 2. 读取协变量
CovariatesImputed <- read.csv("CovariatesImputed.csv", stringsAsFactors = FALSE) %>%
  mutate(across(c(sex, race, smk, alc, cl_med, hbp, dm), as.factor))

# 3. 读取最新的联合分组数据
combined_data <- read.csv("GGIR_combined_group_0301.csv", stringsAsFactors = FALSE)

# 4. 合并并处理暴露分组
df_base <- combined_data %>%
  inner_join(GGIR_selected, by = "eid") %>%
  inner_join(CovariatesImputed, by = "eid") %>%
  mutate(
    # 直接转换为因子
    combined_group = as.factor(combined_group)
  )

# 设定 normal_high 为参考组
df_base$combined_group <- relevel(df_base$combined_group, ref = "normal_high")

cat("联合分组分布情况如下：\n")
print(table(df_base$combined_group))

# ==============================================================================
# 2. 定义分析变量与合并组学数据
# ==============================================================================
met_vars   <- setdiff(colnames(df_met), "eid")
pro_vars   <- setdiff(colnames(df_pro), "eid")
infla_vars <- setdiff(colnames(df_infla), "eid")

cat(sprintf("总计: 代谢物 %d 个, 蛋白质 %d 个, 炎症指标 %d 个\n", 
            length(met_vars), length(pro_vars), length(infla_vars)))

# 定义唯一的协变量公式 
covs_formula <- "age_test + MVPA + season + sex + race + tdi + smk + alc + bmi + fastingtime"

# 合并生成 3 个分析主表
df_merge_met   <- df_met %>% inner_join(df_base, by = "eid") %>% drop_na()
df_merge_pro   <- df_pro %>% inner_join(df_base, by = "eid") %>% drop_na()
df_merge_infla <- df_infla %>% inner_join(df_base, by = "eid") %>% drop_na()

cat(sprintf("有效分析样本量 - 代谢: %d | 蛋白: %d | 炎症: %d\n", 
            nrow(df_merge_met), nrow(df_merge_pro), nrow(df_merge_infla)))

# ==============================================================================
# 3. 设置并行计算核心函数
# ==============================================================================
# 您这边设置了8个核
num_cores <- 8
cat(sprintf("\n启动并行计算环境，使用核心数: %d\n", num_cores))

run_parallel_lm <- function(data, outcome_list, dataset_name, covariates) {
  
  cat(paste0("开始并行拟合 ", dataset_name, " 数据 (共 ", length(outcome_list), " 个指标)...\n"))
  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # 【修复点】：在这里通过 .export 将 covariates 传递进并行后台
  results <- foreach(
    y = outcome_list, 
    .combine = rbind,
    .packages = c("tidyverse", "stats"),
    .export = "covariates" 
  ) %dopar% {
    
    # 构造公式，使用传进来的 covariates 变量
    form <- as.formula(paste0(y, " ~ combined_group + ", covariates))
    
    tryCatch({
      fit <- lm(form, data = data)
      coefs <- summary(fit)$coefficients
      
      # 找出所有表示 combined_group 的行
      target_rows <- grep("^combined_group", rownames(coefs), value = TRUE)
      
      # 遍历 5 个对比组提取结果
      res_list <- lapply(target_rows, function(tr) {
        comp_group <- sub("^combined_group", "", tr)
        
        data.frame(
          Outcome = y,
          N_Analysis = nobs(fit),
          Comparison = paste0(comp_group, " vs normal_high"),
          Beta = coefs[tr, "Estimate"],
          SE = coefs[tr, "Std. Error"],
          P_Value = coefs[tr, "Pr(>|t|)"],
          stringsAsFactors = FALSE
        )
      })
      
      do.call(rbind, res_list)
      
    }, error = function(e) {
      return(NULL)
    })
  }
  
  stopCluster(cl)
  
  # 后处理：按 Comparison 分组进行 FDR 校正
  if (!is.null(results) && nrow(results) > 0) {
    results <- results %>%
      as_tibble() %>%
      group_by(Comparison) %>%
      mutate(P_FDR = p.adjust(P_Value, method = "BH")) %>%
      ungroup() %>%
      arrange(Comparison, P_FDR) %>%
      dplyr::select(Outcome, N_Analysis, Comparison, Beta, SE, P_Value, P_FDR)
  }
  
  cat(paste0(">>> ", dataset_name, " 并行计算与 FDR 校正完成！\n"))
  return(results)
}

# ==============================================================================
# 4. 执行分析并导出独立的 Excel 结果
# ==============================================================================

# 注意这里在调用函数时，把 covs_formula 也传进去了
# --- A. 炎症指标 (Inflammation) ---
res_infla <- run_parallel_lm(df_merge_infla, infla_vars, "Inflammation", covs_formula)
path_infla <- file.path(output_dir, paste0("Linear_Inflam_Combined_Category_", Sys.Date(), ".xlsx"))
write.xlsx(list("Inflammation" = res_infla), file = path_infla, rowNames = FALSE)

# --- B. 代谢组学 (Metabolomics) ---
res_met <- run_parallel_lm(df_merge_met, met_vars, "Metabolomics", covs_formula)
path_met <- file.path(output_dir, paste0("Linear_Metabo_Combined_Category_", Sys.Date(), ".xlsx"))
write.xlsx(list("Metabolomics" = res_met), file = path_met, rowNames = FALSE)

# --- C. 蛋白质组学 (Proteomics) ---
res_pro <- run_parallel_lm(df_merge_pro, pro_vars, "Proteomics", covs_formula)
path_pro <- file.path(output_dir, paste0("Linear_Proteo_Combined_Category_", Sys.Date(), ".xlsx"))
write.xlsx(list("Proteomics" = res_pro), file = path_pro, rowNames = FALSE)

# ==============================================================================
# 5. 完成提示与简单统计
# ==============================================================================
cat("\n======================================================\n")
cat("✅ 所有组学分析已全部完成！\n")
cat("结果文件已分别保存至:\n")
cat(" - 炎症: ", path_infla, "\n")
cat(" - 代谢: ", path_met, "\n")
cat(" - 蛋白: ", path_pro, "\n")
cat("======================================================\n")

cat("\n【显著结果概览 (FDR < 0.05)】\n")
if(!is.null(res_infla)) {
  cat("\n--- 炎症指标 ---\n")
  print(res_infla %>% filter(P_FDR < 0.05) %>% count(Comparison))
}
if(!is.null(res_met)) {
  cat("\n--- 代谢物 ---\n")
  print(res_met %>% filter(P_FDR < 0.05) %>% count(Comparison))
}
if(!is.null(res_pro)) {
  cat("\n--- 蛋白质 ---\n")
  print(res_pro %>% filter(P_FDR < 0.05) %>% count(Comparison))
}