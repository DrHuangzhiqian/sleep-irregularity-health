## =============================================================================
## 睡眠规律性与时长联合暴露 Cox 回归 (分类变量)
## 暴露: combined_group (Reference = normal_high)
## 协变量: age_test + MVPA + season + sex + race + tdi + smk + alc + bmi
## 目标: 提取其余 5 个分类的原始结果，并进行全局 FDR 和 BFI (P*16) 校正
## =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(survival)
  library(doParallel)
  library(foreach)
})

# ==============================================================================
# 1. 路径设置与数据读取
# ==============================================================================
setwd("D:/ukb data/Sleep regularity/Sleep regularity")

cat("读取基础协变量与联合分组数据...\n")

# 1. 基础协变量 (注意：排除了睡眠时长，因为它已经包含在分组里了)
GGIR_selected <- fread("GGIR_selected.csv") %>%
  as_tibble() %>%
  dplyr::select(eid, age_test, MVPA, season) %>%
  mutate(season = as.factor(season))

CovariatesImputed <- read.csv("CovariatesImputed.csv", stringsAsFactors = FALSE) %>%
  mutate(across(c(sex, race, smk, alc, cl_med, dm), as.factor))

# 2. 读取新的联合分组数据
combined_data <- read.csv("GGIR_combined_group_0301.csv", stringsAsFactors = FALSE)

# 3. 合并数据并设置参考组
df_base <- combined_data %>%
  inner_join(GGIR_selected, by = "eid") %>%
  inner_join(CovariatesImputed, by = "eid") %>%
  mutate(
    # 将 combined_group 转换为因子
    combined_group = as.factor(combined_group)
  )

# 【关键修改】：设定 normal_high 为参考组
df_base$combined_group <- relevel(df_base$combined_group, ref = "normal_high")

cat("数据合并完成。联合分组分布情况:\n")
print(table(df_base$combined_group))

# ==============================================================================
# 2. 定义分析变量
# ==============================================================================
# 更新协变量公式 (不含睡眠时长)
common_covs <- "age_test + MVPA + season + sex + race + tdi + smk + alc + bmi"

# 找到所有结局文件
outcome_folder <- "16_outcomes"
outcome_files <- list.files(outcome_folder, pattern = "_outcome_.*\\.tsv$", full.names = TRUE)

if(length(outcome_files) == 0) stop("未找到结局文件，请检查文件夹路径！")

# ==============================================================================
# 3. 设置并行计算环境
# ==============================================================================
num_cores <- 8
cat(sprintf("启动并行计算，使用核心数: %d / %d\n", num_cores, detectCores()))

cl <- makeCluster(num_cores)
registerDoParallel(cl)

# ==============================================================================
# 4. 执行并行循环
# ==============================================================================
cat("开始并行拟合 16 个结局的 Cox 模型...\n")

results_df <- foreach(
  out_file = outcome_files, 
  .combine = rbind, 
  .packages = c("tidyverse", "survival", "readr")
) %dopar% {
  
  disease_name <- str_split(basename(out_file), "_outcome_")[[1]][1]
  
  outcome_df <- readr::read_tsv(out_file, show_col_types = FALSE) %>% 
    dplyr::select(eid, target_status, target_time)
  
  # 合并数据并处理缺失值/随访时间 (维持主分析的 target_time > 0)
  df_analysis <- df_base %>%
    inner_join(outcome_df, by = "eid") %>%
    filter(target_time > 0) %>%
    drop_na()
  
  # 如果发病人数太少，跳过
  if(nrow(df_analysis) < 50 || sum(df_analysis$target_status) < 10) {
    return(NULL)
  }
  
  # 拟合 Cox 模型
  f_formula <- as.formula(paste0("Surv(target_time, target_status) ~ combined_group + ", common_covs))
  
  tryCatch({
    fit <- coxph(f_formula, data = df_analysis)
    
    # 提取模型汇总信息
    sum_fit <- summary(fit)
    coefs <- sum_fit$coefficients
    confs <- sum_fit$conf.int
    
    # 动态匹配由于 combined_group 产生的 5 个虚拟变量行
    target_rows <- grep("^combined_group", rownames(coefs), value = TRUE)
    
    # 循环提取这 5 个对比组的结果
    res_list <- lapply(target_rows, function(tr) {
      # 提取实际的组名，去掉 "combined_group" 前缀
      comp_group <- sub("^combined_group", "", tr)
      
      data.frame(
        Outcome = disease_name,
        N_Analysis = nrow(df_analysis),
        Cases = sum(df_analysis$target_status),
        # 【关键修改】：这里的对比描述改为了 vs normal_high
        Comparison = paste0(comp_group, " vs normal_high"), 
        HR = coefs[tr, "exp(coef)"],
        LCI = confs[tr, "lower .95"],
        UCI = confs[tr, "upper .95"],
        P_Value = coefs[tr, "Pr(>|z|)"]
      )
    })
    
    # 将列表绑定为数据框返回
    do.call(rbind, res_list)
    
  }, error = function(e) {
    return(NULL) # 如果模型不收敛或报错，返回 NULL 不打断循环
  })
}

stopCluster(cl)
cat("并行计算完成！\n")

# ==============================================================================
# 5. 后处理与 多重检验校正 (FDR & Bonferroni)
# ==============================================================================
# 整理最终数据并进行校正
final_results <- results_df %>%
  as_tibble() %>%
  mutate(
    # 1. 全局 FDR 校正
    P_FDR = p.adjust(P_Value, method = "fdr"),
    
    # 2. 按照您的要求计算 BFI: 原始 P * 16 (最大不超过 1.0)
    P_BFI = pmin(P_Value * 16, 1.0)
  ) %>%
  # 排序：先按疾病，再按比较组
  arrange(Outcome, Comparison) %>%
  # 保留高精度纯数字格式
  dplyr::select(
    Outcome, N_Analysis, Cases, Comparison, 
    HR, LCI, UCI, P_Value, P_FDR, P_BFI
  )

# 打印预览
print(head(final_results, 10))

# ==============================================================================
# 6. 保存为 CSV (可以直接用 Excel 打开)
# ==============================================================================
output_dir <- "Cox_Main_sensi_Analysis/Combined_SRI_Duration"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

save_path <- file.path(output_dir, paste0("Combined_Duration_Cate_normal_Results_", Sys.Date(), ".csv"))
write.csv(final_results, save_path, row.names = FALSE)

cat("\n所有联合分组分析结果 (含 FDR 和 BFI 校正) 已保存至:\n", save_path, "\n")