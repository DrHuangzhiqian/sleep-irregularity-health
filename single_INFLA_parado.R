## 炎症指标，单个指标的中介分析

# ==============================================================================
# 0. 加载必要的包
# ==============================================================================
library(tidyverse)
library(data.table)
library(survival)
library(openxlsx) 
library(mediation)  # 核心中介分析包
library(parallel)   # 并行计算
library(foreach)
library(doParallel)

# ==============================================================================
# 2. 筛选中介变量 & 导入数据
# ==============================================================================

# --- 2.1 读取线性回归结果并筛选 (Step 1) ---
lm_result_file <- "D:/ukb data/Sleep regularity/Sleep regularity/Linear_reg/SRI_Omics_All_Results_260126.xlsx" 

cat("正在读取线性回归结果：", lm_result_file, "\n")

# 读取 Inflam_Model2 (校正了睡眠时长)
df_lm_res <- read.xlsx(lm_result_file, sheet = "Inflam_Model2")

# 筛选 FDR < 0.05 的显著炎症指标
significant_vars <- df_lm_res %>%
  filter(FDR < 0.05) %>%
  pull(Outcome) # 提取指标名称

cat("Inflam_Model2 中 FDR < 0.05 的潜在中介变量数量:", length(significant_vars), "\n")

if(length(significant_vars) == 0) stop("没有显著的炎症指标，无法进行中介分析！")

# --- 2.2 导入炎症原始数据 (Step 2) ---
inflam_file_path <- "Blood_inflam_processed_260115.tsv" 

if(!file.exists(inflam_file_path)) {
  stop(paste("找不到炎症数据文件:", inflam_file_path, "请检查路径是否正确。"))
}

# 读取炎症数据
df_inflam <- read_tsv(inflam_file_path, show_col_types = FALSE)

# 取交集：确保筛选出来的名字在原始数据里也能找到
target_mediators <- intersect(significant_vars, colnames(df_inflam))
cat("最终纳入分析的中介变量数量 (交集):", length(target_mediators), "\n")

# 仅保留 eid 和 目标中介变量
df_inflam_subset <- df_inflam %>% 
  dplyr::select(all_of(c("eid", target_mediators)))

# --- 2.3 导入 GGIR 睡眠数据 (含 SRI 和 睡眠时长) ---
GGIR_selected <- fread("HZQ-260121/GGIR_selected.csv") 

GGIR_selected <- GGIR_selected %>%
  dplyr::select(
    eid,
    SleepDurationInSpt_AD_T5A5_mn,   # 协变量：睡眠时长
    SleepRegularityIndex_AD_T5A5_mn, # 暴露：SRI
    age_test,
    MVPA,
    season
  ) %>%
  mutate(season = as.factor(season)) %>%
  # 标准化 SRI (暴露变量)
  mutate(SRI_scaled = as.numeric(scale(SleepRegularityIndex_AD_T5A5_mn)))

# --- 2.4 导入协变量数据 ---
CovariatesImputed <- read.csv("HZQ-260121/CovariatesImputed.csv", stringsAsFactors = FALSE)

CovariatesImputed <- CovariatesImputed %>%
  mutate(across(c(sex, race, smk, alc, cl_med, hbp, dm), as.factor))

# --- 2.5 合并所有底表 ---
df_base <- GGIR_selected %>%
  inner_join(CovariatesImputed, by = "eid") %>%
  inner_join(df_inflam_subset, by = "eid") %>%
  drop_na() # 去除缺失值 (中介分析对NA敏感)

cat("最终纳入分析的样本量 (完整数据):", nrow(df_base), "\n") #87667

# ==============================================================================
# 3. 配置并行计算与模型公式
# ==============================================================================

# 定义暴露变量
exposure_var <- "SRI_scaled"

# 定义协变量公式 (包含 fastingtime 和 睡眠时长)
covariates_string <- "age_test + MVPA + season + sex + race + tdi + smk + alc + bmi + fastingtime + SleepDurationInSpt_AD_T5A5_mn"

# 设置并行核心数 (保留2个核给系统，其余全部占用)
num_cores <- detectCores() - 12
# 防止核心数过少
if(num_cores < 1) num_cores <- 1

cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("已启动并行计算，使用核心数:", num_cores, "\n")

# ==============================================================================
# 4. 循环分析 (结局 -> 中介)
# ==============================================================================

# 获取所有结局文件
outcome_files <- list.files(outcome_folder_name, pattern = "\\.tsv$", full.names = TRUE)

for (outcome_file in outcome_files) {
  
  # --- 4.1 准备当前结局的数据 ---
  raw_name <- tools::file_path_sans_ext(basename(outcome_file))
  outcome_name <- gsub("_outcome_.*", "", raw_name) # 清洗文件名
  
  cat(paste0("\n>>> 正在分析结局: ", outcome_name, " ...\n"))
  
  outcome_df <- read_tsv(outcome_file, show_col_types = FALSE) %>% 
    dplyr::select(eid, target_status, target_time)
  
  # 合并结局数据
  df_analysis <- df_base %>%
    inner_join(outcome_df, by = "eid") %>%
    filter(target_time > 0) # 过滤无效生存时间
  
  # 【修正1：必须转为 data.frame，否则 mediate 函数会报错】
  df_analysis <- as.data.frame(df_analysis)
  
  # --- 4.2 并行遍历所有中介变量 ---
  # 使用 foreach 并行跑每一个 mediator
  
  results_outcome <- foreach(mediator = target_mediators, 
                             .packages = c("mediation", "stats"),
                             .combine = rbind) %dopar% {
                               
    tryCatch({
      # A. 构建 Mediator 模型 (Model M): Linear Regression
      # Mediator ~ Exposure + Covariates
      f_m <- as.formula(paste(mediator, "~", exposure_var, "+", covariates_string))
      model_m <- lm(f_m, data = df_analysis)
      
      # 【关键修正 1】强制将公式对象植入 call 中，防止 "object 'f_m' not found"
      model_m$call$formula <- f_m
      
      # B. 构建 Outcome 模型 (Model Y): Logistic Regression (Status 是 0/1)
      # 逻辑回归：target_status ~ Exposure + Mediator + Covs
      f_y <- as.formula(paste("target_status ~", exposure_var, "+", mediator, "+", covariates_string))
      model_y <- glm(f_y, family = binomial(link = "logit"), data = df_analysis)
      
      # 【关键修正 2】同上，强制植入公式
      model_y$call$formula <- f_y
      
      # C. 执行中介分析 (Mediate Function)
      set.seed(123) 
      # sims=100 用于快速筛选。
      med_out <- mediate(model_m, model_y, 
                         treat = exposure_var, 
                         mediator = mediator, 
                         boot = TRUE,
                         sims = 1000) 
      
      # D. 提取结果
      res_summary <- summary(med_out)
      
      temp_res <- data.frame(
        Disease = outcome_name,
        Mediator = mediator,
        N_Analysis = nrow(df_analysis),
        
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
      
      return(temp_res)
      
    }, error = function(e) {
      return(NULL) # 出错跳过
    })
  } 
  
  # --- 4.3 保存当前结局的结果 ---
  if (!is.null(results_outcome) && nrow(results_outcome) > 0) {
    save_path <- paste0("Results_Mediation_Inflam/Mediation_Infla_", outcome_name, ".csv")
    write.csv(results_outcome, save_path, row.names = FALSE)
    cat(paste("   已保存结果:", save_path, "\n"))
  }
}


# ==============================================================================
# 5. 结束并行与合并最终结果
# ==============================================================================
stopCluster(cl)
cat("\n并行计算结束。\n")

# 合并所有结局的 CSV
all_files <- list.files("Results_Mediation_Inflam", pattern = "\\.csv$", full.names = TRUE)

if(length(all_files) > 0) {
  final_mediation_res <- bind_rows(lapply(all_files, read.csv))
  
  # 导出最终 Excel
  wb <- createWorkbook()
  addWorksheet(wb, "Inflammation Mediation")
  writeData(wb, "Inflammation Mediation", final_mediation_res)
  saveWorkbook(wb, paste0("Results_Mediation_Inflam/Final_Inflam_Mediation_Results_", Sys.Date(), ".xlsx"), overwrite = TRUE)
  
  cat("全部完成！最终汇总文件已生成。\n")
} else {
  cat("警告：没有生成任何结果文件，请检查数据或模型是否出错。\n")
}