library(data.table)
library(dplyr)
library(readr)
library(broom)
library(survival)
library(openxlsx)
library(forcats)
library(tools)
library(parallel)
library(doParallel)
library(foreach)

# ==============================================================
# 协变量定义
# ==============================================================
covariates1 <- c("age_test", "MVPA", "season", "sex", "race", "tdi", "smk", "alc", "bmi")
covariates2 <- c(covariates1, "SleepDurationInSpt_AD_T5A5_mn")

# ==============================================================
# 4. 并行设置
# ==============================================================
num_cores <- 8
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("已启动并行计算，使用核心数:", num_cores, "\n")

# ==============================================================
# 5. 通用函数
# ==============================================================

# ---- 5.1 添加四分位变量 ----
add_sri_quartile <- function(df) {
  df %>%
    mutate(
      sri_q = ntile(SRI, 4),
      sri_q = factor(
        paste0("Q", sri_q),
        levels = c("Q1", "Q2", "Q3", "Q4")
      ),
      sri_q = fct_relevel(sri_q, "Q4")
    )
}

# ---- 5.2 提取 Logistic 结果（OR） ----
tidy_fast_sri_glm <- function(fit) {
  broom::tidy(fit) %>%
    filter(grepl("^sri_qQ[1-3]$", term)) %>%
    mutate(
      contrast = case_when(
        term == "sri_qQ1" ~ "Q1 vs Q4",
        term == "sri_qQ2" ~ "Q2 vs Q4",
        term == "sri_qQ3" ~ "Q3 vs Q4"
      ),
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error),
      P_value = p.value
    ) %>%
    select(contrast, OR, CI_low, CI_high, P_value)
}

# ---- 5.3 提取 Cox 结果（HR） ----
tidy_fast_sri_cox <- function(fit) {
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(grepl("^sri_qQ[1-3]$", term)) %>%
    mutate(
      contrast = case_when(
        term == "sri_qQ1" ~ "Q1 vs Q4",
        term == "sri_qQ2" ~ "Q2 vs Q4",
        term == "sri_qQ3" ~ "Q3 vs Q4"
      ),
      HR = exp(estimate),
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value
    ) %>%
    select(contrast, HR, CI_low, CI_high, P_value)
}

# ---- 5.4 提取连续SRI的P值 ----
tidy_sri_trend_glm <- function(fit, term_name = "SRI") {
  broom::tidy(fit) %>%
    filter(term == term_name) %>%
    transmute(
      contrast = "P trend",
      P_value = p.value
    )
}

tidy_sri_trend_cox <- function(fit, term_name = "SRI") {
  broom::tidy(fit, conf.int = TRUE) %>%
    filter(term == term_name) %>%
    transmute(
      contrast = "P trend",
      P_value = p.value
    )
}

# ---- 5.5 拟合 Logistic ----
fit_glm_model <- function(data, outcome_var, exposure_var, covariates) {
  rhs <- paste(c(exposure_var, covariates), collapse = " + ")
  fml <- as.formula(paste0(outcome_var, " ~ ", rhs))
  glm(fml, family = binomial, data = data)
}

# ---- 5.6 拟合 Cox ----
fit_cox_model <- function(data, time_var, status_var, exposure_var, covariates) {
  rhs <- paste(c(exposure_var, covariates), collapse = " + ")
  fml <- as.formula(paste0("Surv(", time_var, ", ", status_var, ") ~ ", rhs))
  coxph(fml, data = data)
}

# ---- 5.7 统计 N 和 Cases ----
get_prev_summary <- function(df) {
  data.frame(
    N = nrow(df),
    Cases = sum(df$target_prev == 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

get_cox_summary <- function(df) {
  data.frame(
    N = nrow(df),
    Cases = sum(df$target_status == 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# ==============================================================
# 6. 并行遍历所有结局
# ==============================================================

result_list <- foreach(
  outcome_file = outcome_files,
  .packages = c("data.table", "dplyr", "readr", "broom", "survival", "openxlsx", "forcats", "tools")
) %dopar% {
  
  tryCatch({
    
    # ----------------------------------------------------------
    # 6.1 当前结局名称
    # ----------------------------------------------------------
    raw_name <- tools::file_path_sans_ext(basename(outcome_file))
    outcome_name <- gsub("_outcome_.*", "", raw_name)
    
    # 如 death 文件名可能不止一种，可在这里扩展
    is_death <- outcome_name %in% c("death", "all_cause_death")
    
    # ----------------------------------------------------------
    # 6.2 读取结局数据
    # ----------------------------------------------------------
    outcome_df <- read_tsv(outcome_file, show_col_types = FALSE)
    
    required_cols <- if (is_death) {
      c("eid", "target_status", "target_time")
    } else {
      c("eid", "target_baseline", "target_status", "target_time")
    }
    
    missing_cols <- setdiff(required_cols, colnames(outcome_df))
    if (length(missing_cols) > 0) {
      stop(paste0("缺少必要列: ", paste(missing_cols, collapse = ", ")))
    }
    
    outcome_df <- outcome_df %>%
      dplyr::select(all_of(required_cols))
    
    # ----------------------------------------------------------
    # 6.3 合并总数据
    # ----------------------------------------------------------
    df_merge <- GGIR_selected %>%
      inner_join(CovariatesImputed, by = "eid") %>%
      left_join(outcome_df, by = "eid")
    
    wb <- createWorkbook()
    sheet_written <- c()
    
    # 用于收集当前结局所有结果，便于最终总汇总
    all_results_current <- list()
    
    # ==========================================================
    # 6.4 death：只做 incident Cox
    # ==========================================================
    if (is_death) {
      
      df_cox <- df_merge %>%
        filter(!is.na(SRI), !is.na(target_time), !is.na(target_status), target_time >= 0) %>%
        add_sri_quartile()
      
      if (nrow(df_cox) == 0) {
        stop("death结局可用于Cox分析的数据为空")
      }
      
      cox_summary <- get_cox_summary(df_cox)
      
      # ---- incident_cov1 ----
      fit_cox1 <- fit_cox_model(
        data = df_cox,
        time_var = "target_time",
        status_var = "target_status",
        exposure_var = "sri_q",
        covariates = covariates1
      )
      
      res_cox1 <- tidy_fast_sri_cox(fit_cox1) %>%
        mutate(
          Outcome = outcome_name,
          Analysis = "incident",
          Model = "incident_cov1",
          N = cox_summary$N,
          Cases = cox_summary$Cases,
          .before = 1
        )
      
      fit_cox1_trend <- fit_cox_model(
        data = df_cox,
        time_var = "target_time",
        status_var = "target_status",
        exposure_var = "SRI",
        covariates = covariates1
      )
      
      res_cox1_trend <- tidy_sri_trend_cox(fit_cox1_trend) %>%
        mutate(
          Outcome = outcome_name,
          Analysis = "incident",
          Model = "incident_cov1",
          N = cox_summary$N,
          Cases = cox_summary$Cases,
          .before = 1
        )
      
      # ---- incident_cov2 ----
      fit_cox2 <- fit_cox_model(
        data = df_cox,
        time_var = "target_time",
        status_var = "target_status",
        exposure_var = "sri_q",
        covariates = covariates2
      )
      
      res_cox2 <- tidy_fast_sri_cox(fit_cox2) %>%
        mutate(
          Outcome = outcome_name,
          Analysis = "incident",
          Model = "incident_cov2",
          N = cox_summary$N,
          Cases = cox_summary$Cases,
          .before = 1
        )
      
      fit_cox2_trend <- fit_cox_model(
        data = df_cox,
        time_var = "target_time",
        status_var = "target_status",
        exposure_var = "SRI",
        covariates = covariates2
      )
      
      res_cox2_trend <- tidy_sri_trend_cox(fit_cox2_trend) %>%
        mutate(
          Outcome = outcome_name,
          Analysis = "incident",
          Model = "incident_cov2",
          N = cox_summary$N,
          Cases = cox_summary$Cases,
          .before = 1
        )
      
      # 保存sheet
      addWorksheet(wb, "incident_cov1")
      writeData(wb, "incident_cov1", res_cox1)
      sheet_written <- c(sheet_written, "incident_cov1")
      
      addWorksheet(wb, "incident_cov2")
      writeData(wb, "incident_cov2", res_cox2)
      sheet_written <- c(sheet_written, "incident_cov2")
      
      p_trend_inc <- bind_rows(res_cox1_trend, res_cox2_trend)
      addWorksheet(wb, "P_trend_inc")
      writeData(wb, "P_trend_inc", p_trend_inc)
      sheet_written <- c(sheet_written, "P_trend_inc")
      
      all_results_current[["incident_cov1"]] <- res_cox1
      all_results_current[["incident_cov2"]] <- res_cox2
      all_results_current[["P_trend_inc"]] <- p_trend_inc
    }
    
    # ==========================================================
    # 6.5 其他疾病：prevalent + incident
    # ==========================================================
    if (!is_death) {
      
      # ---- Prevalent ----
      df_prev <- df_merge %>%
        mutate(target_prev = target_baseline) %>%
        filter(!is.na(target_prev), !is.na(SRI)) %>%
        add_sri_quartile()
      
      # ---- Incident ----
      df_cox <- df_merge %>%
        filter(target_baseline == 0) %>%
        filter(!is.na(SRI), !is.na(target_time), !is.na(target_status), target_time >= 0) %>%
        add_sri_quartile()
      
      # -------------------- prevalent --------------------
      if (nrow(df_prev) > 0) {
        
        prev_summary <- get_prev_summary(df_prev)
        
        fit_prev1 <- fit_glm_model(
          data = df_prev,
          outcome_var = "target_prev",
          exposure_var = "sri_q",
          covariates = covariates1
        )
        
        res_prev1 <- tidy_fast_sri_glm(fit_prev1) %>%
          mutate(
            Outcome = outcome_name,
            Analysis = "prevalent",
            Model = "prevalent_cov1",
            N = prev_summary$N,
            Cases = prev_summary$Cases,
            .before = 1
          )
        
        fit_prev1_trend <- fit_glm_model(
          data = df_prev,
          outcome_var = "target_prev",
          exposure_var = "SRI",
          covariates = covariates1
        )
        
        res_prev1_trend <- tidy_sri_trend_glm(fit_prev1_trend) %>%
          mutate(
            Outcome = outcome_name,
            Analysis = "prevalent",
            Model = "prevalent_cov1",
            N = prev_summary$N,
            Cases = prev_summary$Cases,
            .before = 1
          )
        
        fit_prev2 <- fit_glm_model(
          data = df_prev,
          outcome_var = "target_prev",
          exposure_var = "sri_q",
          covariates = covariates2
        )
        
        res_prev2 <- tidy_fast_sri_glm(fit_prev2) %>%
          mutate(
            Outcome = outcome_name,
            Analysis = "prevalent",
            Model = "prevalent_cov2",
            N = prev_summary$N,
            Cases = prev_summary$Cases,
            .before = 1
          )
        
        fit_prev2_trend <- fit_glm_model(
          data = df_prev,
          outcome_var = "target_prev",
          exposure_var = "SRI",
          covariates = covariates2
        )
        
        res_prev2_trend <- tidy_sri_trend_glm(fit_prev2_trend) %>%
          mutate(
            Outcome = outcome_name,
            Analysis = "prevalent",
            Model = "prevalent_cov2",
            N = prev_summary$N,
            Cases = prev_summary$Cases,
            .before = 1
          )
        
        addWorksheet(wb, "prevalent_cov1")
        writeData(wb, "prevalent_cov1", res_prev1)
        sheet_written <- c(sheet_written, "prevalent_cov1")
        
        addWorksheet(wb, "prevalent_cov2")
        writeData(wb, "prevalent_cov2", res_prev2)
        sheet_written <- c(sheet_written, "prevalent_cov2")
        
        p_trend_prev <- bind_rows(res_prev1_trend, res_prev2_trend)
        addWorksheet(wb, "P_trend_prev")
        writeData(wb, "P_trend_prev", p_trend_prev)
        sheet_written <- c(sheet_written, "P_trend_prev")
        
        all_results_current[["prevalent_cov1"]] <- res_prev1
        all_results_current[["prevalent_cov2"]] <- res_prev2
        all_results_current[["P_trend_prev"]] <- p_trend_prev
      }
      
      # -------------------- incident --------------------
      if (nrow(df_cox) > 0) {
        
        cox_summary <- get_cox_summary(df_cox)
        
        fit_cox1 <- fit_cox_model(
          data = df_cox,
          time_var = "target_time",
          status_var = "target_status",
          exposure_var = "sri_q",
          covariates = covariates1
        )
        
        res_cox1 <- tidy_fast_sri_cox(fit_cox1) %>%
          mutate(
            Outcome = outcome_name,
            Analysis = "incident",
            Model = "incident_cov1",
            N = cox_summary$N,
            Cases = cox_summary$Cases,
            .before = 1
          )
        
        fit_cox1_trend <- fit_cox_model(
          data = df_cox,
          time_var = "target_time",
          status_var = "target_status",
          exposure_var = "SRI",
          covariates = covariates1
        )
        
        res_cox1_trend <- tidy_sri_trend_cox(fit_cox1_trend) %>%
          mutate(
            Outcome = outcome_name,
            Analysis = "incident",
            Model = "incident_cov1",
            N = cox_summary$N,
            Cases = cox_summary$Cases,
            .before = 1
          )
        
        fit_cox2 <- fit_cox_model(
          data = df_cox,
          time_var = "target_time",
          status_var = "target_status",
          exposure_var = "sri_q",
          covariates = covariates2
        )
        
        res_cox2 <- tidy_fast_sri_cox(fit_cox2) %>%
          mutate(
            Outcome = outcome_name,
            Analysis = "incident",
            Model = "incident_cov2",
            N = cox_summary$N,
            Cases = cox_summary$Cases,
            .before = 1
          )
        
        fit_cox2_trend <- fit_cox_model(
          data = df_cox,
          time_var = "target_time",
          status_var = "target_status",
          exposure_var = "SRI",
          covariates = covariates2
        )
        
        res_cox2_trend <- tidy_sri_trend_cox(fit_cox2_trend) %>%
          mutate(
            Outcome = outcome_name,
            Analysis = "incident",
            Model = "incident_cov2",
            N = cox_summary$N,
            Cases = cox_summary$Cases,
            .before = 1
          )
        
        addWorksheet(wb, "incident_cov1")
        writeData(wb, "incident_cov1", res_cox1)
        sheet_written <- c(sheet_written, "incident_cov1")
        
        addWorksheet(wb, "incident_cov2")
        writeData(wb, "incident_cov2", res_cox2)
        sheet_written <- c(sheet_written, "incident_cov2")
        
        p_trend_inc <- bind_rows(res_cox1_trend, res_cox2_trend)
        addWorksheet(wb, "P_trend_inc")
        writeData(wb, "P_trend_inc", p_trend_inc)
        sheet_written <- c(sheet_written, "P_trend_inc")
        
        all_results_current[["incident_cov1"]] <- res_cox1
        all_results_current[["incident_cov2"]] <- res_cox2
        all_results_current[["P_trend_inc"]] <- p_trend_inc
      }
    }
    
    # ----------------------------------------------------------
    # 6.6 若无任何结果，报错
    # ----------------------------------------------------------
    if (length(sheet_written) == 0) {
      stop("没有任何结果被成功写入工作簿")
    }
    
    # ----------------------------------------------------------
    # 6.7 保存当前结局Excel
    # ----------------------------------------------------------
    output_file <- file.path(output_folder, paste0(outcome_name, "_results.xlsx"))
    saveWorkbook(wb, output_file, overwrite = TRUE)
    
    # 合并当前结局所有结果用于总表
    combined_result <- bind_rows(all_results_current)
    
    list(
      log = data.frame(
        outcome = outcome_name,
        status = "success",
        message = paste0("已保存: ", output_file),
        stringsAsFactors = FALSE
      ),
      result = combined_result
    )
    
  }, error = function(e) {
    list(
      log = data.frame(
        outcome = basename(outcome_file),
        status = "failed",
        message = as.character(e$message),
        stringsAsFactors = FALSE
      ),
      result = NULL
    )
  })
}

# ==============================================================
# 7. 停止并行
# ==============================================================
stopCluster(cl)
cat("并行计算已结束。\n")

# ==============================================================
# 8. 汇总日志
# ==============================================================
result_log_df <- bind_rows(lapply(result_list, function(x) x$log))
print(result_log_df)

write.csv(
  result_log_df,
  file.path(output_folder, "analysis_log.csv"),
  row.names = FALSE
)

# ==============================================================
# 9. 汇总所有成功结果
# ==============================================================
all_results_df <- bind_rows(
  lapply(result_list, function(x) x$result)
)

write.csv(
  all_results_df,
  file.path(output_folder, "all_outcomes_combined_results.csv"),
  row.names = FALSE
)

# 也可额外保存为Excel
wb_all <- createWorkbook()
addWorksheet(wb_all, "all_results")
writeData(wb_all, "all_results", all_results_df)
saveWorkbook(
  wb_all,
  file.path(output_folder, "all_outcomes_combined_results.xlsx"),
  overwrite = TRUE
)

cat("全部结果已汇总保存。\n")
