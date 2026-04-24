library(data.table)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# 如果暴露是 baseline-to-repeat 的 change pattern，那么死亡结局应从 repeat measurement 开始算，而不是从 baseline 开始。

# 1. 把有重复测量的患者全部筛出来
duplicate_patients <- GGIR_selected %>%
  group_by(patient_id) %>%
  filter(n() > 1) %>%
  arrange(patient_id, calendar_date, ID) %>% 
  ungroup()

duplicate_patients_selected <- duplicate_patients %>%
  dplyr::select(
    patient_id,
    calendar_date,
    SleepDurationInSpt_AD_T5A5_mn,
    SleepRegularityIndex_AD_T5A5_mn
  )

patient_sleep_stats <- duplicate_patients %>%
  dplyr::group_by(patient_id) %>%
  dplyr::summarise(
    mean_sleep_duration = mean(SleepDurationInSpt_AD_T5A5_mn, na.rm = TRUE),
    sd_sleep_duration   = sd(SleepDurationInSpt_AD_T5A5_mn, na.rm = TRUE),
    mean_sleep_regularity = mean(SleepRegularityIndex_AD_T5A5_mn, na.rm = TRUE),
    sd_sleep_regularity   = sd(SleepRegularityIndex_AD_T5A5_mn, na.rm = TRUE),
    n_measurements = dplyr::n()
  ) %>%
  dplyr::arrange(patient_id)

sleep_summary <- duplicate_patients_selected %>%
  dplyr::group_by(patient_id) %>%
  dplyr::arrange(calendar_date, .by_group = TRUE) %>%
  dplyr::mutate(measure_order = dplyr::row_number()) %>%
  dplyr::summarise(
    calendar_date_first = dplyr::first(calendar_date),
    calendar_date_last = dplyr::last(calendar_date),
    days_between_first_last = as.numeric(calendar_date_last - calendar_date_first),
    
    SRI_baseline = dplyr::first(SleepRegularityIndex_AD_T5A5_mn),
    duration_baseline = dplyr::first(SleepDurationInSpt_AD_T5A5_mn),
    
    SRI_repeat = {
      x <- SleepRegularityIndex_AD_T5A5_mn[measure_order > 1]
      if (length(x) == 0 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    },
    duration_repeat = {
      x <- SleepDurationInSpt_AD_T5A5_mn[measure_order > 1]
      if (length(x) == 0 || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    },
    .groups = "drop"
  )

sri_q <- quantile(
  GGIR_selected$SleepRegularityIndex_AD_T5A5_mn,
  probs = c(0.25, 0.50, 0.75),
  na.rm = TRUE
)

sleep_summary <- sleep_summary %>%
  dplyr::mutate(
    SRI_baseline_Q = dplyr::case_when(
      is.na(SRI_baseline) ~ NA_integer_,
      SRI_baseline <= sri_q[1] ~ 1L,
      SRI_baseline <= sri_q[2] ~ 2L,
      SRI_baseline <= sri_q[3] ~ 3L,
      SRI_baseline >  sri_q[3] ~ 4L
    ),
    SRI_repeat_Q = dplyr::case_when(
      is.na(SRI_repeat) ~ NA_integer_,
      SRI_repeat <= sri_q[1] ~ 1L,
      SRI_repeat <= sri_q[2] ~ 2L,
      SRI_repeat <= sri_q[3] ~ 3L,
      SRI_repeat >  sri_q[3] ~ 4L
    )
  ) %>%
  dplyr::filter(!is.na(SRI_repeat_Q))

table(sleep_summary$SRI_baseline_Q, useNA = "ifany")
table(sleep_summary$SRI_repeat_Q, useNA = "ifany")


########## 构建SRI_change_pattern ##############
sleep_summary <- sleep_summary %>%
  dplyr::mutate(
    SRI_baseline_subtype = dplyr::case_when(
      is.na(SRI_baseline_Q) ~ NA_character_,
      SRI_baseline_Q == 1L ~ "Irregular",
      SRI_baseline_Q %in% c(2L, 3L) ~ "Medium",
      SRI_baseline_Q == 4L ~ "Regular"
    ),
    SRI_repeat_subtype = dplyr::case_when(
      is.na(SRI_repeat_Q) ~ NA_character_,
      SRI_repeat_Q == 1L ~ "Irregular",
      SRI_repeat_Q %in% c(2L, 3L) ~ "Medium",
      SRI_repeat_Q == 4L ~ "Regular"
    ),
    SRI_baseline_subtype = factor(
      SRI_baseline_subtype,
      levels = c("Regular", "Medium", "Irregular")
    ),
    SRI_repeat_subtype = factor(
      SRI_repeat_subtype,
      levels = c("Regular", "Medium", "Irregular")
    )
  )
library(dplyr)
library(ggplot2)
library(ggalluvial)

# 1. 构建 change pattern 数据
sankey_df <- sleep_summary %>%
  dplyr::filter(
    !is.na(SRI_baseline_subtype),
    !is.na(SRI_repeat_subtype)
  ) %>%
  dplyr::count(SRI_baseline_subtype, SRI_repeat_subtype, name = "Freq") %>%
  dplyr::mutate(
    change_pattern = paste0(SRI_baseline_subtype, " \u2192 ", SRI_repeat_subtype)
  )

# 2. 查看各转变模式人数
pattern_table <- sleep_summary %>%
  dplyr::filter(
    !is.na(SRI_baseline_subtype),
    !is.na(SRI_repeat_subtype)
  ) %>%
  dplyr::mutate(
    change_pattern = paste0(SRI_baseline_subtype, " \u2192 ", SRI_repeat_subtype)
  ) %>%
  dplyr::count(change_pattern, sort = TRUE)

print(pattern_table)

sleep_summary <- sleep_summary %>%
  dplyr::mutate(
    SRI_baseline_score = dplyr::case_when(
      SRI_baseline_subtype == "Irregular" ~ 1L,
      SRI_baseline_subtype == "Medium" ~ 2L,
      SRI_baseline_subtype == "Regular" ~ 3L,
      TRUE ~ NA_integer_
    ),
    SRI_repeat_score = dplyr::case_when(
      SRI_repeat_subtype == "Irregular" ~ 1L,
      SRI_repeat_subtype == "Medium" ~ 2L,
      SRI_repeat_subtype == "Regular" ~ 3L,
      TRUE ~ NA_integer_
    ),
    SRI_change_5cat = dplyr::case_when(
      is.na(SRI_baseline_score) | is.na(SRI_repeat_score) ~ NA_character_,
      SRI_baseline_subtype == "Regular" & SRI_repeat_subtype == "Regular" ~ "Stable regular",
      SRI_repeat_subtype == "Regular" & SRI_repeat_score > SRI_baseline_score ~ "Improved to regular",
      SRI_baseline_subtype == "Medium" & SRI_repeat_subtype == "Medium" ~ "Stable medium",
      SRI_repeat_subtype == "Irregular" & SRI_repeat_score < SRI_baseline_score ~ "Worsened to irregular",
      SRI_baseline_subtype == "Irregular" & SRI_repeat_subtype == "Irregular" ~ "Stable irregular",
      TRUE ~ NA_character_
    ),
    SRI_change_5cat = factor(
      SRI_change_5cat,levels = c(
        "Stable irregular",
        "Stable medium",
        "Stable regular",
        "Worsened to irregular",
        "Improved to regular"
      )
    )
  )

####### 计算全部15个结局 #######
# 1. 获取所有 outcome 文件
outcome_folder_name <- "~/hzq_ukb_SRI/15_outcomes"
outcome_files <- list.files(
  path = outcome_folder_name,
  pattern = "\\.tsv$",
  full.names = TRUE
)

# 2. 设定 change pattern 顺序
sleep_summary <- sleep_summary %>%
  dplyr::mutate(
    SRI_change_5cat = factor(
      SRI_change_5cat,
      levels = c(
        "Stable irregular",
        "Stable medium",
        "Stable regular",
        "Worsened to irregular",
        "Improved to regular"
      )
    )
  )

# 3. 定义单个 outcome 的发病率计算函数
run_single_incidence_5cat <- function(outcome_file) {
  
  outcome_name <- basename(outcome_file) %>%
    stringr::str_remove("_outcome_\\d+\\.tsv$")
  
  outcome_df <- readr::read_tsv(outcome_file, show_col_types = FALSE) %>%
    dplyr::select(eid, target_status, target_time)
  
  res <- sleep_summary %>%
    dplyr::select(eid, SRI_change_5cat, days_between_first_last) %>%
    dplyr::inner_join(outcome_df, by = "eid") %>%
    dplyr::mutate(
      target_time_after_change = target_time
    ) %>%
    dplyr::filter(
      !is.na(SRI_change_5cat),
      !is.na(days_between_first_last),
      !is.na(target_status),
      !is.na(target_time_after_change),
      target_time_after_change > 0
    ) %>%
    dplyr::group_by(SRI_change_5cat) %>%
    dplyr::summarise(
      N = dplyr::n(),
      Cases = sum(target_status == 1, na.rm = TRUE),
      Risk = Cases / N,
      Risk_percent = 100 * Risk,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      outcome = outcome_name,
      Risk_percent = round(Risk_percent, 2),
      label = sprintf("%.2f%%", Risk_percent)
    )
  
  return(res)
}

# 4. 循环汇总所有结局
incidence_5cat_all <- purrr::map_dfr(outcome_files, run_single_incidence_5cat)

# 5. 如果想固定 outcome 显示顺序，可按文件顺序设置
outcome_order <- basename(outcome_files) %>%
  stringr::str_remove("_outcome_\\d+\\.tsv$")

incidence_5cat_all <- incidence_5cat_all %>%
  dplyr::mutate(
    outcome = factor(outcome, levels = outcome_order)
  )

# 6. 查看结果表
print(incidence_5cat_all)