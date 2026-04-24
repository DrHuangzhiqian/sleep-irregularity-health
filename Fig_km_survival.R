# ==============================================================================
# 0. 加载包
# ==============================================================================
library(data.table)
library(dplyr)
library(readr)
library(tidyr)
library(survival)
library(survminer)
library(ggplot2)
library(tibble)
library(scales)
# ==============================================================================
# 4. 合并基础数据并构造 SRI 四分位
# ==============================================================================
df_base <- GGIR_selected %>%
  inner_join(CovariatesImputed, by = "eid") %>%
  mutate(
    SRI_Quartile_num = ntile(SleepRegularityIndex_AD_T5A5_mn, 4),
    SRI_Quartile = factor(
      SRI_Quartile_num,
      levels = c(4, 3, 2, 1),
      labels = c("Q4", "Q3", "Q2", "Q1")
    )
  ) %>%
  as.data.frame()

# ==============================================================================
# 5. 读取结局数据
# ==============================================================================
outcome_file <- "16_outcomes/death_outcome_260116.tsv"

outcome_df <- read_tsv(outcome_file, show_col_types = FALSE) %>%
  dplyr::select(eid, target_status, target_time) %>%
  mutate(
    target_status = as.numeric(target_status),
    target_time = as.numeric(target_time)
  )

# ==============================================================================
# 6. 构建分析数据集
# 如果只是画 KM 曲线，理论上只需要 drop_na(time_months, target_status, SRI_Quartile)
# 这里保留你原来较严格的版本
# ==============================================================================
df_analysis <- df_base %>%
  inner_join(outcome_df, by = "eid") %>%
  filter(!is.na(target_time), target_time > 0) %>%
  mutate(
    time_years  = target_time / 365.25,
    time_months = target_time / 30.4375
  ) %>%
  drop_na(
    time_years,
    time_months,
    target_status,
    SRI_Quartile,
    age_test,
    sex,
    MVPA,
    season,
    race,
    tdi,
    smk,
    alc,
    bmi,
    fastingtime,
    SleepDurationInSpt_AD_T5A5_mn
  ) %>%
  as.data.frame()

cat("分析样本量: ", nrow(df_analysis), "\n")
cat("结局事件数: ", sum(df_analysis$target_status == 1, na.rm = TRUE), "\n")

cat("\nSRI 四分位分布:\n")
print(table(df_analysis$SRI_Quartile, useNA = "ifany"))

# ==============================================================================
# 7. 仅保留 Q1 和 Q4 用于生存曲线
# ==============================================================================
df_surv <- df_analysis %>%
  filter(SRI_Quartile %in% c("Q1", "Q4")) %>%
  mutate(
    SRI_Quartile = factor(SRI_Quartile, levels = c("Q4", "Q1"))
  )

cat("\nQ1/Q4 样本量:\n")
print(table(df_surv$SRI_Quartile))

cat("\nQ1/Q4 事件数:\n")
print(table(df_surv$SRI_Quartile, df_surv$target_status))

cat("\n随访时间（月）概况:\n")
print(summary(df_surv$time_months))

# ==============================================================================
# 8. Kaplan-Meier 生存拟合
# ==============================================================================
fit_km <- survfit(
  Surv(time_months, target_status) ~ SRI_Quartile,
  data = df_surv
)

# ==============================================================================
# 9. 先计算 log-rank P 值，并改成文章风格文字
# ==============================================================================
survdiff_res <- survdiff(
  Surv(time_months, target_status) ~ SRI_Quartile,
  data = df_surv
)

p_value <- 1 - pchisq(survdiff_res$chisq, df = length(survdiff_res$n) - 1)

p_label <- if (p_value < 0.0001) {
  "Log-rank P < 0.0001"
} else {
  paste0("Log-rank P = ", formatC(p_value, format = "f", digits = 4))
}

# ==============================================================================
# 10. 绘制生存曲线
# ==============================================================================
p <- ggsurvplot(
  fit_km,
  data = df_surv,
  conf.int = TRUE,
  conf.int.alpha = 0.12,
  risk.table = FALSE,
  pval = FALSE,            # 关闭默认p值，改用手动注释
  censor = FALSE,
  palette = c("Q4" = "#2b6cb0", "Q1" = "#c83c3c"),
  legend.title = "",
  legend.labs = c("Q4", "Q1"),
  xlab = "Month",
  ylab = "Survival probability (%)",
  break.time.by = 12,
  ggtheme = theme_classic()
)

# ==============================================================================
# 11. 美化 + 百分比坐标轴 + 自定义 p 值
# ==============================================================================
p$plot <- p$plot +
  annotate(
    "text",
    x = 8,
    y = 0.89,
    label = p_label,
    hjust = 0,
    size = 5
  ) +
  scale_x_continuous(
    limits = c(0, 144),
    breaks = seq(0, 144, 12),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    limits = c(0.875, 1.000),
    breaks = c(0.875, 0.900, 0.925, 0.950, 0.975, 1.000),
    labels = label_percent(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  theme(
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.position = c(0.86, 0.84),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# 显示图
print(p)

# ==============================================================================
# 12. 保存图片
# ==============================================================================
ggsave(
  filename = file.path(out_dir, "KM_survival_Q1_vs_Q4_month_CI_LogRankP_percent_87.5to100.png"),
  plot = p$plot,
  width = 8,
  height = 6.5,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "KM_survival_Q1_vs_Q4_month_CI_LogRankP_percent_87.5to100.pdf"),
  plot = p$plot,
  width = 8,
  height = 6.5
)

cat("\n生存曲线已保存到: ", out_dir, "\n")