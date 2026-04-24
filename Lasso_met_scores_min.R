#LASSO筛选特征metabolites，构建signature评分

#####载入重要的包#####
library(data.table)
library(tibble)
library(tidyverse)
library(caret)
library(readr)
library(dplyr)
library(glmnet)
library(ggplot2)
library(openxlsx) #读取线性回归的显著结果


##### 1. 设置文件路径 导入数据#####
setwd("D:/ukb data/Sleep regularity/Sleep regularity/Met_lasso")

## 1）读取线性回归结果并筛选代谢物 (新增步骤)
lm_result_file <- "D:/ukb data/Sleep regularity/Sleep regularity/Linear_reg/SRI_Omics_All_Results_260126.xlsx" 

cat("正在读取线性回归结果：", lm_result_file, "\n")

# 读取 "Metabo_Model2" sheet (对应 Model 2 的结果)
df_lm_res <- read.xlsx(lm_result_file, sheet = "Metabo_Model2")

# 筛选 FDR < 0.05 的显著代谢物
significant_mets <- df_lm_res %>%
  filter(FDR < 0.05) %>%
  pull(Outcome) # 提取代谢物名称列

cat("Model 2 中 FDR < 0.05 的代谢物数量:", length(significant_mets), "\n")

if(length(significant_mets) == 0) stop("没有显著的代谢物，无法进行 LASSO 分析！")


# 2）导入代谢/蛋白
df_met <- read_tsv("Metabolomic_Processed_260113.tsv",show_col_types = FALSE)

# 取交集：确保线性回归筛出来的名字在原始数据里也能找到
met.vars <- intersect(significant_mets, colnames(df_met))

length(met.vars) #记录数量 met186 pro2920

# 合并
df_met_subset <- df_met %>% select(all_of(c("eid", met.vars)))

## 4）合并
df_merge <- GGIR_selected %>%
  inner_join(df_met_subset, by = "eid") %>%
  inner_join(CovariatesImputed, by = "eid") %>%
  na.omit() # LASSO 不允许有 NA，必须剔除缺失值

# 定义基本对象
id.var <- "eid"
outcome.var <- "SRI"

# 定义协变量列表：使用矫正睡眠时长的model2 的协变量
covariates <- c("age_test", "MVPA", "season", "sex", "race", 
                "tdi", "smk", "alc", "bmi", "fastingtime", 
                "SleepDurationInSpt_AD_T5A5_mn")

# 最终检查
cat("分析样本量 (N):", nrow(df_merge), "\n") #89436
cat("结局变量:", outcome.var, "\n")
cat("纳入 LASSO 的代谢物数量:", length(met.vars), "\n")

####### 2.划分训练集和验证集########
set.seed(1234)

trainIndex <- createDataPartition(df_merge$SRI, p = 0.7, list = FALSE)
derive   <- as.data.frame(df_merge[trainIndex, ])
validate <- as.data.frame(df_merge[-trainIndex, ])

## 1) 标准化 代谢物本身已标准化，跳过
derive_scaled <- derive
validate_scaled <- validate

## 2) LASSO 矩阵构建
formula_str <- paste("~", paste(c(covariates, met.vars), collapse = " + "))

x <- model.matrix(
  as.formula(formula_str), data = derive_scaled
)[, -1]   # 去掉截距

y <- derive_scaled$SRI_scaled

# 检查
cat("LASSO 建模使用的 Y (SRI_scaled) 均值:", mean(y), "SD:", sd(y), "\n")
# 这里的均值应接近0，SD接近1

## 3) 精确计算协变量列数
x_colnames <- colnames(x)
# 找出所有协变量生成的列名（包括哑变量）
cov_dummy_names <- colnames(model.matrix(as.formula(paste("~", paste(covariates, collapse="+"))), data=derive_scaled)[, -1])

## 4) 设置惩罚与否
penalty.factor <- rep(1, length(x_colnames))
# 把协变量的 penalty 设为 0
penalty.factor[x_colnames %in% cov_dummy_names] <- 0

# 检查一下
cat("不惩罚的协变量列数:", sum(penalty.factor == 0), "\n")
cat("参与 LASSO 的代谢物列数:", sum(penalty.factor == 1), "\n")


# ===== 3. NA check before LASSO =====
cat("NA in outcome y:", sum(is.na(y)), "\n")
cat("NA in design matrix x:", sum(is.na(x)), "\n")

set.seed(1234)
cvfit <- cv.glmnet(
  x, y,
  alpha = 1,          # LASSO
  family = "gaussian",
  penalty.factor = penalty.factor,
  nfolds = 10
)

best.lambda <- cvfit$lambda.min

#提取系数 只提取蛋白系数，不提取截距，和可能加入的协变量的系数
coef_mat <- coef(cvfit, s = best.lambda)

coef_df <- data.frame(
  metabolite = rownames(coef_mat),
  beta = as.numeric(coef_mat)
) %>%
  filter(metabolite %in% met.vars, beta != 0)

#计算代谢组学评分 TIPS：这里计算出的 score 是为了预测标准化的 SRI
# Derivation score
derive$SRI_met_score <- as.matrix(derive_scaled[, coef_df$metabolite]) %*% coef_df$beta

# Validation score
validate$SRI_met_score <- as.matrix(validate_scaled[, coef_df$metabolite]) %*% coef_df$beta

####### 4. 全队列评分##################
derive_out <- derive %>%
  select(eid, SRI) %>%
  mutate(
    source = "derivation",
    SRI_met_score = as.numeric(derive$SRI_met_score)
  )

validate_out <- validate %>%
  select(eid, SRI) %>%
  mutate(
    source = "validation",
    SRI_met_score = as.numeric(validate$SRI_met_score)
  )

final_scores <- bind_rows(derive_out, validate_out)

#输出tsv文件
fwrite(
  final_scores,
  file = "SRI_metabolomic_signature_scores.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Number of selected metabolites:", nrow(coef_df))

######## 5.保存LASSO筛选出的代谢物及系数 ##################
coef_df_out <- coef_df %>%
  arrange(desc(abs(beta)))

fwrite(
  coef_df_out,
  file = "SRI_metabolomic_signature_coefficients.tsv",
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

coef_df_met <- read_tsv(
  "SRI_metabolomic_signature_coefficients.tsv",
  show_col_types = FALSE
)

## 动态计算图片高度
# 统计代谢物个数
n_mets <- nrow(coef_df_met)

# 设定高度策略：
# 基础高度 2 英寸 (给标题和坐标轴留空)
# 每个代谢物分配 0.15 英寸的高度 (可以根据需要调整这个系数，0.15-0.2 比较合适)
# 设置最小高度为 5，防止只有几个代谢物时图太扁
plot_height <- max(4, 1 + n_mets * 0.1)

cat("检测到", n_mets, "个代谢物，将图片高度设定为:", plot_height, "英寸\n")

#作图
coef_plot <- ggplot(
  coef_df_met, 
  aes(
    x = reorder(metabolite, beta), 
    y = beta, 
    fill = beta
  )
) +
  geom_bar(stat = "identity", width = 0.7) + # 调整柱子宽度，0.7 比较美观
  coord_flip() +  # 横向条形图
  scale_fill_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0
  ) +  # 渐变色：负值蓝色，零白色，正值红色
  labs(
    title = "Selected metabolites and Coefficients", 
    x = "Metabolite", 
    y = "Coefficient Value"
  ) +
  theme_minimal() +
  theme(legend.position = "none", # 移除图例
        axis.text.y = element_text(size = 7))  # 调整y轴变量名大小

ggsave(
  filename = "SRI_metabolic_signature_coefficients.pdf",
  plot = coef_plot,
  width = 7,
  height = plot_height, # 高度动态
  units = "in",
  limitsize = FALSE # 【关键】允许生成超长图片 (防止报错)
)

################################# 6.CV R-square 分布检查 #######################
fit_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5,
  returnResamp = "all"
) #把CV模型结果储存下来

set.seed(1234)

lasso_caret <- caret::train(
  x = x,
  y = y,
  method = "glmnet",
  trControl = fit_control,
  tuneGrid = expand.grid(
    alpha  = 1,                # 影子LASSO
    lambda = cvfit$lambda.min      # 用已经跑过的 lambda 网格
  ),
  metric = "Rsquared",
  penalty.factor = penalty.factor
)

r2_df <- tibble(
  outcome = "SRI",
  r2 = lasso_caret$resample$Rsquared
)

ggplot(r2_df, aes(y = outcome, x = r2)) +
  geom_boxplot() +
  theme_bw() +
  ylab(expression(R^2)) +
  xlab("") +
  labs(
    title = "LASSO model fit",
    subtitle = "R2 from hold-out folds in cross-validation"
  )

ggsave(
  filename = "SRI_met_lasso_r2_cvHoldOuts.pdf",
  device = "pdf",
  height = 4,
  width = 6
)

