#LASSO筛选特征metabolites，构建signature评分

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
