#!/usr/bin/env Rscript

# Template for estimating incremental explained variance from a feature block.
# The script compares a base linear model with a full model that adds candidate
# features, then reports R2, adjusted R2, and the partial F-test.

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    analysis_file = NA_character_,
    feature_list_file = NA_character_,
    outdir = NA_character_,
    id_col = "id",
    outcome_col = NA_character_,
    base_covariates = NA_character_,
    feature_cols = NA_character_,
    feature_name_col = "feature",
    factor_covariates = "",
    missingness_cutoff = 0.30,
    impute_features = TRUE,
    impute_m = 1,
    impute_maxit = 5,
    seed = 42,
    output_prefix = "incremental_r2",
    verbose = FALSE
  )

  bool_value <- function(x) {
    x <- tolower(x)
    if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
    if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
    stop("Invalid logical value: ", x)
  }

  for (arg in args) {
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", kv[1])
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "TRUE"

    if (!key %in% names(opts)) stop("Unknown option: --", key)
    if (is.logical(opts[[key]])) {
      opts[[key]] <- bool_value(val)
    } else if (is.numeric(opts[[key]])) {
      opts[[key]] <- as.numeric(val)
    } else {
      opts[[key]] <- val
    }
  }

  required <- c("analysis_file", "outdir", "outcome_col", "base_covariates")
  missing_required <- required[is.na(unlist(opts[required])) | !nzchar(unlist(opts[required]))]
  if (length(missing_required) > 0) {
    stop("Missing required option(s): --", paste(missing_required, collapse = ", --"))
  }

  if ((is.na(opts$feature_cols) || !nzchar(opts$feature_cols)) &&
      (is.na(opts$feature_list_file) || !nzchar(opts$feature_list_file))) {
    stop("Set either --feature_cols or --feature_list_file.")
  }

  opts$analysis_file <- normalizePath(opts$analysis_file, winslash = "/", mustWork = TRUE)
  if (!is.na(opts$feature_list_file) && nzchar(opts$feature_list_file)) {
    opts$feature_list_file <- normalizePath(opts$feature_list_file, winslash = "/", mustWork = TRUE)
  }
  opts$outdir <- normalizePath(opts$outdir, winslash = "/", mustWork = FALSE)
  opts
}

split_names <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

bt <- function(x) {
  paste0("`", gsub("`", "``", x), "`")
}

read_table_auto <- function(path, select = NULL) {
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  data.table::fread(path, sep = sep, select = select, showProgress = FALSE)
}

read_feature_names <- function(opts) {
  features_from_args <- split_names(opts$feature_cols)
  if (length(features_from_args) > 0) {
    return(unique(features_from_args))
  }

  feature_dt <- read_table_auto(opts$feature_list_file)
  if (!opts$feature_name_col %in% names(feature_dt)) {
    stop("Feature list does not contain column: ", opts$feature_name_col)
  }

  features <- as.character(feature_dt[[opts$feature_name_col]])
  features <- features[!is.na(features) & nzchar(features)]
  unique(features)
}

drop_invalid_features <- function(dt, vars) {
  vars[
    vapply(vars, function(v) {
      x <- dt[[v]]
      x <- x[!is.na(x)]
      length(unique(x)) >= 2
    }, logical(1))
  ]
}

make_feature_dictionary <- function(feature_names, analysis_data, missingness_cutoff) {
  data.table(
    feature = feature_names,
    available_in_analysis_data = feature_names %in% names(analysis_data)
  )[
    ,
    `:=`(
      n_rows = ifelse(
        available_in_analysis_data,
        nrow(analysis_data),
        NA_integer_
      ),
      n_nonmissing = vapply(feature, function(v) {
        if (!v %in% names(analysis_data)) return(NA_integer_)
        sum(!is.na(analysis_data[[v]]))
      }, integer(1)),
      n_missing = vapply(feature, function(v) {
        if (!v %in% names(analysis_data)) return(NA_integer_)
        sum(is.na(analysis_data[[v]]))
      }, integer(1)),
      missing_pct = vapply(feature, function(v) {
        if (!v %in% names(analysis_data)) return(NA_real_)
        mean(is.na(analysis_data[[v]]))
      }, numeric(1))
    )
  ][
    ,
    dropped_reason := fifelse(
      !available_in_analysis_data,
      "missing_column",
      fifelse(
        missing_pct > missingness_cutoff,
        "high_missingness",
        "candidate"
      )
    )
  ]
}

prepare_model_data <- function(analysis_data, outcome_col, base_covariates,
                               factor_covariates, candidate_features,
                               impute_features, impute_m, impute_maxit, seed) {
  model_vars <- unique(c(outcome_col, base_covariates, candidate_features))
  model_data <- copy(analysis_data[, ..model_vars])

  model_data <- model_data[!is.na(get(outcome_col))]
  for (v in base_covariates) {
    model_data <- model_data[!is.na(get(v))]
  }

  for (v in factor_covariates) {
    if (v %in% names(model_data)) {
      model_data[, (v) := droplevels(as.factor(get(v)))]
    }
  }

  for (v in setdiff(candidate_features, factor_covariates)) {
    model_data[, (v) := suppressWarnings(as.numeric(get(v)))]
  }

  candidate_features <- drop_invalid_features(model_data, candidate_features)
  model_vars <- unique(c(outcome_col, base_covariates, candidate_features))
  model_data <- model_data[, ..model_vars]

  if (impute_features && length(candidate_features) > 0) {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("Package 'mice' is required when --impute_features=TRUE.")
    }

    set.seed(seed)
    mice_methods <- mice::make.method(model_data)
    mice_methods[] <- ""
    mice_methods[candidate_features] <- "pmm"

    predictor_matrix <- mice::make.predictorMatrix(model_data)
    predictor_matrix[,] <- 0
    predictor_matrix[candidate_features, model_vars] <- 1
    diag(predictor_matrix) <- 0

    mice_fit <- mice::mice(
      model_data,
      m = impute_m,
      maxit = impute_maxit,
      method = mice_methods,
      predictorMatrix = predictor_matrix,
      printFlag = FALSE,
      seed = seed
    )
    model_data <- as.data.table(mice::complete(mice_fit, 1))
  } else {
    model_data <- model_data[complete.cases(model_data)]
  }

  for (v in factor_covariates) {
    if (v %in% names(model_data)) {
      model_data[, (v) := droplevels(as.factor(get(v)))]
    }
  }

  valid_features <- drop_invalid_features(model_data, candidate_features)
  model_data <- model_data[
    complete.cases(model_data[, c(outcome_col, base_covariates, valid_features), with = FALSE])
  ]

  list(
    data = model_data,
    candidate_features = candidate_features,
    valid_features = valid_features
  )
}

fit_incremental_r2 <- function(model_data, outcome_col, base_covariates, valid_features) {
  base_formula_text <- paste(bt(outcome_col), "~", paste(bt(base_covariates), collapse = " + "))
  full_terms <- c(base_covariates, valid_features)
  full_formula_text <- paste(bt(outcome_col), "~", paste(bt(full_terms), collapse = " + "))

  fit_base <- lm(as.formula(base_formula_text), data = model_data)
  fit_full <- lm(as.formula(full_formula_text), data = model_data)

  base_summary <- summary(fit_base)
  full_summary <- summary(fit_full)
  anova_base_full <- anova(fit_base, fit_full)

  data.table(
    outcome = outcome_col,
    base_model = base_formula_text,
    full_model = paste(bt(outcome_col), "~ base covariates + feature block"),
    n_complete_case = nrow(model_data),
    n_base_covariates = length(base_covariates),
    n_features_used = length(valid_features),
    r2_base = base_summary$r.squared,
    r2_full = full_summary$r.squared,
    delta_r2 = full_summary$r.squared - base_summary$r.squared,
    delta_r2_percent = 100 * (full_summary$r.squared - base_summary$r.squared),
    adj_r2_base = base_summary$adj.r.squared,
    adj_r2_full = full_summary$adj.r.squared,
    delta_adj_r2 = full_summary$adj.r.squared - base_summary$adj.r.squared,
    delta_adj_r2_percent = 100 * (full_summary$adj.r.squared - base_summary$adj.r.squared),
    partial_f_p = anova_base_full$`Pr(>F)`[2],
    note = paste(
      "Incremental explained variance was estimated as full-model R2 minus base-model R2.",
      "Candidate features above the missingness cutoff or without variation were excluded.",
      "Remaining feature missingness was handled according to the imputation setting."
    )
  )
}

main <- function() {
  opts <- parse_args()
  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

  base_covariates <- split_names(opts$base_covariates)
  factor_covariates <- split_names(opts$factor_covariates)
  feature_names <- read_feature_names(opts)

  analysis_data <- read_table_auto(opts$analysis_file)
  required_cols <- unique(c(opts$id_col, opts$outcome_col, base_covariates))
  missing_required <- setdiff(required_cols, names(analysis_data))
  if (length(missing_required) > 0) {
    stop("Analysis file is missing required column(s): ", paste(missing_required, collapse = ", "))
  }

  if (anyDuplicated(analysis_data[[opts$id_col]]) && opts$verbose) {
    warning("Duplicate IDs were detected in the analysis data.")
  }

  feature_dictionary <- make_feature_dictionary(
    feature_names = feature_names,
    analysis_data = analysis_data,
    missingness_cutoff = opts$missingness_cutoff
  )

  candidate_features <- feature_dictionary[
    dropped_reason == "candidate",
    feature
  ]
  if (length(candidate_features) == 0) {
    stop("No candidate features remained after availability and missingness filtering.")
  }

  prepared <- prepare_model_data(
    analysis_data = analysis_data,
    outcome_col = opts$outcome_col,
    base_covariates = base_covariates,
    factor_covariates = factor_covariates,
    candidate_features = candidate_features,
    impute_features = opts$impute_features,
    impute_m = opts$impute_m,
    impute_maxit = opts$impute_maxit,
    seed = opts$seed
  )
  if (length(prepared$valid_features) == 0) {
    stop("No valid candidate features remained after model-data preparation.")
  }

  feature_dictionary[
    feature %in% setdiff(prepared$candidate_features, prepared$valid_features),
    dropped_reason := "no_variation_after_preparation"
  ]
  feature_dictionary[
    feature %in% prepared$valid_features,
    dropped_reason := "used"
  ]
  feature_dictionary[, used_in_full_model := dropped_reason == "used"]

  summary_dt <- fit_incremental_r2(
    model_data = prepared$data,
    outcome_col = opts$outcome_col,
    base_covariates = base_covariates,
    valid_features = prepared$valid_features
  )

  data.table::fwrite(
    summary_dt,
    file.path(opts$outdir, paste0(opts$output_prefix, "_summary.csv"))
  )
  data.table::fwrite(
    feature_dictionary,
    file.path(opts$outdir, paste0(opts$output_prefix, "_feature_dictionary.csv"))
  )

  message("Incremental R2 analysis completed.")
  if (opts$verbose) {
    print(summary_dt)
  }
}

main()
