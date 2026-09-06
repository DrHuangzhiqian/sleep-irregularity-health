#!/usr/bin/env Rscript

# Generic template for high-dimensional mediation screening using HIMA.
# This script is intentionally data-agnostic: all file paths and variable names
# are provided at runtime.

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(HIMA)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    pheno_file = NA_character_,
    mediator_file = NA_character_,
    mediator_list_file = NA_character_,
    outdir = NA_character_,
    id_col = "id",
    exposure_col = NA_character_,
    outcome_col = NA_character_,
    time_col = "",
    status_col = "",
    covariates = "",
    factor_covariates = "",
    mediator_cols = "",
    mediator_name_col = "mediator",
    y_type = "continuous",
    mediator_type = "gaussian",
    penalty = "",
    sigcut = 0.05,
    min_n = 50,
    min_events = 10,
    standardize_exposure = TRUE,
    standardize_mediators = TRUE,
    ncore = 1,
    output_prefix = "hima_mediation",
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

  required <- c("pheno_file", "outdir", "exposure_col", "outcome_col")
  missing_required <- required[is.na(unlist(opts[required])) | !nzchar(unlist(opts[required]))]
  if (length(missing_required) > 0) {
    stop("Missing required option(s): --", paste(missing_required, collapse = ", --"))
  }

  if (!opts$y_type %in% c("continuous", "binary", "survival")) {
    stop("--y_type must be continuous, binary, or survival.")
  }
  if (opts$y_type == "survival" && (!nzchar(opts$time_col) || !nzchar(opts$status_col))) {
    stop("--time_col and --status_col are required when --y_type=survival.")
  }
  if (opts$y_type != "survival" && nzchar(opts$status_col)) {
    opts$outcome_col <- opts$status_col
  }

  opts$pheno_file <- normalizePath(opts$pheno_file, winslash = "/", mustWork = TRUE)
  if (!is.na(opts$mediator_file) && nzchar(opts$mediator_file)) {
    opts$mediator_file <- normalizePath(opts$mediator_file, winslash = "/", mustWork = TRUE)
  }
  if (!is.na(opts$mediator_list_file) && nzchar(opts$mediator_list_file)) {
    opts$mediator_list_file <- normalizePath(opts$mediator_list_file, winslash = "/", mustWork = TRUE)
  }
  opts$outdir <- normalizePath(opts$outdir, winslash = "/", mustWork = FALSE)
  opts
}

split_names <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

read_table_auto <- function(path, select = NULL) {
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  data.table::fread(path, sep = sep, select = select, showProgress = FALSE)
}

as_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

nonconstant_columns <- function(mat) {
  apply(mat, 2, function(x) {
    x <- x[!is.na(x)]
    length(unique(x)) >= 2
  })
}

read_mediator_names <- function(opts, available_names) {
  from_args <- split_names(opts$mediator_cols)
  if (length(from_args) > 0) {
    return(unique(from_args))
  }

  if (!is.na(opts$mediator_list_file) && nzchar(opts$mediator_list_file)) {
    mediator_list <- read_table_auto(opts$mediator_list_file)
    if (!opts$mediator_name_col %in% names(mediator_list)) {
      stop("Mediator list does not contain column: ", opts$mediator_name_col)
    }
    mediators <- as.character(mediator_list[[opts$mediator_name_col]])
    mediators <- mediators[!is.na(mediators) & nzchar(mediators)]
    return(unique(mediators))
  }

  reserved <- c(
    opts$id_col,
    opts$exposure_col,
    opts$outcome_col,
    opts$time_col,
    opts$status_col,
    split_names(opts$covariates)
  )
  setdiff(available_names, reserved[nzchar(reserved)])
}

prepare_data <- function(opts) {
  covariates <- split_names(opts$covariates)
  factor_covariates <- split_names(opts$factor_covariates)

  pheno <- read_table_auto(opts$pheno_file)
  required_pheno <- unique(c(opts$id_col, opts$exposure_col, opts$outcome_col, opts$time_col, opts$status_col, covariates))
  required_pheno <- required_pheno[nzchar(required_pheno)]
  missing_pheno <- setdiff(required_pheno, names(pheno))
  if (length(missing_pheno) > 0) {
    stop("Phenotype file is missing column(s): ", paste(missing_pheno, collapse = ", "))
  }

  pheno[, (opts$id_col) := as.character(get(opts$id_col))]

  if (!is.na(opts$mediator_file) && nzchar(opts$mediator_file)) {
    mediator_dt <- read_table_auto(opts$mediator_file)
    if (!opts$id_col %in% names(mediator_dt)) {
      stop("Mediator file is missing ID column: ", opts$id_col)
    }
    mediator_dt[, (opts$id_col) := as.character(get(opts$id_col))]
    dat <- merge(pheno, mediator_dt, by = opts$id_col, all = FALSE, sort = FALSE)
  } else {
    dat <- copy(pheno)
  }

  mediator_names <- read_mediator_names(opts, names(dat))
  missing_mediators <- setdiff(mediator_names, names(dat))
  if (length(missing_mediators) > 0) {
    warning("Dropping missing mediator column(s): ", paste(missing_mediators, collapse = ", "))
  }
  mediator_names <- intersect(mediator_names, names(dat))
  if (length(mediator_names) == 0) {
    stop("No mediator columns were available.")
  }

  dat[, (opts$exposure_col) := as_num(get(opts$exposure_col))]
  if (opts$standardize_exposure) {
    dat[, (opts$exposure_col) := as.numeric(scale(get(opts$exposure_col)))]
  }

  if (opts$y_type == "survival") {
    dat[, (opts$time_col) := as_num(get(opts$time_col))]
    dat[, (opts$status_col) := as.integer(get(opts$status_col))]
  } else if (opts$y_type == "binary") {
    dat[, (opts$outcome_col) := as.integer(get(opts$outcome_col))]
  } else {
    dat[, (opts$outcome_col) := as_num(get(opts$outcome_col))]
  }

  for (v in factor_covariates) {
    if (v %in% names(dat)) {
      dat[, (v) := droplevels(as.factor(get(v)))]
    }
  }

  for (v in setdiff(covariates, factor_covariates)) {
    if (v %in% names(dat)) {
      dat[, (v) := as_num(get(v))]
    }
  }

  for (v in mediator_names) {
    dat[, (v) := as_num(get(v))]
  }

  if (opts$standardize_mediators) {
    for (v in mediator_names) {
      dat[, (v) := as.numeric(scale(get(v)))]
    }
  }

  needed <- unique(c(opts$exposure_col, opts$outcome_col, opts$time_col, opts$status_col, covariates, mediator_names))
  needed <- needed[nzchar(needed)]
  dat <- dat[complete.cases(dat[, ..needed])]

  if (nrow(dat) < opts$min_n) {
    stop("Insufficient complete cases after filtering.")
  }
  if (opts$y_type %in% c("binary", "survival")) {
    n_events <- sum(dat[[if (opts$y_type == "survival") opts$status_col else opts$outcome_col]] == 1, na.rm = TRUE)
    if (n_events < opts$min_events) {
      stop("Insufficient outcome events after filtering.")
    }
  }

  list(
    data = dat,
    covariates = covariates,
    factor_covariates = factor_covariates,
    mediator_names = mediator_names
  )
}

make_hima_inputs <- function(prepared, opts) {
  dat <- prepared$data
  covariates <- prepared$covariates
  mediator_names <- prepared$mediator_names

  mediator_matrix <- as.matrix(dat[, ..mediator_names])
  storage.mode(mediator_matrix) <- "double"
  mediator_keep <- nonconstant_columns(mediator_matrix)
  mediator_matrix <- mediator_matrix[, mediator_keep, drop = FALSE]
  mediator_map <- data.table(
    mediator = mediator_names[mediator_keep],
    hima_id = make.names(mediator_names[mediator_keep], unique = TRUE)
  )
  colnames(mediator_matrix) <- mediator_map$hima_id
  if (ncol(mediator_matrix) == 0) {
    stop("No nonconstant mediator columns remained.")
  }

  data_pheno <- data.frame(hima_exposure = dat[[opts$exposure_col]], check.names = FALSE)
  if (opts$y_type == "survival") {
    data_pheno$hima_time <- dat[[opts$time_col]]
    data_pheno$hima_status <- dat[[opts$status_col]]
  } else {
    data_pheno$hima_outcome <- dat[[opts$outcome_col]]
  }

  covar_map <- data.table(hima_covariate = character(), original_covariate = character())
  covar_names <- character()
  if (length(covariates) > 0) {
    covar_formula <- as.formula(paste("~", paste(covariates, collapse = " + ")))
    covar_matrix <- model.matrix(covar_formula, data = dat)
    covar_matrix <- covar_matrix[, colnames(covar_matrix) != "(Intercept)", drop = FALSE]
    covar_keep <- nonconstant_columns(covar_matrix)
    original_covariates <- colnames(covar_matrix)[covar_keep]
    covar_matrix <- covar_matrix[, covar_keep, drop = FALSE]
    covar_names <- make.names(original_covariates, unique = TRUE)
    colnames(covar_matrix) <- covar_names
    covar_map <- data.table(
      hima_covariate = covar_names,
      original_covariate = original_covariates
    )
    data_pheno <- cbind(data_pheno, as.data.frame(covar_matrix, check.names = FALSE))
  }

  rhs <- c("hima_exposure", covar_names)
  formula <- if (opts$y_type == "survival") {
    as.formula(paste("Surv(hima_time, hima_status) ~", paste(rhs, collapse = " + ")))
  } else {
    as.formula(paste("hima_outcome ~", paste(rhs, collapse = " + ")))
  }

  list(
    formula = formula,
    data_pheno = data_pheno,
    mediator_matrix = mediator_matrix,
    mediator_map = mediator_map,
    covariate_map = covar_map
  )
}

hima_to_table <- function(fit, mediator_map) {
  if (is.null(fit) || is.null(fit$ID) || length(fit$ID) == 0) {
    return(data.table())
  }

  res <- data.table(
    hima_id = fit$ID,
    alpha = as_num(fit$alpha),
    beta = as_num(fit$beta),
    indirect_effect = as_num(fit[["alpha*beta"]]),
    relative_importance_pct = as_num(fit$rimp),
    p_value = as_num(fit[["p-value"]])
  )
  merge(res, mediator_map, by = "hima_id", all.x = TRUE, sort = FALSE)
}

main <- function() {
  opts <- parse_args()
  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

  prepared <- prepare_data(opts)
  hima_inputs <- make_hima_inputs(prepared, opts)

  penalty <- opts$penalty
  if (!nzchar(penalty)) {
    penalty <- if (opts$y_type == "survival") "DBlasso" else "MCP"
  }

  fit <- tryCatch(
    HIMA::hima(
      formula = hima_inputs$formula,
      data.pheno = hima_inputs$data_pheno,
      data.M = hima_inputs$mediator_matrix,
      mediator.type = opts$mediator_type,
      penalty = penalty,
      scale = FALSE,
      sigcut = opts$sigcut,
      Y.type = opts$y_type,
      parallel = opts$ncore > 1,
      ncore = opts$ncore,
      verbose = opts$verbose
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    warning("HIMA did not return selected mediators: ", fit$message)
    results <- data.table()
  } else {
    results <- hima_to_table(fit, hima_inputs$mediator_map)
  }
  if (nrow(results) == 0) {
    results <- data.table(
      mediator = character(),
      hima_id = character(),
      alpha = numeric(),
      beta = numeric(),
      indirect_effect = numeric(),
      relative_importance_pct = numeric(),
      p_value = numeric()
    )
  }

  run_summary <- data.table(
    y_type = opts$y_type,
    mediator_type = opts$mediator_type,
    penalty = penalty,
    sigcut = opts$sigcut,
    n_complete_case = nrow(prepared$data),
    n_mediators_tested = ncol(hima_inputs$mediator_matrix),
    n_mediators_selected = nrow(results),
    n_covariate_terms = nrow(hima_inputs$covariate_map)
  )

  data.table::fwrite(
    results,
    file.path(opts$outdir, paste0(opts$output_prefix, "_selected_mediators.csv"))
  )
  data.table::fwrite(
    hima_inputs$mediator_map,
    file.path(opts$outdir, paste0(opts$output_prefix, "_mediator_map.csv"))
  )
  data.table::fwrite(
    hima_inputs$covariate_map,
    file.path(opts$outdir, paste0(opts$output_prefix, "_covariate_map.csv"))
  )
  data.table::fwrite(
    run_summary,
    file.path(opts$outdir, paste0(opts$output_prefix, "_run_summary.csv"))
  )

  message("HIMA mediation screening completed.")
  if (opts$verbose) {
    print(run_summary)
  }
}

main()
