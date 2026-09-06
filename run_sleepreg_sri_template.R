#!/usr/bin/env Rscript

# Template runner for calculating sleep regularity index from GGIR outputs.
# Required GGIR layout: each selected output directory should contain
# results/QC/part4_nightsummary_sleep_full.csv.

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    ggir_parent = NA_character_,
    outdir = NA_character_,
    min_sri_days = 5,
    sidef = "T5A5",
    use_naps = TRUE,
    use_waso = TRUE,
    use_miscal = TRUE,
    use_ggir_nonwear = TRUE,
    use_custom_nonwear = FALSE,
    nonwear_in_ggir_sleep = FALSE,
    write_swv = FALSE,
    write_raster = FALSE,
    ggir_nonwear_stage = "final",
    install_sleepreg = FALSE,
    output_pattern = "^output_",
    exclude_pattern = "",
    list_only = FALSE,
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
    key <- kv[1]
    val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "TRUE"
    key <- gsub("-", "_", key)

    if (!key %in% names(opts)) stop("Unknown option: --", key)
    if (is.logical(opts[[key]])) {
      opts[[key]] <- bool_value(val)
    } else if (is.numeric(opts[[key]])) {
      opts[[key]] <- as.numeric(val)
    } else {
      opts[[key]] <- val
    }
  }

  if (is.na(opts$ggir_parent) || !nzchar(opts$ggir_parent)) {
    stop("Set --ggir_parent to the parent directory containing GGIR output folders.")
  }
  if (is.na(opts$outdir) || !nzchar(opts$outdir)) {
    stop("Set --outdir to the directory where SRI outputs should be written.")
  }

  opts$ggir_parent <- normalizePath(opts$ggir_parent, winslash = "/", mustWork = TRUE)
  opts$outdir <- normalizePath(opts$outdir, winslash = "/", mustWork = FALSE)
  opts
}

ensure_sleepreg <- function(install_sleepreg = FALSE) {
  if (!requireNamespace("sleepreg", quietly = TRUE)) {
    if (!install_sleepreg) {
      stop(
        "Package 'sleepreg' is not installed. Re-run with --install_sleepreg=TRUE ",
        "or install it manually before running this script."
      )
    }
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    remotes::install_github(
      "vincentvanhees/sleepreg",
      ref = "upgradeGGIRversion",
      upgrade = "never",
      dependencies = TRUE
    )
  }
  suppressPackageStartupMessages(library(sleepreg))
}

find_ggir_output_dirs <- function(ggir_parent, output_pattern = "^output_", exclude_pattern = "") {
  dirs <- list.dirs(ggir_parent, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[grepl(output_pattern, basename(dirs))]
  if (!is.na(exclude_pattern) && nzchar(exclude_pattern)) {
    dirs <- dirs[!grepl(exclude_pattern, basename(dirs))]
  }

  required <- file.path(dirs, "results", "QC", "part4_nightsummary_sleep_full.csv")
  dirs[file.exists(required)]
}

make_fastconfig_sri_function <- function(selected_output_dirs) {
  src <- deparse(sleepreg:::SRI_from_GGIR)
  src[1] <- paste0("SRI_from_GGIR_fastconfig <- ", src[1])

  fast_config_reader <- paste(
    "conf <- tryCatch({",
    "  config_lines <- readLines(confloc, warn = FALSE)",
    "  desired_line <- config_lines[grepl('^desiredtz,', config_lines)]",
    "  desired_tz <- ''",
    "  if (length(desired_line) > 0) {",
    "    desired_df <- utils::read.csv(text = paste(c('argument,value,context', desired_line[1]), collapse = '\\n'), stringsAsFactors = FALSE)",
    "    desired_tz <- desired_df$value[1]",
    "    if (is.na(desired_tz)) desired_tz <- ''",
    "  }",
    "  data.frame(argument = 'desiredtz', value = desired_tz, context = 'params_general', stringsAsFactors = FALSE)",
    "}, error = function(e) {",
    "  data.frame(argument = 'desiredtz', value = '', context = 'params_general', stringsAsFactors = FALSE)",
    "})",
    sep = "\n"
  )

  patched_config <- sub("conf <- read.csv\\(confloc\\)", fast_config_reader, src)
  if (identical(patched_config, src)) {
    stop("Could not patch sleepreg::SRI_from_GGIR config reader; package source may have changed.")
  }

  patched <- patched_config
  patched <- sub(
    "dir_list <- list.dirs\\(outputdir, recursive = FALSE\\)",
    "dir_list <- .selected_dir_list",
    patched
  )
  patched <- sub(
    "study_list <- list.dirs\\(outputdir, recursive = FALSE, full.names = FALSE\\)",
    "study_list <- .selected_study_list",
    patched
  )
  patched <- sub(
    "dir_list <- dir_list\\[grep\\(\"output_\", dir_list\\)\\]",
    "dir_list <- dir_list",
    patched
  )
  patched <- sub(
    "study_list <- study_list\\[grep\\(\"output_\", study_list\\)\\]",
    "study_list <- study_list",
    patched
  )

  env <- new.env(parent = .GlobalEnv)
  env$.selected_dir_list <- selected_output_dirs
  env$.selected_study_list <- basename(selected_output_dirs)
  eval(parse(text = patched), envir = env)
  env$SRI_from_GGIR_fastconfig
}

check_output <- function(sri_file, verbose = FALSE) {
  if (!file.exists(sri_file) || file.info(sri_file)$size == 0) {
    stop("SRI.csv was not created or is empty.")
  }

  if (verbose) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
      install.packages("data.table", repos = "https://cloud.r-project.org")
    }
    s <- data.table::fread(sri_file, nrows = 0)
    message("Output columns: ", paste(names(s), collapse = ", "))
  }

  invisible(TRUE)
}

main <- function() {
  opts <- parse_args()
  dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(opts$outdir, "run_sleepreg_sri.log")
  log_con <- file(log_file, open = "at")
  sink(log_con, split = TRUE)
  sink(log_con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)

  message("Starting sleep regularity index calculation.")
  message("minSRIdays: ", opts$min_sri_days)
  message("SIdef: ", opts$sidef)
  message("write SWV: ", opts$write_swv)
  message("write raster: ", opts$write_raster)
  message("output pattern: ", opts$output_pattern)

  output_dirs <- find_ggir_output_dirs(
    ggir_parent = opts$ggir_parent,
    output_pattern = opts$output_pattern,
    exclude_pattern = opts$exclude_pattern
  )
  if (length(output_dirs) == 0) {
    stop("No valid GGIR output directories found.")
  }
  message("Valid GGIR output directories found.")
  if (opts$list_only) {
    message("list_only=TRUE, stopping before SRI calculation.")
    return(invisible(output_dirs))
  }

  ensure_sleepreg(opts$install_sleepreg)
  if (opts$verbose) {
    message("sleepreg version: ", as.character(utils::packageVersion("sleepreg")))
    if (requireNamespace("GGIR", quietly = TRUE)) {
      message("GGIR version: ", as.character(utils::packageVersion("GGIR")))
    }
  }

  sri_fun <- make_fastconfig_sri_function(output_dirs)
  sri_fun(
    outputdir = opts$ggir_parent,
    alloutdir = opts$outdir,
    use.naps = opts$use_naps,
    use.WASO = opts$use_waso,
    use.miscal = opts$use_miscal,
    use.GGIRnonwear = opts$use_ggir_nonwear,
    use.customnonwear = opts$use_custom_nonwear,
    nonWearInGGIRsleep = opts$nonwear_in_ggir_sleep,
    wr.SWV = opts$write_swv,
    wr.raster = opts$write_raster,
    minSRIdays = opts$min_sri_days,
    SIdef = opts$sidef,
    GGIR_nonwear_detection_stage = opts$ggir_nonwear_stage
  )

  sri_file <- file.path(opts$outdir, "SRI.csv")
  check_output(sri_file, verbose = opts$verbose)
  message("Sleep regularity index calculation completed.")
  message("SRI output file created.")
}

main()
