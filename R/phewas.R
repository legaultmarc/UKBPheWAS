library(here)
library(doParallel)  # nolint
library(foreach)

# Get covariates and samples to analyse
# Get list of codes to test
# Define the minimum number of cases
# Define the code levels
#   By default only observed codes but allow parent codes to be tested at the
#     - block level
#     - level 1 (e.g. I23.8 -> I23)
# Option for primary only or primary + secondary
# Option to include death records
#
# Note always match children codes.


#' A configuration class to parametrize PheWAS analyses.
#'
#' @export Configuration
#' @exportClass Configuration
Configuration <- setClass(
  "Configuration",
  slots = list(
    include_secondary_hospit = "logical",
    include_death_records = "logical",
    min_num_cases = "numeric",
    ncpus = "numeric",
    voi_filter = "function",
    xs = "data.frame",
    model_rhs = "character",
    output_prefix = "character",
    use_fastglm = "logical"
  )
)


#' Construction function to create configuration objects for pheWAS analyses.
#'
#' This does mostly parameter checking.
#'
#' @export
create_configuration <- function(
    include_secondary_hospit = TRUE,
    include_death_records = TRUE,
    min_num_cases = 50,
    ncpus = 3,
    voi_idx = NULL,
    voi_name = NULL,
    voi_predicate = NULL,
    xs = NULL,
    model_rhs = "",
    output_prefix = "phewas",
    use_fastglm = FALSE
  ) {

  # Check that some sort of variable of interest filtering has been provided.
  n_voi_filters <- sum(
    !sapply(list(voi_idx, voi_name, voi_predicate), is.null)
  )

  if (n_voi_filters == 0) {
    warning(paste0(
      "No variable of interest filtering has been provided. All parameter ",
      "results will be saved which will require a lot of memory. ",
      "You can provide 'voi_idx', 'voi_name' or 'voi_predicate' to control ",
      "which parameters are reported."
    ))

    voi_filter <- function(i, row) TRUE
  }

  else if (n_voi_filters > 1) {
    stop(paste0(
      "More than one variable of interest filters have been provided. ",
      "For complex filtering, only define the 'voi_predicate'."
    ))
  }

  else {
    # There is only one kind of filtering that was provided.
    if (!is.null(voi_idx))
      voi_filter <- variable_of_interest_index(voi_idx)

    else if (!is.null(voi_name))
      voi_filter <- variable_of_interest_name(voi_name)

    else if (!is.null(voi_predicate))
      voi_filter <- voi_predicate

    else
      stop("Unexpected error (predicates).")
  }

  Configuration(
    include_secondary_hospit = include_secondary_hospit,
    include_death_records = include_death_records,
    min_num_cases = min_num_cases,
    ncpus = ncpus,
    voi_filter = voi_filter,
    xs = xs,
    model_rhs = model_rhs,
    output_prefix = output_prefix,
    use_fastglm = use_fastglm
  )
}


#' Logging function that includes a timestamp.
log_message <- function(message) {
  cat(paste0(format(Sys.time(), "%Y-%m-%d at %X"), " -- ",  message, "\n"))
}


#' Predicate generator function for variable of interest filtering based on
#' parameter names.
variable_of_interest_name <- function(name) {
  predicate <- function(i, row) {
    row$term == name
  }
  predicate
}


#' Predicate generator function for variable of interest filtering based on
#' parameter indices.
variable_of_interest_index <- function(idx=2) {
  predicate <- function(i, row) {
    i == idx
  }
  predicate
}


#' Logs the configuration and welcome message.
log_config <- function(configuration, biomarker = FALSE) {
  log_message("Starting UKBPheWAS\n")
  cat("Running analysis with the following parameters:\n")

  if (!biomarker) {
    cat("    - Logistic regression of HES and death records data\n\n")
    cat(paste0(
      "    - Include secondary hospitalization codes: ",
      configuration@include_secondary_hospit, "\n"
    ))
    cat(paste0(
      "    - Include death records: ",
      configuration@include_death_records, "\n"
    ))
    cat(paste0(
      "    - Minimum number of cases for inclusion: ",
      configuration@min_num_cases, "\n"
    ))
    cat(paste0(
      "    - Use of fastglm package: ",
      configuration@use_fastglm, "\n\n"
    ))

  }
  else {
    cat("    - Linear regression of biomarkers\n\n")
  }

  cat(paste0(
    "    - Number of CPUs to use: ",
    configuration@ncpus, "\n"
  ))

}


#' Main function to run a pheWAS on biomarker data.
#'
#' @param con A database connection to the Biomarker UK Biobank database.
#' @param configuration A configuration object.
#'
#' @export
runBiomarkerPheWAS <- function(con, configuration) {

  log_config(configuration, biomarker = TRUE)

  log_message("Running analysis on biomarkers...")
  results <- internal_run_biomarker_phewas(configuration, con)
  write.csv(
    results,
    paste0(configuration@output_prefix, "_biomarkers.csv")
  )
  cat("Done!\n\n")

  results
}


#' Main function to run a pheWAS on hospit and/or death data.
#'
#' This function takes a configuration object and dispatches to the relevant
#' function.
#'
#' @param con A database connection to the UK Biobank (StatGen format).
#' @param configuration A configuration object.
#'
#' @export
runPheWAS <- function(con, configuration) {

  log_config(configuration)

  # We get all cases in memory and we do the filtering in R.
  all_cases <- get_full_records(
    con,
    configuration@include_secondary_hospit,
    configuration@include_death_records
  )

  # We exclude codes starting with S-Z as they correspond to special codings or
  # external causes that are likely irrelevant in the context of a pheWAS.
  excluded_chapters <- c("S", "T", "U", "V", "W", "X", "Y", "Z")

  all_cases <- all_cases[
    !(substr(all_cases$diag_icd10, 1, 1) %in% excluded_chapters),
  ]

  # There are two ways of aggregating codes prior to analysis.
  # 1. By blocks (i.e. chunks of codes)
  # 2. By 3 characters codes (e.g. I20)

  # By default I will do both as well as no pre-processing and write the
  # results separately.
  log_message("Running analysis for blocks of ICD10 codes...")
  blocks_results <- run_block_phewas(configuration, con, all_cases)
  write.csv(blocks_results, paste0(configuration@output_prefix, "_blocks.csv"))
  cat("Done!\n\n")

  log_message("Running analysis for 3 character ICD10 codes...")
  three_char_icd10_results <- run_3_char_phewas(configuration, con, all_cases)
  write.csv(
    three_char_icd10_results,
    paste0(configuration@output_prefix, "_three_chars.csv")
  )
  cat("Done!\n\n")

  log_message("Running analysis for naive ICD10 codes (no pre-processing)...")
  naive_results <- run_naive_phewas(configuration, con, all_cases)
  write.csv(
    naive_results,
    paste0(configuration@output_prefix, "_naive.csv")
  )
  cat("Done!\n\n")

  log_message("All analyses completed")

  list(
    blocks_results = blocks_results,
    three_char_icd10_results = three_char_icd10_results,
    naive_results = naive_results
  )

}


#' Helper function to run a pheWAS on blocks of ICD10 codes.
#'
#' @param configuration A pheWAS configuration object.
#' @param con A connection to the DB.
#' @param all_cases A dataframe with all of the EMR or death record data.
#'
#' @import here
#' @import parallel
#' @import doParallel
#' @import foreach
run_block_phewas <- function(configuration, con, all_cases) {
  # Read the block metadata.
  blocks <- read.csv(here("../../data/icd10/icd10_blocks.csv"))

  cl <- makeCluster(configuration@ncpus)
  registerDoParallel(cl)  # nolint

  results <- foreach(
    i = 1:nrow(blocks),
    .combine = "rbind",
    .packages = c("UKBPheWAS", "broom")
  ) %dopar% {

    block <- blocks[i, ]

    # Define cases.
    cases <- all_cases[
      sapply(all_cases$diag_icd10, function(code) {
        UKBPheWAS::code_in_range(code, block$left, block$right)
      }),
    ]["eid"]

    # Do the regression.
    cur_result <- internal_do_logistic(cases, configuration)

    # Add information on the current block.
    if (!is.null(cur_result)) {
      cur_result <- cbind(
        data.frame(outcome_description = block$block), cur_result
      )
    }

    cur_result

  }

  stopCluster(cl)

  results
}


#' Run a 3 character ICD10 code pheWAS.
#'
#' @param configuration A pheWAS configuration object.
#' @param con A connection to the DB.
#' @param all_cases A dataframe with all of the EMR or death record data.
#'
#' @import parallel
#' @import doParallel
#' @import foreach
run_3_char_phewas <- function(configuration, con, all_cases) {

  all_cases$short_code <- substr(all_cases$diag_icd10, 1, 3)
  three_char_codes <- unique(all_cases$short_code)

  cl <- makeCluster(configuration@ncpus)
  registerDoParallel(cl)  # nolint

  results <- foreach(
    code = three_char_codes,
    .combine = "rbind",
    .packages = c("broom")
  ) %dopar% {

    # Find all cases that match the 3 character ICD10 code.
    cases <- all_cases[all_cases$short_code == code, ]["eid"]

    # Make sure we got a data frame with a single column named "eid".
    if (names(cases) != "eid")
      stop("Error identifying cases for 3 character ICD10 code.")

    cur_result <- internal_do_logistic(cases, configuration)

    # Remember the 3 character code used.
    if (!is.null(cur_result)) {
      cur_result <- cbind(data.frame(outcome_description = code), cur_result)
    }

    cur_result

  }

  stopCluster(cl)

  results
}


#' Run a naive pheWAS
#'
#' @param configuration A pheWAS configuration object.
#' @param con A connection to the DB.
#' @param all_cases A dataframe with all of the EMR or death record data.
#'
#' @import parallel
#' @import doParallel
#' @import foreach
run_naive_phewas <- function(configuration, con, all_cases) {

  cl <- makeCluster(configuration@ncpus)
  registerDoParallel(cl)  # nolint

  results <- foreach(
    code = unique(all_cases$diag_icd10),
    .combine = "rbind",
    .packages = c("broom")
  ) %dopar% {

    # We compare up to the length of the current code.
    cur_code_length <- nchar(code)

    cases <- all_cases[
      substr(all_cases$diag_icd10, 1, cur_code_length) == code,
    ]["eid"]

    # Make sure we got a data frame with a single column named "eid".
    if (names(cases) != "eid")
      stop("Error identifying cases for the naive pheWAS.")

    cur_result <- internal_do_logistic(cases, configuration)

    # Remember the code used.
    if (!is.null(cur_result)) {
      cur_result <- cbind(data.frame(outcome_description = code), cur_result)
    }

    cur_result

  }

  stopCluster(cl)

  results
}


#' Run a pheWAS on biomarker traits.
#'
#' @param configuration A pheWAS configuration object.
#' @param con A connection to the biomarker DB.
#'
#' @import parallel
#' @import doParallel
#' @import foreach
internal_run_biomarker_phewas <- function(configuration, con) {

  # Get biomarkers.
  data <- get_all_biomarkers(con)

  variable_ids <- unique(data$variable_id)

  cl <- makeCluster(configuration@ncpus)
  registerDoParallel(cl)  # nolint

  results <- foreach(
    var = variable_ids,
    .combine = "rbind",
    .packages = c("broom")
  ) %dopar% {

    # Filter the data.
    cur <- data[data$variable_id == var, ]

    cur_result <- internal_do_linear(cur[, c("sample_id", "val")],
                                     configuration)

    if (!is.null(cur_result)) {
      cur_result <- cbind(
        # Add info on the outcome to the results.
        data.frame(
          outcome_description = cur[1, "description"],
          ukb_variable_id = var
        ),
        cur_result
      )
    }
    else {
      warning(paste0("Linear regression failed for variable ", var))
    }

    cur_result

  }

  stopCluster(cl)

  results
}


#' Utility function for row-wise filtering of results data frame.
apply_predicate_to_results <- function(results, predicate) {
  results <- results[
    sapply(1:nrow(results), function(i) {
      predicate(i, results[i, ])
    }),
  ]

  if (nrow(results) == 0) {
    warning("No results selected by variable of interest filter.")
  }

  results
}


#' Do the logistic regression.
#'
#' Do the actual regression from the current cases and the prepared covariates
#' dataframe.
#'
#' @param cur_cases A column data frame of sample IDs to consider to be cases.
#' @param configuration A configuration object containing information on the
#'                      analysis as well as the covariables matrix.
#'
#' @seealso fastglm
internal_do_logistic <- function(cur_cases, configuration) {

  # Check if there are enough cases to test.
  if (nrow(cur_cases) < configuration@min_num_cases || nrow(cur_cases) == 0) {
    return(NULL)
  }

  cur_cases$case <- 1

  # Join the cases with the covariables.
  # We assume the first column in the XS matrix is sample IDs.
  # This is a right outer join where we assume that samples not in the
  # cur_cases df are controls.
  df <- merge(
    cur_cases, configuration@xs,
    by.x = "eid", by.y = names(configuration@xs)[1],
    all.x = FALSE, all.y = TRUE
  )

  # Assume individuals not in the "cases" DF are controls.
  df[is.na(df$case), "case"] <- 0

  # Drop the missing values.
  df <- df[complete.cases(df), ]

  # Prepare the model.
  formula <- prepare_formula("case", configuration)

  if (configuration@use_fastglm) {
    # Use the LLT method as it is one of the fastest.
    # Other choices are: 1:qr; 3:LDLT
    # See doc: https://github.com/jaredhuling/fastglm
    # For installation and details.
    fit <- fastglm::fastglm(
      x = model.matrix(formula, data = df),
      y = df$case,
      data = df,
      family = binomial(),
      method = 2
    )
  }
  else {
    fit <- glm(formula, data = df, family = "binomial")
  }

  results <- format_logistic_fit(fit, df, configuration@use_fastglm)

  # Apply the predicate to filter the results.
  apply_predicate_to_results(results, configuration@voi_filter)

}

#' A "copy" of broom::tidy for fastglm objects.
#'
#' The goal is to achieve results close enough or identical to 
#' format_logistic_fit so that the use of fastglm is mostly transparent.
tidy_fastglm <- function(fit) {
  # This reproduces the results from broom::tidy
  df <- data.frame(
    term = names(fit$coefficients),
    estimate = fit$coefficients,
    std.error = fit$se
  )
  df$statistic <- df$estimate / df$std.error
  df$p.value <- 2 * pnorm(-abs(df$statistic))

  df
}


#' Format a GLM fit object.
#'
#' This function calculates some extra statistics and formats the fit object
#' as a dataframe.
#'
#' @import broom
format_logistic_fit <- function(fit, data, fastglm_format) {
  if (fastglm_format) {

    df <- tidy_fastglm(fit)

    # Manually calculate the Wald CI.
    ci <- data.frame(
      ci_95_low = df$estimate + qnorm(0.05 / 2) * df$std.error,
      ci_95_high = df$estimate - qnorm(0.05 / 2) * df$std.error
    )

    # We use complete cases before so this should be safe to estimate the
    # number of observations.
    df$nobs <- nrow(data)
  }

  else {
    df <- broom::tidy(fit)
    ci <- as.data.frame(confint.default(fit))
    names(ci) <- c("ci_95_low", "ci_95_high")
    df$nobs <- nobs(fit)
  }

  df$n_cases <- sum(data$case == 1)
  df$n_controls <- sum(data$case == 0)
  df$prevalence <- df$n_cases / df$nobs

  # Add the 95% CI
  df <- cbind(df, ci)

  df
}


prepare_formula <- function(outcome_label, configuration) {
  if (configuration@model_rhs == "") {
    cols <- names(configuration@xs)

    formula <- paste0(
      outcome_label, " ~ ",
      paste0(cols[2:length(cols)], collapse = " + ")
    )
  }

  else {
    formula <- paste0(outcome_label, " ~ ", configuration@model_rhs)
  }

  as.formula(formula)

}


#' Do the linear regression.
#'
#' Do the linear regression for the current biomarker.
#'
#' @param y A data.frame with two columns: sample_id and value (assumed to be
#'          provided in that order).
#' @param configuration A configuration object containing information on the
#'                      analysis as well as the covariables matrix.
#'
#' @import broom
internal_do_linear <- function(y, configuration) {
  df <- merge(
    y, configuration@xs,
    by.x = names(y)[1],
    by.y = names(configuration@xs)[1],
    all.x = FALSE, all.y = TRUE
  )

  df <- df[complete.cases(df), ]

  if (nrow(df) == 0) {
    return(NULL)
  }

  formula <- prepare_formula("val", configuration)

  fit <- lm(formula, data = df)

  results <- broom::tidy(fit)

  # Add the CI.
  ci <- as.data.frame(confint.default(fit))
  names(ci) <- c("ci_95_low", "ci_95_high")
  results <- cbind(results, ci)

  apply_predicate_to_results(results, configuration@voi_filter)

}
