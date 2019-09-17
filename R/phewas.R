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
    model_rhs = "character"
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
    model_rhs = ""
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
    model_rhs = model_rhs
  )
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


#' Main function to run a pheWAS.
#'
#' This function takes a configuration object and dispatches to the relevant
#' function.
#'
#' @param con A database connection to the UK Biobank (StatGen format).
#' @param configuration A configuration object.
#'
#' @export
runPheWAS <- function(con, configuration) {

  codes <- get_icd10_with_at_least_n_cases(
    con,
    configuration@min_num_cases,
    configuration@include_secondary_hospit,
    configuration@include_death_records
  )

  # We get all cases in memory and we do the filtering in R.
  all_cases <- get_full_records(
    con,
    configuration@include_secondary_hospit,
    configuration@include_death_records
  )

  # There are two ways of aggregating codes prior to analysis.
  # 1. By blocks (i.e. chunks of codes)
  # 2. By 3 characters codes (e.g. I20)

  # By default I will do both as well as no pre-processing and write the
  # results separately.
  ### blocks_results <- run_block_phewas(configuration, con, all_cases, codes)

  ### list(blocks_results=blocks_results)
  three_char_icd10_results <- run_3_char_phewas(configuration, con, all_cases,
                                                codes)

  list(three_char_icd10_results = three_char_icd10_results)
}


#' Helper function to run a pheWAS on blocks of ICD10 codes.
#'
#' @param configuration A pheWAS configuration object.
#' @param con A connection to the DB.
#' @param all_cases A dataframe with all of the EMR or death record data.
#' @param codes All ICD10 codes that were found in the DB with enough cases.
#'
#' @import here
#' @import parallel
#' @import doParallel
#' @import foreach
run_block_phewas <- function(configuration, con, all_cases, codes) {
  # Read the block metadata.
  blocks <- read.csv(here("../../data/icd10/icd10_blocks.csv"))

  cl <- makeCluster(configuration@ncpus)
  registerDoParallel(cl)  # nolint

  results <- foreach(
    i = 1:nrow(blocks),
    .combine = "rbind",
    .packages = c("UKBPheWAS", "broom")
  ) %dopar% {

    # Disable linting because 'i' is not found by the linter even though it
    # gets defined by dopar.
    block <- blocks[i, ]  # nolint

    # Check to see if some included diagnostic codes are in the current
    # iteration's block.
    matching_codes <- codes[
      sapply(codes, function(code) {
        UKBPheWAS::code_in_range(code, block$left, block$right)
      })
    ]

    cur_result <- NULL

    if (length(matching_codes) == 0) {
      print(paste0("Ignoring ", block$block, " because of insufficient cases"))
    }

    else {
      # Define cases.
      cases <- all_cases[all_cases$diag_icd10 %in% matching_codes, ]["eid"]

      # Do the regression.
      cur_result <- internal_do_logistic(cases, configuration@xs,
                                         configuration@model_rhs)

      # Filter the results to keep only relevant rows.
      cur_result <- apply_predicate_to_results(cur_result,
                                               configuration@voi_filter)

      # Add information on the current block.
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
#' @param codes All ICD10 codes that were found in the DB with enough cases.
#'
#' @import parallel
#' @import doParallel
#' @import foreach
run_3_char_phewas <- function(configuration, con, all_cases, codes) {
  three_char_codes <- unique(
    sapply(codes, function(code) {
      substr(code, 1, 3)
    })
  )

  cl <- makeCluster(configuration@ncpus)
  registerDoParallel(cl)  # nolint

  results <- foreach(
    code = three_char_codes,
    .combine = "rbind",
    .packages = c("broom")
  ) %dopar% {

    # Find all cases that match the 3 character ICD10 code.
    all_cases$short_code <- substr(all_cases$diag_icd10, 1, 3)
    cases <- all_cases[all_cases$short_code == code, ]["eid"]

    # Make sure we got a data frame with a single column named "eid".
    if (names(cases) != "eid")
      stop("Error identifying cases for 3 character ICD10 code.")

    cur_result <- internal_do_logistic(
      cases,
      configuration@xs,
      configuration@model_rhs
    )

    # Keep only the variable of interest.
    cur_result <- apply_predicate_to_results(cur_result,
                                             configuration@voi_filter)

    # Remember the 3 character code used.
    cur_result <- cbind(
      data.frame(outcome_description = code), cur_result
    )

    cur_result

  }

  stopCluster(cl)

  results
}


#' Utility function for row-wise filtering of results data frame.
apply_predicate_to_results <- function(results, predicate) {
  results[
    sapply(1:nrow(results), function(i) {
      predicate(i, results[i, ])
    }),
  ]
}


#' Do the logistic regression.
#'
#' Do the actual regression from the current cases and the prepared covariates
#' dataframe.
#'
#' @param cur_cases A column data frame of sample IDs to consider to be cases.
#' @param xs_df A matrix of covariates to be included in the multivariate
#'              regression.
internal_do_logistic <- function(cur_cases, xs_df, model_rhs) {

  cur_cases$case <- 1

  # Join the cases with the covariables.
  # We assume the first column in the XS matrix is sample IDs.
  # This is a right outer join where we assume that samples not in the
  # cur_cases df are controls.
  df <- merge(
    cur_cases, xs_df,
    by.x = "eid", by.y = names(xs_df)[1],
    all.x = FALSE, all.y = TRUE
  )

  df[is.na(df$case), "case"] <- 0

  # Drop the missing values.
  df <- df[complete.cases(df), ]

  # Prepare the model.
  if (is.null(model_rhs)) {
    # Treat all columns as covariables except for the sample_id column which is
    # assumed to be the first.
    cols <- names(xs_df)
    formula <- paste(
      "case ~ ",
      paste0(cols[2:length(cols)], collapse = " + ")
    )
  }
  else {
    formula <- paste0("case ~ ", model_rhs)
  }

  formula <- as.formula(formula)

  fit <- glm(formula, data = df, family = "binomial")

  format_fit(fit, df)
}


#' Format a GLM fit object.
#'
#' This function calculates some extra statistics and formats the fit object
#' as a dataframe.
#'
#' @import broom
format_fit <- function(fit, data) {
  df <- broom::tidy(fit)
  df$nobs <- nobs(fit)
  df$n_cases <- sum(data$case == 1)
  df$n_controls <- sum(data$case == 0)
  df$prevalence <- df$n_cases / df$nobs

  # Add the 95% CI
  ci <- confint.default(fit)
  colnames(ci) <- c("ci_95_low", "ci_95_high")
  df <- cbind(df, ci)

  df
}
