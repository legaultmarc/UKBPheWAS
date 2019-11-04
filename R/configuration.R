#' @export DoBinary
DoBinary <- function(
  include_death_records = TRUE,
  include_secondary_hospit = TRUE,
  min_num_cases = 50,
  use_fastglm = TRUE
) {
  binary_config <- list(
    include_secondary_hospit = include_secondary_hospit,
    include_death_records = include_death_records,
    min_num_cases = min_num_cases,
    use_fastglm = use_fastglm,
    callback = do_logistic
  )

  class(binary_config) <- "BinaryConfig"

  binary_config
}


#' @export DoBinaryLRT
DoBinaryLRT <- function(
  augmented_variables,

  include_death_records = TRUE,
  include_secondary_hospit = TRUE,
  min_num_cases = 50
) {
  binary_config <- list(
    include_death_records = include_death_records,
    include_secondary_hospit = include_secondary_hospit,
    min_num_cases = min_num_cases,

    augmented_variables = augmented_variables,

    callback = do_logistic_LRT_test
  )

  class(binary_config) <- "BinaryConfig"

  binary_config
}


#' @export SkipBinary
SkipBinary <- function() {
  binary_config <- list()
  class(binary_config) <- c("BinaryConfig", "Skip")

  binary_config
}


#' @export DoLinear
DoLinear <- function() {
  linear_config <- list(callback = do_linear)
  class(linear_config) <- "LinearConfig"

  linear_config
}


#' @export DoFTest
DoFTest <- function(augmented_variables) {
  linear_config <- list(
    augmented_variables = augmented_variables,
    callback = do_lm_F_test
  )
  class(linear_config) <- "LinearConfig"

  linear_config
}

#' @export SkipLinear
SkipLinear <- function() {
  linear_config <- list()
  class(linear_config) <- c("LinearConfig", "Skip")

  linear_config
}


#' Helper function to create Configuration objects.
#'
#' @export phewas_configuration
phewas_configuration <- function(
  db_password = NULL,
  limit = NULL,
  ncpus = 10,
  xs = NULL,
  model_rhs = "",
  output_prefix = "phewas",
  continuous_variables_path = "/data/projects/uk_biobank/data/reports/SGR-2094/",

  voi_idx = NULL,
  voi_name = NULL,
  results_filter = NULL
) {

  # Check that some sort of variable of interest filtering has been provided.
  n_voi_filters <- sum(
    !sapply(list(voi_idx, voi_name, results_filter), is.null)
  )

  if (n_voi_filters == 0) {
    warning(paste0(
      "No variable of interest filtering has been provided. All parameter ",
      "results will be saved which will require a lot of memory. ",
      "You can provide 'voi_idx', 'voi_name' or 'results_filter' to control ",
      "which parameters are reported."
    ))

    voi_filter <- function(i, row) TRUE
  }

  else if (n_voi_filters > 1) {
    stop(paste0(
      "More than one variable of interest filters have been provided. ",
      "For complex filtering, only define the 'results_filter'."
    ))
  }

  else {
    # There is only one kind of filtering that was provided.
    if (!is.null(voi_idx))
      results_filter <- variable_of_interest_index(voi_idx)

    else if (!is.null(voi_name))
      results_filter <- variable_of_interest_name(voi_name)

    else if (!is.null(results_filter))
      results_filter <- results_filter

    else
      stop("Unexpected error (predicates).")
  }

  # Check that some exogenous variables have been provided.
  if (nrow(xs) < 2) {
    stop(paste0("Invalid regressors: '", xs, "'. The xs matrix needs ",
                "at least two rows."))
  }

  if (nchar(model_rhs) == 0) {
    stop(paste0("Expected a right hand side for the formula specifying ",
                "exogenous variables. Got: '", model_rhs, "'"))
  }

  # Check that a database password has been provided.
  if (is.null(db_password)) {
    warning(
      "No database password provided to connect to the UK Biobank database"
    )
  }

  config <- list(
    db_password = db_password,
    limit = limit,
    ncpus = ncpus,
    xs = xs,
    model_rhs = model_rhs,
    output_prefix = output_prefix,
    continuous_variables_path = continuous_variables_path,
    binary_configuration = DoBinary(),
    linear_configuration = DoLinear(),
    results_filter = results_filter
  )

  class(config) <- "Configuration"

  config
}


should_do_linear <- function(configuration) {
  !is(configuration$linear_configuration, "Skip")
}


should_do_binary <- function(configuration) {
  !is(configuration$binary_configuration, "Skip")
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
