#' Do the logistic regression.
#'
#' Do the actual regression from the current cases and the prepared covariates
#' dataframe.
#'
#' @seealso fastglm
do_logistic <- function(configuration, data) {
  binary_conf <- configuration$binary_configuration

  # Check if there are enough cases to test.
  n_cases <- nrow(data$y)
  if (n_cases < binary_conf$min_num_cases || n_cases == 0) {
    return(NULL)
  }

  data$y$case <- 1

  # Join the cases with the covariables.
  # We assume the first column in the XS matrix is sample IDs.
  # This is a right outer join where we assume that samples not in the
  # cur_cases df are controls.
  df <- merge(
    data$y, configuration$xs,
    by.x = "sample_id", by.y = names(configuration$xs)[1],
    all.x = FALSE, all.y = TRUE
  )

  # Assume individuals not in the "cases" DF are controls.
  df[is.na(df$case), "case"] <- 0

  # Drop the missing values.
  df <- df[complete.cases(df), ]

  # Prepare the model.
  formula <- prepare_formula("case", configuration)

  if (binary_conf$use_fastglm) {
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

  format_logistic_fit(fit, df, binary_conf$use_fastglm)

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

  # Strip rownames which are redundant with term.
  rownames(df) <- NULL

  df$n_cases <- sum(data$case == 1)
  df$n_controls <- sum(data$case == 0)
  df$prevalence <- df$n_cases / df$nobs

  # Add the 95% CI
  df <- cbind(df, ci)

  df
}


prepare_formula <- function(outcome_label, configuration) {
  if (configuration$model_rhs == "") {
    cols <- names(configuration$xs)

    formula <- paste0(
      outcome_label, " ~ ",
      paste0(cols[2:length(cols)], collapse = " + ")
    )
  }

  else {
    formula <- paste0(outcome_label, " ~ ", configuration$model_rhs)
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
    y, configuration$xs,
    by.x = names(y)[1],
    by.y = names(configuration$xs)[1],
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

  results

}
