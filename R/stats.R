#' Do the logistic regression.
#'
#' Do the actual regression from the current cases and the prepared covariates
#' dataframe.
#'
#' This function is a callback for data generators.
#'
#' @seealso fastglm
do_logistic <- function(configuration, data, exclude = NULL) {
  binary_conf <- configuration$binary_configuration

  # Check if there are enough cases to test.
  n_cases <- nrow(data$y)
  if (n_cases < binary_conf$min_num_cases || n_cases == 0) {
    return(NULL)
  }

  # Join the cases with the covariables.
  # This also sets case status. Individuals are assumed control if not cases.
  df <- join_y_x_from_cases(data$y, configuration$xs)

  # Set excluded individuals if needed. This is used to cleanup controls
  # (remove individuals with evidence of disease) for cancer codes.
  if (!is.null(exclude)) {
    df[df$sample_id %in% exclude, "case"] <- NA
  }

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

  formatted <- format_logistic_fit(fit, df, binary_conf$use_fastglm)

  # Add information on the outcome.
  formatted <- cbind(
    data.frame(outcome_id = data$id, outcome_label = data$label),
    formatted
  )

  formatted

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


prepare_formula <- function(outcome_label, configuration, as.char = FALSE) {
  # If no RHS specified take everything in xs (except sample_id)
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

  if (as.char) {
    return(formula)
  }
  else {
    return(as.formula(formula))
  }

}


#' Do the linear regression.
do_linear <- function(configuration, data) {
  y <- data$y
  df <- inner_join_y_x(y, configuration$xs)

  if (nrow(df) == 0) {
    return(NULL)
  }

  formula <- prepare_formula(names(y)[2], configuration)

  fit <- lm(formula, data = df)

  results <- broom::tidy(fit)

  # Add the CI and information on the outcome
  ci <- as.data.frame(confint.default(fit))
  names(ci) <- c("ci_95_low", "ci_95_high")

  results <- cbind(
    data.frame(outcome_id = data$id, outcome_label = data$label),
    results,
    ci
  )

  rownames(results) <- NULL

  results

}


do_lm_F_test <- function(configuration, data) {
  df <- inner_join_y_x(data$y, configuration$xs)
  if (nrow(df) == 0) { return(NULL) }
  f <- gof(df, names(data$y)[2], configuration, "linear")

  return(data.frame(
    outcome_id = data$id,
    outcome_label = data$label,
    rss_base = f[1, 2],
    rss_augmented = f[2, 2],
    sum_of_sq = f[2, 4],
    F_stat = f[2, 5],
    p = f[2, 6]
  ))

}


do_logistic_LRT_test <- function(configuration, data, exclude = NULL) {
  y <- data$y
  df <- join_y_x_from_cases(y, configuration$xs)

  # Set excluded individuals if needed.
  if (!is.null(exclude)) {
    df[df$sample_id %in% exclude, "case"] <- NA
  }

  n_cases <- sum(df[, "case"] == 1)
  n_controls <- sum(df[, "case"] == 0)
  n_excl_from_ctrls <- sum(is.na(df[, "case"]))

  if (nrow(df) == 0) { return(NULL) }
  lrt <- gof(df, "case", configuration, "binary")

  return(data.frame(
    outcome_id = data$id,
    outcome_label = data$label,
    n_cases = n_cases,
    n_controls = n_controls,
    n_excl_from_ctrls = n_excl_from_ctrls,
    prevalence = n_cases / (n_cases + n_controls),
    resid_deviance_base = lrt[1, 2],
    resid_deviance_augmented = lrt[2, 2],
    deviance = lrt[2, 4],
    p = lrt[2, 5]
  ))
}


gof <- function(df, outcome, configuration, mode) {

  if (mode == "linear") {
    aug <- configuration$linear_configuration$augmented_variables
    fitter <- lm
    stat_test <- anova
  }

  else if (mode == "binary") {
    aug <- configuration$binary_configuration$augmented_variables
    fitter <- function(...) {
      glm(..., family = "binomial")
    }

    stat_test <- function(...) {
      anova(..., test = "LRT")
    }
  }
  else {
    stop("Unexpected mode for GOF.")
  }

  base_model_formula <- prepare_formula(
    outcome, configuration, as.char = TRUE
  )

  augmented_model_formula <- paste0(
    base_model_formula, " + ",
    paste(aug, collapse = " + ")
  )

  base_model_formula <- as.formula(base_model_formula)
  augmented_model_formula <- as.formula(augmented_model_formula)

  fit_base <- fitter(base_model_formula, data = df)
  fit_aug <- fitter(augmented_model_formula, data = df)

  stat_test(fit_base, fit_aug)

}


join_y_x <- function(y, x, all.x, all.y) {
  if (ncol(y) != 2) {
    stop(paste0(
      "Expected y to be a matrix of n x 2 where the first column is the ",
      "sample_id."
    ))
  }

  df <- merge(
    y, x,
    by.x = names(y)[1],
    by.y = names(x)[1],
    all.x = all.x, all.y = all.y
  )

  if (nrow(df) == 0) {
    warning("Join with covariables had no overlap.")
  }

  return(df)

}


inner_join_y_x <- function(y, x) {
  df <- join_y_x(y, x, all.x = TRUE, all.y = TRUE)
  df[complete.cases(df), ]
}


join_y_x_from_cases <- function(cases, x) {
  cases$case <- 1
  df <- join_y_x(cases, x, all.x = FALSE, all.y = TRUE)
  df[is.na(df$case), "case"] <- 0

  return(df)
}
