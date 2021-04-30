library(fastglm)

args <- commandArgs(trailingOnly=TRUE)

worker_id <- args[1]
dealer_addr <- args[2]
monitor_addr <- args[3]


fit_fastglm <- function(data_matrix) {
  fastglm(
    # We add the intercept to the x matrix.
    x = cbind(1, data_matrix[, -1]),
    y = data_matrix[, 1],
    family = binomial(),
    method = 2
  )
}


fastglm_anova <- function(base, augmented) {
  d_dev <- base$deviance - augmented$deviance
  d_df <- base$df.residual - augmented$df.residual
  p <- 1 - pchisq(d_dev, d_df)

  list(
    resid_deviance_base = base$deviance,
    resid_deviance_augmented = augmented$deviance,
    deviance = d_dev,
    p = p
  )
}


logistic_deviance_diff_test_worker <- function(worker_id, ...) {
  # Read the analysis configuration object.
  conf <- fromJSON(file = "analysis_configuration.json")

  # Read the covariables.
  covars <- get_xs(conf)

  # Read the XPCs.
  xpcs <- read.csv(conf$binary_conf$xpcs_path)
  covars <- merge(covars, xpcs, by = "sample_id")

  # For now, sample_id needs to be a string. Eventually, we will infer this
  # better.
  covars$sample_id <- as.character(covars$sample_id)

  # Create a header file.
  cat(paste0("variable_id,analysis_type,sex_subset,",
             "n_cases,n_controls,",
             "n_excluded_from_controls,",
             "deviance_base,n_params_base,",
             "deviance_augmented,n_params_augmented,",
             "deviance_diff,p\n"),
      file="header_summary.csv")

  cat("R: Opened output file.\n")
  output_file <- file(
    paste0("results_worker_", worker_id, "_summary.csv"), open = "wt"
  )

  model_file <- file(
    paste0("results_worker_", worker_id, "_model.json"), open = "wt"
  )

  # Check if sex is in the independant variables.
  # This is important later when the variability in the sex variable is
  # tested after filtering.
  base_cols <- char_to_terms(conf$model_rhs)
  aug_cols <- conf$binary_conf$augmented_variables
  base_aug_cols <- c(base_cols, aug_cols)

  # Callback for when data is to be processed from the queue.
  do.work <- function(metadata, data) {
    data <- deserialize(data)

    cases <- as.character(data[data[, "y"] == 1, "sample_id"])
    to_exclude <- as.character(data[is.na(data[, "y"]), "sample_id"])

    min_num_cases <- conf$binary_conf$min_num_cases

    covars$y <- as.numeric(covars$sample_id %in% cases)
    covars[covars$sample_id %in% to_exclude, "y"] <- NA

    n_excluded <- sum(is.na(covars[, "y"]))

    aug_model_rhs <- paste(base_aug_cols, collapse = " + ")
    aug_formula <- as.formula(paste0(
      "y ~ ", aug_model_rhs
    ))

    aug_data_matrix <- model.frame(aug_formula, data = covars)
    keep <- complete.cases(aug_data_matrix)

    if (sum(keep) < min_num_cases) {
      return()
    }

    aug_data_matrix <- aug_data_matrix[keep, ]
    col_drop_li <- drop_columns_with_no_variance(aug_data_matrix)
    aug_data_matrix <- col_drop_li$mat
    dropped_cols <- col_drop_li$dropped_cols

    if (length(dropped_cols) != 0) {
      base_cols <- base_cols[!(base_cols %in% dropped_cols)]
      base_aug_cols <- base_aug_cols[!(base_aug_cols %in% dropped_cols)]
    }

    n_cases <- sum(aug_data_matrix[, "y"] == 1)
    n_controls <- sum(aug_data_matrix[, "y"] == 0)

    aug_data_matrix <- as.matrix(aug_data_matrix)

    fit_base <- fit_fastglm(aug_data_matrix[, c("y", base_cols)])
    fit_aug <- fit_fastglm(aug_data_matrix)

    lrt <- fastglm_anova(fit_base, fit_aug)

    line <- paste(
      metadata$variable_id,
      metadata$analysis_type,
      ifelse(is.null(metadata$sex_subset), "BOTH", metadata$sex_subset),
      n_cases,
      n_controls,
      n_excluded,
      lrt$resid_deviance_base,
      length(base_cols) + 1,  # Add 1 for intercept.
      lrt$resid_deviance_augmented,
      length(base_aug_cols) + 1,
      lrt$deviance,
      lrt$p,
      sep=","
    )

    writeLines(line, con = output_file)

    # We will write the model coefficients in case they are needed.
    infer_df <- data.frame(
      variable = names(fit_aug$coefficients),
      beta = fit_aug$coefficients,
      se = fit_aug$se,
      z = fit_aug$coefficients / fit_aug$se
    )
    infer_df$p <- 2 * pnorm(-abs(infer_df$z))
    infer_df$nlog10p <- -2 * pnorm(-abs(infer_df$z), log.p = T) / log(10)

    infer_df[1, "variable"] <- "intercept"

    cat(toJSON(
      list(
        variable_id = metadata$variable_id,
        analysis_type = metadata$analysis_type,
        model_fit = infer_df
      )
    ), file = model_file)
    cat("\n", file = model_file)

  }

  # Start the main loop.
  Worker(worker_id, ..., callback=do.work)

  cat("R: closing output file.\n")
  close(output_file)

}

# do.call(linear_f_test_worker, list(worker_id, dealer_addr, monitor_addr))
logistic_deviance_diff_test_worker(worker_id, dealer_addr, monitor_addr)
