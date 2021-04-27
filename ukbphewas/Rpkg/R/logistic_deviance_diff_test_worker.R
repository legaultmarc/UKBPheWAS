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

  # If there is a defined subset, we subset right away to accelerate analyses.
  if (!is.null(conf$subset)) {
    conf$subset <- as.character(conf$subset)
    cat(paste0(
      "R: Subsetting ", length(conf$subset), " individuals.\n"
    ))
    covars <- covars[covars$sample_id %in% conf$subset, ]
  }

  # Create a header file.
  cat(paste0("variable_id,analysis_type,n_cases,n_controls,",
             "n_excluded_from_controls,deviance_base,",
             "deviance_augmented,deviance_diff,p\n"),
      file="header.csv")

  cat("R: Opened output file.\n")
  output_file <- file(
    paste0("results_worker_", worker_id, ".csv"),
    open = "wt"
  )

  # Check if sex is in the independant variables.
  # This is important later when the variability in the sex variable is
  # tested after filtering.
  base_cols <- char_to_terms(conf$model_rhs)
  aug_cols = conf$binary_conf$augmented_variables

  sex_in_indep_variables <- (
    !is.null(conf$sex_column) && (conf$sex_column %in% c(base_cols, aug_cols))
  )

  # If the analysis is sex stratified, remove the sex from the independant
  # variables.
  if (conf$sex_stratified && sex_in_indep_variables) {
    base_cols  <- base_cols[base_cols != conf$sex_column]
    aug_cols  <- aug_cols[aug_cols != conf$sex_column]

    sex_in_indep_variables <- FALSE
  }

  # Callback for when data is to be processed from the queue.
  do.work <- function(metadata, data) {
    # Single column sample_id
    data <- deserialize(data)
    cases <- as.character(data[data[, "y"] == 1, "sample_id"])
    to_exclude <- as.character(data[is.na(data[, "y"]), "sample_id"])

    min_num_cases <- conf$binary_conf$min_num_cases

    covars$y <- as.numeric(covars$sample_id %in% cases)
    covars[covars$sample_id %in% to_exclude, "y"] <- NA

    n_excluded <- sum(is.na(covars[, "y"]))

    aug_model_rhs <- paste(c(base_cols, aug_cols), collapse = " + ")

    aug_formula <- as.formula(paste0(
      "y ~ ", aug_model_rhs
    ))

    aug_data_matrix <- model.frame(aug_formula, data = covars)
    keep <- complete.cases(aug_data_matrix)

    if (sum(keep) < min_num_cases) {
      return()
    }

    aug_data_matrix <- as.matrix(aug_data_matrix[keep, ])

    # if the sex column is defined and included as an independant variable.
    # We check that there is variability remaining after subsetting.
    if (sex_in_indep_variables) {
      n_unique_sex <- length(unique(aug_data_matrix[, conf$sex_column]))

      if (n_unique_sex <= 1) {
        # There is no variability in the sex (only men or women are being
        # analyzed.
        # So we remove this column from the matrix and from base_cols.
        aug_data_matrix <- aug_data_matrix[,
          (colnames(aug_data_matrix) != conf$sex_column)
        ]

        base_cols <- base_cols[base_cols != conf$sex_column]

        cat("R: Removing sex column because there is no variability.\n")
      }
    }

    n_cases <- sum(aug_data_matrix[, "y"] == 1)
    n_controls <- sum(aug_data_matrix[, "y"] == 0)

    fit_base <- fit_fastglm(aug_data_matrix[, c("y", base_cols)])
    fit_aug <- fit_fastglm(aug_data_matrix)

    lrt <- fastglm_anova(fit_base, fit_aug)

    line <- paste(
      metadata$variable_id,
      metadata$analysis_type,
      n_cases,
      n_controls,
      n_excluded,
      lrt$resid_deviance_base,
      lrt$resid_deviance_augmented,
      lrt$deviance,
      lrt$p,
      sep=","
    )

    writeLines(line, con = output_file)

  }

  # Start the main loop.
  Worker(worker_id, ..., callback=do.work)

  cat("R: closing output file.\n")
  close(output_file)

}

# do.call(linear_f_test_worker, list(worker_id, dealer_addr, monitor_addr))
logistic_deviance_diff_test_worker(worker_id, dealer_addr, monitor_addr)
