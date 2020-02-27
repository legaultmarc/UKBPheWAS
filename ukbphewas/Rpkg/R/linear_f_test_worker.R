# R worker that listens to a ZMQ socket and runs analyses.
#
# They will all have the same command line interface
# All of the data distribution will be done using socket including an
# initialization and a high throughput phase.
#
# Shared covariables will be serialized to disk to be passed more easily.
# 
# What R needs:
# fit1 <- lm(y ~ covars)
# fit2 <- lm(y ~ covars + gene_pcs)
# to_json(anova(fit1, fit2))
#
# y is taken from ZMQ
# covars is taken from a serialized file.
# gene_pcs as well
# source("workers/worker.R")
library(arrow)

args <- commandArgs(trailingOnly=TRUE)

worker_id <- args[1]
dealer_addr <- args[2]
monitor_addr <- args[3]

linear_f_test_worker <- function(worker_id, ...) {
  # Read the analysis configuration object.
  conf <- fromJSON(file = "analysis_configuration.json")

  # Read the covariables.
  covars <- get_xs(conf)

  # Create a header file.
  cat(paste0("variable_id,analysis_type,n_samples,rss_base,rss_augmented,",
             "sum_of_sq,F_stat,p\n"),
      file="header.csv")

  cat("R: Opened output file.\n")
  output_file <- file(
    paste0("results_worker_", worker_id, ".csv"),
    open = "wt"
  )

  # Check if sex is in the independant variables and if the analysis is
  # stratified by sex.
  base_cols <- char_to_terms(conf$model_rhs)
  aug_cols <- conf$linear_conf$augmented_variables

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

  base_aug_cols <- c(base_cols, aug_cols)

  # Callback for when data is to be processed from the queue.
  do.work <- function(metadata, data) {
    data <- as.data.frame(read_table(data))

    # We explicitly treat the sample_ids as strings.
    data$sample_id <- as.character(data$sample_id)

    # We find indexers for the lhs and rhs.
    overlapping_samples <- intersect(covars$sample_id, data$sample_id)

    if (length(overlapping_samples) == 0) {
      cat("R: No overlapping samples, skipping.\n")
    }

    covars_idx <- match(overlapping_samples, covars$sample_id)
    data_idx <- match(overlapping_samples, data$sample_id)

    # Fit both models.
    outcome <- names(data)[names(data) != "sample_id"]

    m <- cbind(data[data_idx, outcome], covars[covars_idx, base_aug_cols])
    m <- as.matrix(m[complete.cases(m), ])

    # if the sex column is defined and included as an independant variable.
    # We check that there is variability remaining after subsetting.
    if (sex_in_indep_variables) {
      n_unique_sex <- length(unique(m[, conf$sex_column]))
      print(n_unique_sex)

      if (n_unique_sex <= 1) {
        # There is no variability in the sex (only men or women are being
        # analyzed.
        # So we remove this column from the matrix and from base_cols.
        m <- m[,
          (colnames(aug_data_matrix) != conf$sex_column)
        ]

        base_cols <- base_cols[base_cols != conf$sex_column]
        base_aug_cols <- base_aug_cols[base_aug_cols != conf$sex_column]

        cat("R: Removing sex column because there is no variability.\n")
      }
    }

    fit_base <- lm(m[, 1] ~ m[, base_cols])
    fit_aug <- lm(m[, 1] ~ m[, base_aug_cols])

    f <- anova(fit_base, fit_aug)

    # return(data.frame(
    #   outcome_id = data$id,
    #   outcome_label = data$label,
    #   rss_base = f[1, 2],
    #   rss_augmented = f[2, 2],
    #   sum_of_sq = f[2, 4],
    #   F_stat = f[2, 5],
    #   p = f[2, 6]
    # ))
    line <- paste(
      metadata$variable_id,
      metadata$analysis_type,
      nrow(m),  # n_samples
      f[1, 2],  # rss base
      f[2, 2],  # rss aug
      f[2, 4],  # ssq
      f[2, 5],  # F
      f[2, 6],  # p
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
linear_f_test_worker(worker_id, dealer_addr, monitor_addr)
