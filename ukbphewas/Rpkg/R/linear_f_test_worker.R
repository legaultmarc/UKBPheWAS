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

args <- commandArgs(trailingOnly=TRUE)

worker_id <- args[1]
dealer_addr <- args[2]
monitor_addr <- args[3]


linear_f_test_worker <- function(worker_id, ...) {
  # Read the analysis configuration object.
  conf <- fromJSON(file = "analysis_configuration.json")

  # Read the covariables.
  covars <- get_xs(conf)

  # Read the XPCs.

  xpcs <- read.csv(
    conf$linear_conf$xpcs_path,
    colClasses=read_csv_filter_columns(
      conf$linear_conf$xpcs_path, 
      c("sample_id", conf$linear_conf$augmented_variables)
    )
  )
  covars <- merge(covars, xpcs, by = "sample_id")

  # For now, sample_id needs to be a string. Eventually, we will infer this
  # better.
  covars$sample_id <- as.character(covars$sample_id)

  # Create a header file.
  cat(paste0("variable_id,analysis_type,n_samples,rss_base,n_params_base,",
             "rss_augmented,n_params_aug,sum_of_sq,F_stat,p\n"),
      file="header_summary.csv")

  cat("R: Opened output file.\n")
  output_file <- file(
    paste0("results_worker_", worker_id, "_summary.csv"), open = "wt"
  )
  model_file <- file(
    paste0("results_worker_", worker_id, "_model.json"), open = "wt"
  )

  base_cols <- char_to_terms(conf$model_rhs)
  aug_cols <- conf$linear_conf$augmented_variables
  base_aug_cols <- c(base_cols, aug_cols)

  # Callback for when data is to be processed from the queue.
  do.work <- function(metadata, data) {
    data <- deserialize(data)

    # We find indexers for the lhs and rhs.
    if (class(covars$sample_id) != class(data$sample_id)) {
      cat(paste0("Covars sample_id class: ", class(covars$sample_id), "\n"))
      cat(paste0("Data generator sample_id class: ", class(data$sample_id), "\n"))
      stop("sample_id type mismatch between covars and data generators.")
    }
    overlapping_samples <- intersect(covars$sample_id, data$sample_id)

    if (length(overlapping_samples) == 0) {
      cat("R: No overlapping samples, skipping.\n")
    }

    covars_idx <- match(overlapping_samples, covars$sample_id)
    data_idx <- match(overlapping_samples, data$sample_id)

    # Fit both models.
    outcome <- names(data)[names(data) != "sample_id"]

    m <- cbind(data[data_idx, outcome], covars[covars_idx, base_aug_cols])
    m <- m[complete.cases(m), ]

    # Find columns on the right hand side that have no variance after
    # joining.
    col_drop_li <- drop_columns_with_no_variance(m)
    m <- col_drop_li$mat
    dropped_cols <- col_drop_li$dropped_cols

    if (length(dropped_cols) != 0) {
      base_cols <- base_cols[!(base_cols %in% dropped_cols)]
      base_aug_cols <- base_aug_cols[!(base_aug_cols %in% dropped_cols)]
    }

    ylab <- colnames(m)[1]

    # I'd rather use the formula API so that the column names are nicer in the
    # output.
    base_formula <- as.formula(paste0(
      ylab, " ~ ", paste0(base_cols, collapse = " + ")
    ))
    fit_base <- lm(base_formula, data = m)

    aug_formula <- as.formula(paste0(
      ylab, " ~ ", paste0(base_aug_cols, collapse = " + ")
    ))
    fit_aug <- lm(aug_formula, data = m)

    f <- anova(fit_base, fit_aug)

    # Save the parameters as json.
    infer_df <- as.data.frame(summary(fit_aug)$coefficients)
    names(infer_df) <- c("beta", "se", "t", "p")
    infer_df <- cbind(data.frame(term = rownames(infer_df)), infer_df)
    rownames(infer_df) <- NULL

    # Calculate nlog10p approximate using chi2
    infer_df$nlog10p <- pchisq(
      (infer_df$beta / infer_df$se) ** 2,
      1,
      lower.tail = F,
      log.p = T
    ) / -log(10)

    line <- paste(
      metadata$variable_id,
      metadata$analysis_type,
      nrow(m),  # n_samples
      f[1, 2],  # rss base
      length(base_cols) + 1,  # n parameters base (+1 for intercept)
      f[2, 2],  # rss aug
      length(base_aug_cols) + 1,  # n parameters base (+1 for intercept)
      f[2, 4],  # ssq
      f[2, 5],  # F
      f[2, 6],  # p
      sep=","
    )

    writeLines(line, con = output_file)

    # Coefficients for the augmented model columns.
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
  Worker(worker_id, ..., callback = do.work)

  cat("R: closing output file.\n")
  close(output_file)
  close(model_file)

}

# do.call(linear_f_test_worker, list(worker_id, dealer_addr, monitor_addr))
linear_f_test_worker(worker_id, dealer_addr, monitor_addr)
