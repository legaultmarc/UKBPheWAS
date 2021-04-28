args <- commandArgs(trailingOnly=TRUE)

worker_id <- args[1]
dealer_addr <- args[2]
monitor_addr <- args[3]

continuous_descriptive_statistics_worker <- function(worker_id, ...) {
  # Read the analysis configuration object.
  conf <- fromJSON(file = "analysis_configuration.json")

  # Read the covariables.
  covars <- get_xs(conf)

  # Write header file.
  cat(paste0("variable_id,analysis_type,n_full_dataset,n_analysis_dataset\n"),
      file="header.csv")

  cat("R: Opened output file.\n")
  output_file <- file(
    paste0("results_worker_", worker_id, ".csv"),
    open = "wt"
  )

  do.work <- function(metadata, data) {
    data <- deserialize(data)

    n_full_dataset <- nrow(data)

    # Join with the xs.
    data <- merge(data, covars, by = "sample_id")
    data <- data[complete.cases(data), ]

    n_analysis_dataset <- nrow(data)

    line <- paste(
      metadata$variable_id,
      metadata$analysis_type,
      n_full_dataset,
      n_analysis_dataset,
      sep = ","
    )

    writeLines(line, con = output_file)

  }

  # Start the main loop.
  Worker(worker_id, ..., callback=do.work)

  cat("R: closing output file.\n")
  close(output_file)

}

# do.call(linear_f_test_worker, list(worker_id, dealer_addr, monitor_addr))
continuous_descriptive_statistics_worker(worker_id, dealer_addr, monitor_addr)
