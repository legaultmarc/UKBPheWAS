args <- commandArgs(trailingOnly=TRUE)

worker_id <- args[1]
dealer_addr <- args[2]
monitor_addr <- args[3]

continuous_descriptive_statistics_worker <- function(worker_id, ...) {
  # Read the analysis configuration object.
  conf <- fromJSON(file = "analysis_configuration.json")

  # Read the covariables.
  covars <- get_xs(conf)

  do.work <- function(metadata, data) {
    data <- deserialize(data)

    # Join with the xs.
    data <- merge(data, covars, by = "sample_id")
    data <- data[complete.cases(data), ]

    cat(paste0("I work hard ", metadata$variable_id, ", ", nrow(data), "\n"))

  }

  # Start the main loop.
  Worker(worker_id, ..., callback=do.work)

  cat("R: closing output file.\n")
  # close(output_file)

}

# do.call(linear_f_test_worker, list(worker_id, dealer_addr, monitor_addr))
continuous_descriptive_statistics_worker(worker_id, dealer_addr, monitor_addr)
