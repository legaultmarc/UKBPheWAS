library(doParallel)
library(foreach)


# Main method to run pheWAS
runPheWAS2 <- function(configuration) {

  log_configuration(configuration)

  raw <- extract_raw_data(configuration)

  cl <- makeCluster(configuration$ncpus)
  registerDoParallel(cl)

  # Models for binary outcomes.
  if (should_do_binary(configuration)) {
    # 3 character codes logistic.
    results <- generator_icd10_three_chars(configuration, raw, do_logistic, cl,
                                           limit = 20)
    results <- clean_and_save("3chars", results, configuration)
  }
  else {
    cat("SKIPPING models for binary outcomes.\n")
  }

  if (!is(configuration$linear_configuration, "Skip")) {
    # TODO
  }
  else {
    cat("SKIPPING models for linear outcomes.\n")
  }

  stopCluster(cl)

}


log_configuration <- function(config) {

  print(str(config))

}
