library(doParallel)
library(foreach)


#' Main method to run pheWAS
#'
#' @export runPheWAS2
runPheWAS2 <- function(configuration, raw_cache = NULL) {

  log_configuration(configuration)

  limit = 20

  if (is.null(raw_cache)) {
    raw <- extract_raw_data(configuration)
  }
  else {
    cat("Using cached raw data.\n")
    raw <- raw_cache
  }

  cl <- makeCluster(configuration$ncpus)
  registerDoParallel(cl)

  # Models for binary outcomes.
  if (should_do_binary(configuration)) {

    if(!is(configuration$binary_configuration, "BinaryConfig")) {
      stop(paste0(
        "Invalid binary configuration of class: ",
        class(configuration$binary_configuration), ". ",
        "Expected BinaryConfig"
      ))
    }

    # 3 character codes logistic.
    cat("Running analysis based on 3 character codes...\n")
    results <- generator_icd10_three_chars(
      configuration, raw, configuration$binary_configuration$callback, cl,
      limit = limit
    )
    results <- clean_and_save("3chars", results, configuration)
    cat("DONE!\n")

    # Raw codes
    cat("Running analysis based on raw ICD10 codes...\n")
    results <- generator_icd10_raw(
      configuration, raw, configuration$binary_configuration$callback, cl,
      limit = limit
    )
    results <- clean_and_save("raw", results, configuration)
    cat("DONE!\n")

    # ICD10 blocks
    cat("Running analysis based on ICD10 blocks...\n")
    results <- generator_icd10_blocks(
      configuration, raw, configuration$binary_configuration$callback, cl,
      limit = limit
    )
    results <- clean_and_save("blocks", results, configuration)
    cat("DONE!\n")

  }

  else {
    cat("SKIPPING models for binary outcomes.\n")
  }

  # Models for linear outcomes.
  if (!is(configuration$linear_configuration, "Skip")) {

    if(!is(configuration$linear_configuration, "LinearConfig")) {
      stop(paste0(
        "Invalid linear configuration of class: ",
        class(configuration$linear_configuration), ". ",
        "Expected LinearConfig"
      ))
    }

    cat("Running analysis for continuous variables...\n")
    results <- generator_standardized_continuous(
      configuration, raw, configuration$linear_configuration$callback, cl,
      limit = limit
    )
    results <- clean_and_save("linear", results, configuration)
    cat("DONE\n")

  }

  else {
    cat("SKIPPING models for linear outcomes.\n")
  }

  stopCluster(cl)

}


log_configuration <- function(config) {

  print(str(config))

}
