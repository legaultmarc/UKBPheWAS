library(doParallel)
library(parallel)
library(foreach)


bin_data_generators <- list(
  self_reported_diseases = list(
    gen = generator_self_reported_diseases
  ),
  cv_endpoints = list(
    gen = generator_cv_endpoints
  ),
  icd10_three_chars = list(
    gen = generator_icd10_three_chars
  ),
  icd10_blocks = list(
    gen = generator_icd10_blocks
  ),
  icd10_raw = list(
    gen = generator_icd10_raw
  )
)


#' Main method to run pheWAS
#'
#' @import parallel
#' @import doParallel
#' @export runPheWAS
runPheWAS <- function(configuration, build_cache = FALSE, raw_cache = NULL) {

  # Log the configuration summarizing the analysis.
  log_configuration(configuration)

  # Load or create the data cache if needed.
  if (is.null(raw_cache)) {
    raw <- extract_raw_data(configuration)

    if (build_cache) {
      saveRDS(raw, file = "_UKBPheWAS_cache.rds.gz", compress = "gzip")
    }
  }
  else {
    cat("Using cached raw data.\n")
    raw <- raw_cache
  }

  cl <- makeForkCluster(configuration$ncpus)
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

    # Skip the raw analysis of ICD10 codes if needed.
    # For now we can only "easily" skip this because it is the heaviest
    # computationally and the noisiest in terms of interpretations.
    if (configuration$binary_configuration$skip_icd10_raw) {
      bin_data_generators <- bin_data_generators[
        names(bin_data_generators) != "icd10_raw"
      ]
    }

    for (gen_name in names(bin_data_generators)) {

      cat(paste0("Running analysis for ", gen_name, "\n"))

      results <- bin_data_generators[[gen_name]]$gen(
        configuration, raw, configuration$binary_configuration$callback, cl,
        limit = configuration$limit
      )

      cat("Got results, clean and save\n")

      results <- clean_and_save(gen_name, results, configuration)

      remove(results)
      gc(full = TRUE)

      cat("DONE!\n")

    }

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
      limit = configuration$limit
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

  if (!is.null(config$db_password)) {
    config$db_password <- "**********"
  }

  print(str(config))

}
