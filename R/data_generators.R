# Others are raw and blocks.
# Generators call callbacks for every result and should implement the
# parallelism.
generator_icd10_three_chars <- function(
  configuration, raw_data, callback, cl, limit=NULL
) {
  # Data frame with "eid" and "diag_icd10"
  data <- raw_data$diseases

  # Get list of 3 char codes.
  codes <- unique(substr(data$diag_icd10, 1, 3))

  if (!is.null(limit)) {
    warning(paste0(
      "Limiting the number of 3 character codes to ", limit, ". This should ",
      "not be used routinely, mostly useful for testing."
    ))
    codes <- codes[1:limit]
  }

  icd10 <- UKBPheWAS::icd10

  results <- foreach (
    code = codes,
    .combine = "rbind",
    .packages = c("UKBPheWAS", "broom")
  ) %dopar% {

    # Identify cases.
    # It's more efficient to only represent cases, so we send that to the
    # callback who should code controls if needed.
    cases <- unique(data[substr(data$diag_icd10, 1, 3) == code, "eid"])

    # Try to find a label in the metadata.
    label <- icd10[icd10$coding_key == code, "meaning"]
    if (identical(label, character(0))) {
        label <- code
    }

    cur_data <- list(
      id = code,
      label = label,
      y = data.frame(sample_id = cases)
    )

    callback(configuration, cur_data)

  }

  return(results)
}


generator_icd10_raw <- function(
  configuration, raw_data, callback, cl, limit=NULL
) {
  # Data frame with "eid" and "diag_icd10"
  data <- raw_data$diseases

  # Get list of 3 char codes.
  codes <- unique(data$diag_icd10)

  if (!is.null(limit)) {
    warning(paste0(
      "Limiting the number of raw ICD10 codes to ", limit, ". This should ",
      "not be used routinely, mostly useful for testing."
    ))
    codes <- codes[1:limit]
  }

  icd10 <- UKBPheWAS::icd10

  results <- foreach (
    code = codes,
    .combine = "rbind",
    .packages = c("UKBPheWAS", "broom")
  ) %dopar% {
    code_len <- nchar(code)
    cases <- unique(data[substr(data$diag_icd10, 1, code_len) == code, "eid"])

    label <- icd10[icd10$coding_key == code, "meaning"]
    if (identical(label, character(0))) {
        label <- code
    }

    cur_data <- list(
      id = code,
      label = label,
      y = data.frame(sample_id = cases)
    )

    callback(configuration, cur_data)

  }

  return(results)
}


generator_icd10_blocks <- function(
  configuration, raw_data, callback, cl, limit=NULL
) {
  # Data frame with "eid" and "diag_icd10"
  data <- raw_data$diseases

  blocks <- UKBPheWAS::icd10_blocks

  if (!is.null(limit)) {
    warning(paste0(
      "Limiting the number of raw ICD10 odes to ", limit, ". This should ",
      "not be used routinely, mostly useful for testing."
    ))
  }

  max_idx <- min(limit, nrow(blocks))

  results <- foreach (
    i = 1:max_idx,
    .combine = "rbind",
    .packages = c("UKBPheWAS", "broom")
  ) %dopar% {

    block <- blocks[i, ]

    cases <- data[
      sapply(data$diag_icd10, function(code) {
        code_in_range(code, block$left, block$right)
      }),
    ]["eid"]

    cur_data <- list(
      id = paste0(block$left, "-", block$right),
      label = block$block,
      y = data.frame(sample_id = cases)
    )

    callback(configuration, cur_data)

  }

  return(results)

}


# All generators should yield (id,label,y data frame [sample_id, y])
generator_standardized_continuous <- function(
  configuration, raw_data, callback, cl, limit=NULL
) {

  n_pheno <- nrow(raw_data$continuous_metadata)

  if (!is.null(limit)) {
    warning(paste0(
      "Limiting the number of continuous phenotypes to ", limit, ". This should ",
      "not be used routinely, mostly useful for testing."
    ))
    n_pheno <- limit
  }

  results <- foreach (
    i = 1:n_pheno,
    .combine = "rbind",
    .packages = c("UKBPheWAS", "broom")
  ) %dopar% {

    # Current metadata row.
    meta <- raw_data$continuous_metadata[i, ]

    # Extract data.
    y <- raw_data$continuous[
      raw_data$continuous$variable == paste0("v", meta$ukbphewas_id),
      c("sample_id", "value")
    ]

    names(y) <- c("sample_id", "y")

    # Standardize
    y$y <- (y$y - mean(y$y)) / sd(y$y)

    out <- list(
      id = paste0("cont_v", meta$ukbphewas_id),
      label = meta$variable,
      y = y
    )

    callback(configuration, out)

  }
  
  return(results)

}
