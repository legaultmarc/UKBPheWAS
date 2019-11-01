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

  # Get ICD10 metadata to have meaningful labels.
  icd10 <- read.csv(
    here("../../data/icd10/icd10_meaning.csv"),
    stringsAsFactors = FALSE
  )

  # foreach (
  #   code = codes,
  #   .combine = "rbind",
  #   .packages = c("UKBPheWAS", "broom")
  # ) %dopar% {
  results <- NULL
  for (code in codes) {

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

    # callback(configuration, cur_data)
    results <- rbind(results, callback(configuration, cur_data))

  }
  results

}


# All generators should yield (id,label,y data frame [sample_id, y])
# generator_standardized_continuous <- function(raw_data) {
# 
#   i <- 1
# 
#   # Generator object of standardized continuous variables.
#   yielder <- function() {
# 
#     if (i > nrow(raw_data$continuous_metadata)) {
#       return(NULL)
#     }
# 
#     # Current metadata row.
#     meta <- raw_data$continuous_metadata[i, ]
# 
#     # Extract data.
#     y <- raw_data$continuous[
#       raw_data$continuous$variable == paste0("v", meta$ukbphewas_id),
#       c("sample_id", "value")
#     ]
# 
#     names(y) <- c("sample_id", "y")
# 
#     # Standardize
#     y$y <- (y$y - mean(y$y)) / sd(y$y)
# 
#     out <- list(
#       id = paste0("cont_v", meta$ukbphewas_id),
#       label = meta$variable,
#       y = y
#     )
# 
#     i <<- i + 1
# 
#     return(out)
# 
#   }
#   
#   return(yielder)
# 
# }
