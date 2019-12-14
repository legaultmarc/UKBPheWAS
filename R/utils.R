library(RPostgres)


#' Get all hospitalization and/or death records from the UKB databse.
#'
#' @import RPostgres
get_full_records <- function(
  con,
  include_secondary_hospit,
  include_death_records
) {

  sql <- paste0(
    "select distinct ",
    " eid, diag_icd10 ",
    internal_get_from_stmt(include_secondary_hospit, include_death_records),
    "where diag_icd10 is not null"
  )

  query <- dbSendQuery(con, sql)
  res <- dbFetch(query)
  dbClearResult(query)

  # Make sure diag_icd10 are strings.
  res$diag_icd10 <- as.character(res$diag_icd10)

  res
}


#' Get curated cv_endpoints.
get_cv_endpoints <- function(con) {

  query <- dbSendQuery(con, query_get_cv_endpoints_)
  res <- dbFetch(query)
  dbClearResult(query)

  res

}


#' Create a SQL string with an appropriate FROM statement to obtain cases.
#'
#' Specifically the generated table will have two columns: eid and diag_icd10
#' from the appropriate HESIN table (primary or primary+secondary) and with
#' or without the death records.
#'
#' @param include_secondary_hospit Include the secondary hospitalization codes.
#' @param include_death_records Include the causes of death.
#'
#' @return A string containing the FROM statement
internal_get_from_stmt <- function(
  include_secondary_hospit,
  include_death_records
) {

  # We use a different table for the hospitalization data depending on the
  # use of primary or primary+secondary hospitalization codes.
  hes_table <- ifelse(
    include_secondary_hospit,
    "recurrent_events.full_hesin",
    "hesin"
  )

  sql <- NULL
  if (include_death_records) {
    # We union on the causes of death if requested.
    sql <- paste0(
      "from ( ",
      "  select eid, diag_icd10 from ", hes_table,
      "  union ",
      "  select sample_id, cause from cv_endpoints.death_records ",
      ") t "
    )
  }

  else {
    # Otherwise we don't use the cause of death.
    sql <- paste0(
      "from ", hes_table, " "
    )
  }

  return(sql)
}

#' Get the age and sex from the same dataset as the HCN4 paper.
#'
#' @param con A connection to the UK Biobank
#'
#' @return A dataframe of covariables.
#' @export
get_covariables_hcn4_project <- function(con) {
  query <- dbSendQuery(
    con, "select sample_id, male, age from legaultm_hcn4_paper.final_dataset"
  )
  res <- dbFetch(query)
  dbClearResult(query)

  res
}


get_cancer_data <- function(con) {
  # Get cancer ICD10 codes.
  query <- dbSendQuery(
    con,
    paste0("select distinct sample_id as eid, value as diag_icd10 ",
           "from variable_categorical where variable_id = 40006")
  )
  cases_info <- dbFetch(query)
  dbClearResult(query)

  # Also find individuals with no ICD10 that self-report a cancer or that
  # have a cancer ICD9 code and return them as cancer_excl_from_controls.
  query <- dbSendQuery(
    con,
    paste0(
      "select distinct sample_id ",
      "from variable_categorical where variable_id = 40013 ",  # ICD9
      "union ",
      "select sample_id ",
      "from variable_categorical where variable_id = 20001"  # Self-reported
    )
  )

  excl_from_ctrls = dbFetch(query)$sample_id
  dbClearResult(query)

  excl_from_ctrls <- excl_from_ctrls[
    !(excl_from_ctrls %in% cases_info$eid)
  ]

  warning(paste0(
    "Recorded ", length(excl_from_ctrls), " individuals to remove from ",
    "cancer controls (because of ICD9 code or self-report)."
  ))

  list(
    cancer_case_data = cases_info,
    cancer_excl_from_controls = excl_from_ctrls
  )
}


extract_raw_data <- function(configuration) {

  con <- get_ukb_connection(configuration$db_password)

  data <- list()

  # Extract information on diseases if required.
  if (should_do_binary(configuration)) {
    diseases <- get_full_records(
      con,
      configuration$binary_configuration$include_secondary_hospit,
      configuration$binary_configuration$include_death_records
    )

    # We exclude some chapters corresponding to diseases that are
    # 'external' or cancer as we rely on the cancer registry instead.
    excluded_chapters <- c("C", "D", "S", "T", "U", "V", "W", "X", "Y", "Z")
    diseases <- diseases[
      !(substr(diseases$diag_icd10, 1, 1) %in% excluded_chapters),
    ]

    data$diseases <- diseases

    # Also extract information on cancer data from the cancer registry.
    cancer <- get_cancer_data(con)

    # Merge the ICD10 codes.
    data$diseases <- rbind(data$diseases, cancer$cancer_case_data)

    # But also remember the exclusions from controls.
    data$cancer_excl_from_controls <- cancer$cancer_excl_from_controls

    # Last piece of data is for the manually defined cv_endpoints:
    data$cv_endpoints <- get_cv_endpoints(con)

  }

  # Extract information on continuous traits.
  if (should_do_linear(configuration)) {
    data$continuous <- read.csv(
      paste0(
        configuration$continuous_variables_path,
        "/transformed_qt_traits.csv.gz"
      ),
      stringsAsFactors = FALSE
    )

    data$continuous_metadata <- read.csv(
      paste0(
        configuration$continuous_variables_path,
        "/transforms_metadata.csv"
      ),
      stringsAsFactors = FALSE
    )
  }

  dbDisconnect(con)

  data

}


# Keep only rows of interest from output and save to disk.
clean_and_save <- function(analysis_label, results, configuration) {
  if (!is.null(configuration$results_filter)) {
    results <- results[
      sapply(1:nrow(results), function(i) {
        configuration$results_filter(i, results[i, ])
      }),
    ]
  }

  if (nrow(results) == 0) {
    warning("No results selected by variable of interest filter.")
  }

  else {
    filename <- paste0(configuration$output_prefix, "_", analysis_label, ".csv")
    cat(paste0("WRITING results file: ", filename, "\n"))

    write.csv(results, filename, row.names = FALSE)
  }
}
