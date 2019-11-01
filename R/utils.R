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


extract_raw_data <- function(
  configuration,
  continuous_variables_path = "/data/projects/uk_biobank/data/reports/SGR-2094/"
) {

  con <- get_ukb_connection(configuration$db_password)

  data <- list()

  # Extract information on diseases if required.
  if (should_do_binary(configuration)) {
    data$diseases <- get_full_records(
      con,
      configuration$binary_configuration$include_secondary_hospit,
      configuration$binary_configuration$include_death_records
    )
  }

  if (should_do_linear(configuration)) {
    data$continuous <- read.csv(
      paste0(continuous_variables_path, "/transformed_qt_traits.csv.gz"),
      stringsAsFactors = FALSE
    )

    data$continuous_metadata <- read.csv(
      paste0(continuous_variables_path, "/transforms_metadata.csv"),
      stringsAsFactors = FALSE
    )
  }

  dbDisconnect(con)

  data

}


# Keep only rows of interest from output and save to disk.
clean_and_save <- function(analysis_label, results, configuration) {
  print(str(results))
  # TODO configuration$results_filter undefined.
  # probably change to results instead of voi
  results <- results[
    sapply(1:nrow(results), function(i) {
      configuration$results_filter(i, results[i, ])
    }),
  ]

  if (nrow(results) == 0) {
    warning("No results selected by variable of interest filter.")
  }

  else {
    write.csv(
      results,
      paste0(configuration$output_prefix, "_", analysis_label, ".csv")
    )
  }
}
