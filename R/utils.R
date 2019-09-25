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


#' Get all biomarker data in the long format (averaged over visits if needed).
#'
#' @import RPostgres
get_all_biomarkers <- function(con) {

  sql <- paste0(
    "select ",
    "  var.sample_id, ",
    "  var.variable_id, ",
    "  string_agg(distinct meta.description, ', ') as description, ",
    "  avg(var.value) as val ",
    "from variable_float var right outer join ( ",
    "  select * from variable_metadata ",
    "  where data_type='float' and coding_id is NULL ",
    ") meta on var.variable_id=meta.variable_id ",
    "group by var.variable_id, var.sample_id; "
  )

  query <- dbSendQuery(con, sql)
  res <- dbFetch(query)
  dbClearResult(query)

  res

}
