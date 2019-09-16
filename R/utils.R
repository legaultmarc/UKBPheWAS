library(RPostgres)

#' Get a vec of icd10 codes with at least n occurences in hospitalization data.
#'
#' This relies on the "full hesin" table.
#'
#' @param n The minimum number of cases for inclusion.
#'
#' @return A list of ICD10 codes.
#'
#' @export
get_icd10_with_at_least_n_cases <- function(
  con,
  n,
  include_secondary_hospit=T,
  include_death_records=T
)
{

  sql <- paste0(
    "select ",
    " diag_icd10 ",
    internal_get_from_stmt(include_secondary_hospit, include_death_records),
    "where diag_icd10 is not null ",
    "group by diag_icd10 ",
    "having count(distinct eid) >= $1 "
  )

  query <- dbSendQuery(con, sql)
  dbBind(query, list(n))

  df <- dbFetch(query)
  dbClearResult(query)
  df$diag_icd10
}


get_full_records <- function(
  con,
  include_secondary_hospit,
  include_death_records
) {

  sql <- paste0(
    "select ",
    " eid, diag_icd10 ",
    internal_get_from_stmt(include_secondary_hospit, include_death_records),
    "where diag_icd10 is not null"
  )

  query <- dbSendQuery(con, sql)
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
)
{

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
