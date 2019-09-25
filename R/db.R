library(DBI)

#' Get a connection to the UK Biobank database.
#'
#' @return A database connection object
#'
#' @examples
#' \dontrun{
#' get_ukb_connection("my_password")
#' }
#'
#' @import RPostgres
#'
#' @export
get_ukb_connection <- function(password) {
  drv <- RPostgres::Postgres()
  con <- dbConnect(drv, dbname = "uk_biobank_v2",  # nolint
                   host = "sripl-sg-stk01.statgen.local",
                   user = "uk_biobank_user",
                   password = password)

  con
}


#' Get a connection to the UK Biobank biomarkers temporary DB.
#'
#' @return A database connection object
#'
#' @examples
#' \dontrun{
#' get_biomarker_connection("my_password")
#' }
#'
#' @import RPostgres
#'
#' @export
get_biomarker_connection <- function(password) {
  drv <- RPostgres::Postgres()
  con <- dbConnect(drv, dbname = "tmp_uk_biobank_lipids",
                   host = "sripl-sg-stk01.statgen.local",
                   user = "uk_biobank_user",
                   password = password)

  con
}
