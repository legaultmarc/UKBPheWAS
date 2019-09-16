library(DBI)

#' Get a connection to the UK Biobank database.
#'
#' @return A database connection object
#'
#' @examples
#' ukb_connect()
#'
#' @import RPostgres
#'
#' @export
get_ukb_connection <- function(password) {
  drv <- RPostgres::Postgres()
  con <- dbConnect(drv, dbname="uk_biobank_v2",
                   host="sripl-sg-stk01.statgen.local", user="uk_biobank_user",
                   password=password)

  con
}
