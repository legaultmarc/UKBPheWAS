#' Strip dots from the ICD10 code (e.g. I20.0 -> I20).
#'
#' @param icd10 The code.
#'
#' @return An integer with the level of the code as described above.
strip_icd <- function(icd10) {
  gsub("\\.", "", icd10)
}

#' Infer the level of an ICD10 code.
#'
#' Level 0 is the chapter (e.g. I)
#' Level 1 is for the next two numbers. (e.g. I20)
#' Level 2 is for codes that are more precise. (e.g. I23.8)
#'
#' @param icd10 The code.
#'
#' @return An integer with the level of the code as described above.
#'
#' @export
infer_icd10_level <- function(icd10) {
  # Strip dots.
  icd10 <- strip_icd(icd10)

  # Map length of the resulting vector to level.
  n <- nchar(icd10)

  if (n == 1) {
    return(0)
  }

  if (n == 3) {
    return(1)
  }

  if (n > 3) {
    return(2)
  }

  # Invalid code.
  return(NULL)
}

#' ICD10 code to group.
#'
#' The groups are sets of contiguous ICD10 codes. For example, I20-I25 
#' consists of ischemic heart diseases.
#' They need to be hard-coded somewhere so that codes can be grouped
#' automatically during a pheWAS.
#'
#' @param icd10 The code.
#'
#' @return A group name and a unique string representing the range as a list
#'         e.g. list(range="I20-I25", "Ischemic heart diseases")
#'
#' @export
get_icd10_group <- function(icd10) {
  groups <- list(

  )
}

#' Comparison operator for ICD10 codes.
#'
#' Utility function for ICD10 ordering so that it is possible to quickly
#' check if a code is within a range of codes.
#' This function is aware of the code structure (e.g. so that say B02 would be
#' greater than A05).
#' If a code is a subcode of another one (e.g. I20.0 and I20), the return value
#' will be 0.
#'
#' @param a The first ICD10 code.
#' @param b The second ICD10 code.
#'
#' @return -1 if a < b, 0 if a == b and 1 if a > b
compare_codes <- function(a, b) {

  # Strip dots from ICD10 codes.
  a <- strip_icd(a)
  b <- strip_icd(b)

  # Compare based on chapters (first character).
  chap_a <- substr(a, 1, 1)
  chap_b <- substr(b, 1, 1)

  if (chap_a < chap_b) return(-1)
  if (chap_a > chap_b) return(1)

  # Chapters are the same we need to compare codes.
  # We first truncate to the level of the broadest code.
  # At the same time, we remove the chapter (hence the -1 in the indexing).
  min_length <- min(nchar(a), nchar(b))
  a <- substr(a, 2, min_length)
  b <- substr(b, 2, min_length)

  if (a == b) return(0)
  if (a < b)
    return(-1)
  else
    return(1)

}

#' Tests if a given code is in an inclusive range of codes.
#'
#' @param code The code to test.
#' @param left The left bound of the inclusive range.
#' @param right The right bound of the inclusive range.
#'
#' @return T if the code is in the range, F otherwise
#'
#' @export
code_in_range <- function(code, left, right) {

  gte_than_left <- compare_codes(left, code) <= 0
  lte_than_right <- compare_codes(code, right) <= 0

  gte_than_left && lte_than_right

}
