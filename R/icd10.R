#' Check if an array of ICD10 codes fall within a range.
#'
#' @export
codes_in_range <- function(codes, left, right) {

    # Strip dots.
    strip_dots <- function(x) {
        gsub("\\.", "", x)
    }

    left <- strip_dots(left)
    right <- strip_dots(right)

    stopifnot(length(left) == length(right))
    block_bound_length <- nchar(left)

    # Set dots.
    parse_icd10 <- function(x, toLength = NULL) {
        x <- strip_dots(x)

        if (!is.null(toLength)) {
          x <- substr(x, 1, toLength)
        }

        chapter <- substr(x, 1, 1)
        num <- substring(x, 2)

        tens <- substr(num, 1, 2)
        decimal <- substring(num, 3)

        num = as.numeric(paste0(tens, ".", decimal))

        list(chapter = chapter, num = num)
    }

    codes <- parse_icd10(codes, block_bound_length)
    left <- parse_icd10(left)
    right <- parse_icd10(right)

    # Greater (or equal) than left bound.
    gt_left <- (
        # Either chapter is greater
        (codes$chapter > left$chapter) |

        # Or chapter is the same, but number is greater or equal.
        ((codes$chapter == left$chapter) & (codes$num >= left$num))
    )

    # Less than (or equal) than right bound.
    lt_right <- (
        (codes$chapter < right$chapter) |
        ((codes$chapter == right$chapter) & (codes$num <= right$num))
    )

    gt_left & lt_right
}

#' Check if ICD10 code is a cancer code.
is_cancer_code <- function(code) {
  chapter <- substr(code, 1, 1)
  num <- as.numeric(substr(code, 2, 3))

  chapter == "C" || (chapter == "D" && num <= 48)
}
