#' Incomplete beta function B(z; a, b)
#' @noRd
iBeta <- function(z, a, b) {
   pbeta(z, a, b) * beta(a, b)
}
