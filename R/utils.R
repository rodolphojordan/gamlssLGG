#' Internal function for the log-normalizing constant
#'
#' @param nu Shape parameter
#' @keywords internal
#' @export
logC_LGG <- function(nu) {
  log(abs(nu)) -
    lgamma(nu^(-2)) +
    nu^(-2) * log(nu^(-2))
}
