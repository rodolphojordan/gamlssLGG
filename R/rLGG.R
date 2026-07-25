#' Random Generation for the Log-Generalized Gamma Distribution
#'
#' @param n Number of observations
#' @param mu Location parameter
#' @param sigma Scale parameter (positive)
#' @param nu Shape parameter
#'
#' @return Vector of random numbers
#' @export
#'
#' @examples
#' rLGG(10, mu = 0, sigma = 1, nu = 0.5)
rLGG <- function(n,
                 mu = 0,
                 sigma = 1,
                 nu = 0.5) {
  u <- runif(n)
  qLGG(
    u,
    mu = mu,
    sigma = sigma,
    nu = nu
  )
}
