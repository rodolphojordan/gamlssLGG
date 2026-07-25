#' Log-Generalized Gamma Quantile Function
#'
#' @param p Vector of probabilities
#' @param mu Location parameter
#' @param sigma Scale parameter (positive)
#' @param nu Shape parameter
#' @param lower.tail Logical; if TRUE, probabilities are P[X <= x]
#' @param log.p Logical; if TRUE, probabilities are given as log
#'
#' @return Quantile values
#' @export
#'
#' @examples
#' qLGG(0.5, mu = 0, sigma = 1, nu = 0.5)
#' qLGG(c(0.25, 0.5, 0.75), mu = 0, sigma = 1, nu = 0.5)
qLGG <- function(p,
                 mu = 0,
                 sigma = 1,
                 nu = 0.5,
                 lower.tail = TRUE,
                 log.p = FALSE) {

  sigma <- pmax(sigma, 1e-10)

  if (all(abs(nu) < 1e-06)) {
    return(
      qnorm(
        p,
        mean = mu,
        sd = sigma,
        lower.tail = lower.tail,
        log.p = log.p
      )
    )
  }

  if (log.p)
    p <- exp(p)

  if (!lower.tail)
    p <- 1 - p

  z <- qgamma(
    p,
    shape = nu^(-2),
    rate = 1,
    lower.tail = (nu > 0)
  )

  mu + (sigma / nu) * log(nu^2 * z)
}
