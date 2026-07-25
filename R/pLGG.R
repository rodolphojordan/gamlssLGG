#' Log-Generalized Gamma Distribution Function
#'
#' @param q Vector of quantiles
#' @param mu Location parameter
#' @param sigma Scale parameter (positive)
#' @param nu Shape parameter
#' @param lower.tail Logical; if TRUE, probabilities are P[X <= x]
#' @param log.p Logical; if TRUE, probabilities are given as log
#'
#' @return Cumulative distribution function values
#' @export
#'
#' @examples
#' pLGG(0, mu = 0, sigma = 1, nu = 0.5)
#' pLGG(c(-1, 0, 1), mu = 0, sigma = 1, nu = 0.5)
pLGG <- function(q,
                 mu = 0,
                 sigma = 1,
                 nu = 0.5,
                 lower.tail = TRUE,
                 log.p = FALSE) {

  sigma <- pmax(sigma, 1e-10)

  if (all(abs(nu) < 1e-06)) {
    return(
      pnorm(
        q,
        mean = mu,
        sd = sigma,
        lower.tail = lower.tail,
        log.p = log.p
      )
    )
  }

  eta <- (nu / sigma) * (q - mu)
  eta <- pmin(eta, 700)
  eta <- pmax(eta, -700)

  z <- (1 / nu^2) * exp(eta)

  pgamma(
    z,
    shape = nu^(-2),
    rate = 1,
    lower.tail = ifelse(nu > 0,
                        lower.tail,
                        !lower.tail),
    log.p = log.p
  )
}
