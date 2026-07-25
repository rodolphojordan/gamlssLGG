#' Log-Generalized Gamma Density Function
#'
#' @param y Vector of quantiles
#' @param mu Location parameter
#' @param sigma Scale parameter (positive)
#' @param nu Shape parameter
#' @param log Logical; if TRUE, probabilities are given as log
#'
#' @return Density values
#' @export
#'
#' @examples
#' dLGG(0, mu = 0, sigma = 1, nu = 0.5)
#' dLGG(0:2, mu = 0, sigma = 1, nu = 0.5, log = TRUE)
dLGG <- function(y,
                 mu = 0,
                 sigma = 1,
                 nu = 0.5,
                 log = FALSE) {

  sigma <- pmax(sigma, 1e-10)

  # Caso limite: nu -> 0 converge para Normal(mu,sigma²)
  if (any(abs(nu) < 1e-06)) {
    return(
      dnorm(
        y,
        mean = mu,
        sd = sigma,
        log = log
      )
    )
  }

  # Constante da densidade
  logC <- logC_LGG(nu)

  # Log-densidade
  eta <- (nu / sigma) * (y - mu)
  eta <- pmin(eta, 700)
  eta <- pmax(eta, -700)

  z <- exp(eta)

  logfy <-
    logC -
    log(sigma) +
    (y - mu) / (nu * sigma) -
    (1 / nu^2) * z

  if (log) {
    return(logfy)
  }

  exp(logfy)
}
