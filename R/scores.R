#' Analytical Score for mu
#'
#' @param y Observations
#' @param mu Location parameter
#' @param sigma Scale parameter
#' @param nu Shape parameter
#'
#' @return Score values
#' @export
score_mu <- function(y, mu, sigma, nu) {
  z <- exp((nu / sigma) * (y - mu))
  (z - 1) / (nu * sigma)
}

#' Analytical Score for sigma
#'
#' @param y Observations
#' @param mu Location parameter
#' @param sigma Scale parameter
#' @param nu Shape parameter
#'
#' @return Score values
#' @export
score_sigma <- function(y, mu, sigma, nu) {
  w <- y - mu
  z <- exp((nu / sigma) * w)
  -1 / sigma + w * (z - 1) / (nu * sigma^2)
}

#' Analytical Score for nu
#'
#' @param y Observations
#' @param mu Location parameter
#' @param sigma Scale parameter
#' @param nu Shape parameter
#'
#' @return Score values
#' @export
score_nu <- function(y, mu, sigma, nu) {
  w <- y - mu
  z <- exp((nu / sigma) * w)

  # Derivada de log(C(nu))
  dlogC <- 1 / nu + 2 * nu^(-3) * (digamma(nu^(-2)) - log(nu^(-2)) - 1)

  dlogC - w / (nu^2 * sigma) + 2 * z / nu^3 - w * z / (nu^2 * sigma)
}
