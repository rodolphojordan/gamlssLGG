#' Approximate second derivative with respect to mu
#'
#' Computes the approximate second derivative of the LGG
#' log-density with respect to the location parameter.
#'
#' @param y Vector of observations.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#' @param nu Shape parameter.
#'
#' @return Numeric vector.
#'
#' @export
d2ldm2 <- function(y, mu, sigma, nu) {

  -(score_mu(y, mu, sigma, nu)^2)

}


#' Approximate second derivative with respect to sigma
#'
#' @param y Vector of observations.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#' @param nu Shape parameter.
#'
#' @return Numeric vector.
#'
#' @export
d2ldd2 <- function(y, mu, sigma, nu) {

  -(score_sigma(y, mu, sigma, nu)^2)

}


#' Approximate second derivative with respect to nu
#'
#' @param y Vector of observations.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#' @param nu Shape parameter.
#'
#' @return Numeric vector.
#'
#' @export
d2ldv2 <- function(y, mu, sigma, nu) {

  -(score_nu(y, mu, sigma, nu)^2)

}


#' Approximate mixed derivative between mu and sigma
#'
#' @export
d2ldmdd <- function(y, mu, sigma, nu) {

  rep(0, length(y))

}


#' Approximate mixed derivative between mu and nu
#'
#' @export
d2ldmdv <- function(y, mu, sigma, nu) {

  rep(0, length(y))

}


#' Approximate mixed derivative between sigma and nu
#'
#' @export
d2ldddv <- function(y, mu, sigma, nu) {

  rep(0, length(y))

}
