#' GAMLSS Family for the Log-Generalized Gamma Distribution
#'
#' @param mu.link Link function for mu
#' @param sigma.link Link function for sigma
#' @param nu.link Link function for nu
#'
#' @return A gamlss.family object
#' @export
#'
#' @import gamlss
#' @import gamlss.dist
#'
#' @examples
#' \dontrun{
#' library(gamlss)
#' set.seed(123)
#' n <- 300
#' x <- runif(n)
#' mu <- 1 + 2*x
#' y <- rLGG(n, mu = mu, sigma = 0.5, nu = 0.4)
#' dados <- data.frame(y, x)
#' fit <- gamlss(y ~ x, family = LGG(), data = dados, trace = TRUE)
#' summary(fit)
#' }
LGG <- function(mu.link = "identity",
                sigma.link = "log",
                nu.link = "identity") {

  # Links dos parâmetros
  mstats <- checklink(
    "mu.link",
    "LGG",
    substitute(mu.link),
    c("identity", "log", "inverse", "own")
  )

  dstats <- checklink(
    "sigma.link",
    "LGG",
    substitute(sigma.link),
    c("log", "identity", "inverse", "own")
  )

  vstats <- checklink(
    "nu.link",
    "LGG",
    substitute(nu.link),
    c("identity", "log", "inverse", "own")
  )

  structure(
    list(
      family = c("LGG", "Log-Generalized Gamma"),
      parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE),
      nopar = 3,
      type = "Continuous",

      mu.link = as.character(substitute(mu.link)),
      sigma.link = as.character(substitute(sigma.link)),
      nu.link = as.character(substitute(nu.link)),

      mu.linkfun = mstats$linkfun,
      sigma.linkfun = dstats$linkfun,
      nu.linkfun = vstats$linkfun,

      mu.linkinv = mstats$linkinv,
      sigma.linkinv = dstats$linkinv,
      nu.linkinv = vstats$linkinv,

      mu.dr = mstats$mu.eta,
      sigma.dr = dstats$mu.eta,
      nu.dr = vstats$mu.eta,

      # Scores analíticos
      dldm = function(y, mu, sigma, nu) {
        score_mu(y, mu, sigma, nu)
      },

      dldd = function(y, mu, sigma, nu) {
        score_sigma(y, mu, sigma, nu)
      },

      dldv = function(y, mu, sigma, nu) {
        score_nu(y, mu, sigma, nu)
      },

      # Hessiana aproximada
      d2ldm2 = function(y, mu, sigma, nu) {
        -(score_mu(y, mu, sigma, nu)^2)
      },

      d2ldd2 = function(y, mu, sigma, nu) {
        -(score_sigma(y, mu, sigma, nu)^2)
      },

      d2ldv2 = function(y, mu, sigma, nu) {
        -(score_nu(y, mu, sigma, nu)^2)
      },

      d2ldmdd = function(y, mu, sigma, nu) {
        rep(0, length(y))
      },

      d2ldmdv = function(y, mu, sigma, nu) {
        rep(0, length(y))
      },

      d2ldddv = function(y, mu, sigma, nu) {
        rep(0, length(y))
      },

      # Deviance
      G.dev.incr = function(y, mu, sigma, nu, ...) {
        -2 * dLGG(y, mu, sigma, nu, log = TRUE)
      },

      # Resíduos quantílicos
      rqres = expression(
        rqres(
          pfun = "pLGG",
          type = "Continuous",
          y = y,
          mu = mu,
          sigma = sigma,
          nu = nu
        )
      ),

      # Valores iniciais
      mu.initial = expression(
        mu <- rep(mean(y), length(y))
      ),

      sigma.initial = expression(
        sigma <- rep(sd(y), length(y))
      ),

      nu.initial = expression(
        nu <- rep(0.5, length(y))
      ),

      # Restrições
      mu.valid = function(mu) TRUE,

      sigma.valid = function(sigma) all(sigma > 0),

      nu.valid = function(nu) TRUE,

      y.valid = function(y) all(is.finite(y))
    ),
    class = c("gamlss.family", "family")
  )
}
