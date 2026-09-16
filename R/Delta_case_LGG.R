#' Delta matrix for case-weight perturbation in the LGG model
#'
#' Computes the Delta matrix associated with case-weight perturbation
#' for a fitted LGG GAMLSS model.
#'
#' The matrix is constructed from the derivatives of the individual
#' log-likelihood contributions with respect to the complete vector of
#' regression coefficients. It is partitioned into the parametric
#' coefficients of mu, the coefficients of all smooth terms in mu,
#' and the parameters of sigma and nu.
#'
#' The current implementation assumes the links
#' \code{identity} for \code{mu}, \code{log} for \code{sigma},
#' and \code{identity} for \code{nu}.
#'
#' For case-weight perturbation, the resulting matrix has the form
#'
#' \deqn{
#' \Delta =
#' [\Delta_\beta,\Delta_\gamma,\Delta_\sigma,\Delta_\nu],
#' }
#'
#' where \code{beta} denotes the parametric coefficients of the
#' \code{mu} predictor and \code{gamma} denotes the coefficients of
#' all smooth terms in the \code{mu} predictor.
#'
#' @param fit A fitted \code{gamlss} model with family \code{LGG}.
#'
#' @return A matrix containing the Delta matrix for case-weight
#'   perturbation. The columns are ordered as follows:
#'   \code{mu} parametric coefficients, smooth-term coefficients,
#'   \code{sigma} coefficients, and \code{nu} coefficients.
#'
#' @importFrom stats fitted
#'
#' @export

Delta_case_LGG <- function(fit) {

  # ============================================================
  # 1. Verificações
  # ============================================================

  if (!inherits(fit, "gamlss")) {
    stop("fit deve ser um objeto da classe 'gamlss'.")
  }

  fam <- fit$family[1]

  if (fam != "LGG") {
    stop("A função foi desenvolvida para a família LGG.")
  }


  # ============================================================
  # 2. Resposta
  # ============================================================

  y <- model.response(model.frame(fit))

  n <- length(y)


  # ============================================================
  # 3. Valores ajustados
  # ============================================================

  mu_hat <- fitted(fit, what = "mu")

  # Mantidos como escalares, conforme a formulação atual.
  # Isso é apropriado quando sigma e nu possuem apenas intercepto.
  sigma_hat <- fitted(fit, what = "sigma")[1]

  nu_hat <- fitted(fit, what = "nu")[1]


  # ============================================================
  # 4. Matriz paramétrica de mu
  # ============================================================

  X_mu_full <- fit$mu.x

  beta_mu_all <- coef(fit, "mu")

  # Remove as colunas auxiliares associadas aos termos suaves
  ind_param_mu <- !is.na(beta_mu_all)

  X_mu_param <- X_mu_full[
    ,
    ind_param_mu,
    drop = FALSE
  ]


  # ============================================================
  # 5. Matriz dos termos suaves
  # ============================================================

  terSmo <- getSmo(fit)

  n_smooth <- length(terSmo$smooth)

  X_smooth_list <- vector(
    "list",
    n_smooth
  )

  if (n_smooth > 0) {

    X_model_smooth <- model.matrix(terSmo)

    for (j in seq_len(n_smooth)) {

      sm <- terSmo$smooth[[j]]

      X_smooth_list[[j]] <- X_model_smooth[
        ,
        sm$first.para:sm$last.para,
        drop = FALSE
      ]
    }

    # Junta todos os termos suaves
    N <- do.call(
      cbind,
      X_smooth_list
    )

  } else {

    N <- NULL
  }


  # ============================================================
  # 6. Índices dos parâmetros
  # ============================================================

  p_beta <- ncol(X_mu_param)

  p_gamma <- if (is.null(N)) {
    0
  } else {
    ncol(N)
  }

  idx_beta <- seq_len(p_beta)

  if (p_gamma > 0) {

    idx_gamma <- p_beta + seq_len(p_gamma)

  } else {

    idx_gamma <- integer(0)
  }


  # ============================================================
  # 7. Resíduo padronizado
  # ============================================================

  epsilon <- (y - mu_hat) / sigma_hat


  # ============================================================
  # 8. zeta_nu
  # ============================================================

  z_nu <- (1 / nu_hat) +
    (2 / nu_hat^3) *
    (
      digamma(1 / nu_hat^2) +
        2 * log(abs(nu_hat)) -
        1
    )


  # ============================================================
  # 9. Fator comum aos blocos beta e gamma
  # ============================================================

  fator <- (1 / (nu_hat * sigma_hat)) *
    (
      -1 +
        exp(nu_hat * epsilon)
    )


  # ============================================================
  # 10. Delta_beta
  # ============================================================

  Delta_beta <- sweep(
    X_mu_param,
    1,
    fator,
    "*"
  )


  # ============================================================
  # 11. Delta_gamma
  # ============================================================

  if (p_gamma > 0) {

    Delta_gamma <- sweep(
      N,
      1,
      fator,
      "*"
    )

  } else {

    Delta_gamma <- matrix(
      numeric(0),
      nrow = n,
      ncol = 0
    )
  }


  # ============================================================
  # 12. Delta_sigma
  # ============================================================

  Delta_sigma <- -1 / sigma_hat -
    epsilon / (nu_hat * sigma_hat) +
    (
      epsilon /
        (nu_hat * sigma_hat)
    ) *
    exp(nu_hat * epsilon)


  # ============================================================
  # 13. Delta_nu
  # ============================================================

  Delta_nu <- z_nu -
    epsilon / nu_hat^2 +
    (2 / nu_hat^3) *
    exp(nu_hat * epsilon) -
    (1 / nu_hat^2) *
    epsilon *
    exp(nu_hat * epsilon)


  # ============================================================
  # 14. Matriz Delta completa
  # ============================================================

  Delta <- cbind(
    Delta_beta,
    Delta_gamma,
    Delta_sigma,
    Delta_nu
  )


  # ============================================================
  # 15. Nomes das colunas
  # ============================================================

  names_beta <- colnames(X_mu_param)

  names_gamma <- if (p_gamma > 0) {
    colnames(N)
  } else {
    character(0)
  }

  names_sigma <- paste0(
    "sigma.",
    names(coef(fit, "sigma"))
  )

  names_nu <- paste0(
    "nu.",
    names(coef(fit, "nu"))
  )

  colnames(Delta) <- c(
    names_beta,
    names_gamma,
    names_sigma,
    names_nu
  )


  # ============================================================
  # 16. Retorno
  # ============================================================

  return(Delta)
}
