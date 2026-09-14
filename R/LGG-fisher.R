#' Penalized observed Fisher information for the LGG model
#'
#' Computes the observed Fisher information matrix and its penalized
#' version for a fitted LGG GAMLSS model.
#'
#' The current implementation assumes the links
#' \code{identity} for \code{mu}, \code{log} for \code{sigma},
#' and \code{identity} for \code{nu}.
#'
#' Smooth terms are extracted from the GAMLSS object through
#' \code{fit$mu.coefSmo[[k]]$smooth[[1]]$S}.
#'
#' @param fit A fitted \code{gamlss} model with family \code{LGG}.
#' @param data Data frame used to fit the model.
#' @param penalized Logical. If \code{TRUE}, returns the penalized
#'   observed Fisher information.
#'
#' @return A list containing:
#' \describe{
#'   \item{I_obs}{Unpenalized observed Fisher information matrix.}
#'   \item{P}{Penalty matrix.}
#'   \item{I_obs_pen}{Penalized observed Fisher information matrix.}
#'   \item{S}{List of penalty matrices for the smooth terms.}
#'   \item{lambda}{Smoothing parameters.}
#'   \item{X_mu}{Complete design matrix for the mu predictor.}
#'   \item{H_beta}{List of observation-specific Hessian contributions.}
#' }
#'
#' @importFrom stats fitted
#' @importFrom mgcv PredictMat
#'
#' @export

observed_fisher_pen_LGG <- function(fit, data, penalized = TRUE) {

  # ============================================================
  # 1. Verificacoes
  # ============================================================

  if (!inherits(fit, "gamlss")) {
    stop("'fit' deve ser um objeto ajustado pela funcao gamlss().")
  }

  if (missing(data)) {
    stop("'data' deve ser fornecido.")
  }

  # ============================================================
  # 2. Verificar os links
  # ============================================================

  if (fit$mu.link != "identity") {
    stop("A funcao assume link identity para mu.")
  }

  if (fit$sigma.link != "log") {
    stop("A funcao assume link log para sigma.")
  }

  if (fit$nu.link != "identity") {
    stop("A funcao assume link identity para nu.")
  }

  # ============================================================
  # 3. Valores ajustados
  # ============================================================

  mu_hat <- fitted(fit, what = "mu")
  sigma_hat <- fitted(fit, what = "sigma")
  nu_hat <- fitted(fit, what = "nu")

  n <- length(mu_hat)

  # ============================================================
  # 4. Matriz de projeto da parte parametrica
  #
  # ============================================================

  X_all <- as.matrix(fit$mu.x)[,1:(length(fit$mu.coefficients)
                                       -length(fit$mu.coefSm))]

  beta_mu <- coef(fit, what = "mu")

  beta_names <- names(beta_mu)

  # ============================================================
  # 5. Identificar os coeficientes dos termos suaves
  # ============================================================

  smooth_coef_names <- character(0)

  if (!is.null(fit$mu.coefSmo) &&
      length(fit$mu.coefSmo) > 0) {

    for (k in seq_along(fit$mu.coefSmo)) {

      sm <- fit$mu.coefSmo[[k]]
      smo <- sm$smooth[[1]]

      X_k <- mgcv::PredictMat(
        smo,
        data
      )

      K_k <- ncol(X_k)

      label_k <- smo$label

      if (is.null(label_k) || is.na(label_k)) {
        label_k <- paste0("s", k)
      }

      smooth_coef_names <- c(
        smooth_coef_names,
        paste0(label_k, seq_len(K_k))
      )
    }
  }

  # ============================================================
  # 6. Identificar somente os coeficientes parametricos
  # ============================================================

  param_names <- beta_names[
    !beta_names %in% smooth_coef_names
  ]

  # Remover possiveis colunas auxiliares ga(~s(...))
  param_names <- param_names[
    param_names %in% colnames(X_all)
  ]

  if (length(param_names) == 0) {
    stop(
      "Nao foi possivel identificar os termos parametricos ",
      "da formula de mu."
    )
  }

  X_param <- X_all[
    ,
    param_names,
    drop = FALSE
  ]

  # ============================================================
  # 7. Construir a matriz completa X_mu
  # ============================================================

  X_mu <- X_param

  smooth_info <- list()

  if (!is.null(fit$mu.coefSmo) &&
      length(fit$mu.coefSmo) > 0) {

    for (k in seq_along(fit$mu.coefSmo)) {

      sm <- fit$mu.coefSmo[[k]]
      smo <- sm$smooth[[1]]

      # Matriz de projeto do termo suave
      X_k <- mgcv::PredictMat(
        smo,
        data
      )

      # ========================================================
      # MATRIZES DE PENALIZACAO
      #
      # Estrutura utilizada:
      #
      # fit$mu.coefSmo[[k]]$smooth[[1]]$S
      # ========================================================

      S_k <- smo$S

      # Parametros de suavizacao
      lambda_k <- sm$sp

      # Verificacao
      if (length(S_k) != length(lambda_k)) {

        if (length(S_k) == 1 &&
            length(lambda_k) >= 1) {

          lambda_k <- lambda_k[1]

        } else {

          stop(
            "O numero de matrizes de penalizacao e ",
            "parametros de suavizacao nao coincide ",
            "para o termo suave ", k, "."
          )
        }
      }

      # Adicionar matriz do termo suave
      X_mu <- cbind(
        X_mu,
        X_k
      )

      # Guardar informacoes
      smooth_info[[k]] <- list(
        X = X_k,
        S = S_k,
        lambda = lambda_k,
        label = smo$label
      )
    }
  }

  # ============================================================
  # 8. Dimensao do vetor de parametros
  # ============================================================

  p_mu <- ncol(X_mu)

  # Um intercepto para sigma e um para nu
  p <- p_mu + 2

  ind_mu <- seq_len(p_mu)

  ind_sigma <- p_mu + 1

  ind_nu <- p_mu + 2

  # ============================================================
  # 9. Inicializar a informacao de Fisher observada
  # ============================================================

  I_obs <- matrix(
    0,
    nrow = p,
    ncol = p
  )

  # ============================================================
  # 10. Guardar as Hessianas individuais
  # ============================================================

  H_beta <- vector(
    "list",
    n
  )

  # ============================================================
  # 11. Calculo da informacao observada
  # ============================================================

  for (i in seq_len(n)) {

    y_i <- fit$y[i]

    mu_i <- mu_hat[i]

    sigma_i <- sigma_hat[i]

    nu_i <- nu_hat[i]

    # ----------------------------------------------------------
    # Hessiana em relacao a (mu, sigma, nu)
    # ----------------------------------------------------------

    H_theta <- matrix(
      c(
        d2ldm2(
          y_i,
          mu_i,
          sigma_i,
          nu_i
        ),

        d2ldmdd(
          y_i,
          mu_i,
          sigma_i,
          nu_i
        ),

        d2ldmdv(
          y_i,
          mu_i,
          sigma_i,
          nu_i
        ),

        d2ldmdd(
          y_i,
          mu_i,
          sigma_i,
          nu_i
        ),

        d2ldd2(
          y_i,
          mu_i,
          sigma_i,
          nu_i
        ),

        d2ldddv(
          y_i,
          mu_i,
          sigma_i,
          nu_i
        ),

        d2ldmdv(
          y_i,
          mu_i,
          sigma_i,
          nu_i
        ),

        d2ldddv(
          y_i,
          mu_i,
          sigma_i,
          nu_i
        ),

        d2ldv2(
          y_i,
          mu_i,
          sigma_i,
          nu_i
        )
      ),
      nrow = 3,
      ncol = 3,
      byrow = TRUE
    )

    # ----------------------------------------------------------
    # Matriz de derivadas dos parametros em relacao aos
    # preditores lineares
    #
    # mu    = eta_mu
    # sigma = exp(eta_sigma)
    # nu    = eta_nu
    # ----------------------------------------------------------

    A_i <- matrix(
      0,
      nrow = 3,
      ncol = p
    )

    # mu
    A_i[
      1,
      ind_mu
    ] <- X_mu[i, ]

    # sigma
    A_i[
      2,
      ind_sigma
    ] <- sigma_i

    # nu
    A_i[
      3,
      ind_nu
    ] <- 1

    # ----------------------------------------------------------
    # Hessiana em relacao aos coeficientes
    # ----------------------------------------------------------

    H_beta_i <- t(A_i) %*%
      H_theta %*%
      A_i

    # ----------------------------------------------------------
    # Termo adicional devido ao link log de sigma
    # ----------------------------------------------------------

    s_sigma <- score_sigma(
      y_i,
      mu_i,
      sigma_i,
      nu_i
    )

    H_beta_i[
      ind_sigma,
      ind_sigma
    ] <-
      H_beta_i[
        ind_sigma,
        ind_sigma
      ] +
      s_sigma * sigma_i

    # ----------------------------------------------------------
    # Armazenar Hessiana individual
    # ----------------------------------------------------------

    H_beta[[i]] <- H_beta_i

    # ----------------------------------------------------------
    # Informacao observada:
    #
    # I_obs = - sum(H_i)
    # ----------------------------------------------------------

    I_obs <- I_obs - H_beta_i
  }

  # ============================================================
  # 12. Construir a matriz de penalizacao P
  # ============================================================

  P <- matrix(
    0,
    nrow = p,
    ncol = p
  )

  S_list <- list()

  lambda_list <- list()

  # Primeiro coeficiente dos termos suaves
  start_smooth <- ncol(X_param) + 1

  # ============================================================
  # 13. Inserir as penalizacoes dos termos suaves
  # ============================================================

  if (length(smooth_info) > 0) {

    for (k in seq_along(smooth_info)) {

      info_k <- smooth_info[[k]]

      X_k <- info_k$X

      S_k <- info_k$S

      lambda_k <- info_k$lambda

      K_k <- ncol(X_k)

      # Indices dos coeficientes do termo suave
      ind_k <- start_smooth:
        (start_smooth + K_k - 1)

      # --------------------------------------------------------
      # Caso com uma matriz de penalizacao
      # --------------------------------------------------------

      if (length(S_k) == 1) {

        P[
          ind_k,
          ind_k
        ] <-
          P[
            ind_k,
            ind_k
          ] +
          lambda_k[1] * S_k[[1]]

      } else {

        # ------------------------------------------------------
        # Caso com mais de uma matriz de penalizacao
        # ------------------------------------------------------

        if (length(lambda_k) != length(S_k)) {

          stop(
            "Para o termo suave ", k,
            ", o numero de parametros de suavizacao ",
            "nao coincide com o numero de matrizes S."
          )
        }

        P_k <- matrix(
          0,
          nrow = K_k,
          ncol = K_k
        )

        for (j in seq_along(S_k)) {

          P_k <- P_k +
            lambda_k[j] * S_k[[j]]
        }

        P[
          ind_k,
          ind_k
        ] <-
          P[
            ind_k,
            ind_k
          ] +
          P_k
      }

      # Guardar
      S_list[[k]] <- S_k

      lambda_list[[k]] <- lambda_k

      # Proximo termo suave
      start_smooth <-
        start_smooth + K_k
    }
  }

  # ============================================================
  # 14. Informacao de Fisher observada penalizada
  # ============================================================

  if (penalized) {

    I_obs_pen <-
      I_obs + P

  } else {

    I_obs_pen <-
      I_obs
  }

  # ============================================================
  # 15. Nomes dos coeficientes
  # ============================================================

  smooth_names <- character(0)

  if (length(smooth_info) > 0) {

    for (k in seq_along(smooth_info)) {

      K_k <- ncol(
        smooth_info[[k]]$X
      )

      label_k <-
        smooth_info[[k]]$label

      if (is.null(label_k) ||
          is.na(label_k)) {

        label_k <- paste0(
          "s",
          k
        )
      }

      smooth_names <- c(
        smooth_names,
        paste0(
          label_k,
          seq_len(K_k)
        )
      )
    }
  }

  coef_names <- c(
    colnames(X_param),
    smooth_names,
    "sigma.(Intercept)",
    "nu.(Intercept)"
  )

  # ============================================================
  # 16. Atribuir nomes
  # ============================================================

  dimnames(I_obs) <- list(
    coef_names,
    coef_names
  )

  dimnames(P) <- list(
    coef_names,
    coef_names
  )

  dimnames(I_obs_pen) <- list(
    coef_names,
    coef_names
  )

  # ============================================================
  # 17. Retorno
  # ============================================================

  return(
    list(
      I_obs = I_obs,
      P = P,
      I_obs_pen = I_obs_pen,
      S = S_list,
      lambda = lambda_list,
      X_mu = X_mu,
      X_param = X_param,
      H_beta = H_beta,
      smooth_info = smooth_info
    )
  )
}
