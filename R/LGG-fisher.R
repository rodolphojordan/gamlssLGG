.observed_fisherLGG <- function(fit,
                                    penalized = TRUE,
                                    tol = 1e-8) {

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
  # 2. Valores ajustados
  # ============================================================

  y <- model.response(model.frame(fit))

  mu <- fitted(fit, what = "mu")
  sigma <- fitted(fit, what = "sigma")
  nu <- fitted(fit, what = "nu")

  n <- length(y)


  # ============================================================
  # 3. Estrutura dos suavizadores
  # ============================================================

  terSmo <- getSmo(fit)

  d <- length(terSmo$sp)


  # ============================================================
  # 4. Matriz paramétrica de mu
  # ============================================================

  X_mu_full <- fit$mu.x

  # Mantém somente as colunas paramétricas.

  beta_mu_all <- coef(fit, "mu")

  ind_param_mu <- !is.na(beta_mu_all)

  X_mu_param <- X_mu_full[, ind_param_mu, drop = FALSE]

  beta_mu_param <- beta_mu_all[ind_param_mu]


  # ============================================================
  # 5. Matrizes dos suavizadores
  # ============================================================

  X_mu_smooth <- NULL

  smooth_info <- vector("list", d)

  if (d > 0) {

    for (j in seq_len(d)) {

      sm <- terSmo$smooth[[j]]

      # Matriz de desenho da spline
      X_smooth_j <- model.matrix(terSmo)[,
                                         sm$first.para:sm$last.para,
                                         drop = FALSE]

      # Número de coeficientes
      q_j <- ncol(X_smooth_j)

      X_mu_smooth <- cbind(
        X_mu_smooth,
        X_smooth_j
      )

      # Penalização
      S_j <- sm$S

      # Parâmetro de suavização
      lambda_j <- terSmo$sp[j]

      smooth_info[[j]] <- list(
        X = X_smooth_j,
        S = S_j,
        lambda = lambda_j,
        q = q_j
      )
    }
  }


  # ============================================================
  # 6. Matriz completa do preditor de mu
  # ============================================================

  if (!is.null(X_mu_smooth)) {

    X_mu <- cbind(
      X_mu_param,
      X_mu_smooth
    )

  } else {

    X_mu <- X_mu_param
  }


  # ============================================================
  # 7. Matrizes de sigma e nu
  # ============================================================

  X_sigma <- fit$sigma.x
  X_nu <- fit$nu.x


  # ============================================================
  # 8. Número total de parâmetros
  # ============================================================

  p_mu <- ncol(X_mu)
  p_sigma <- ncol(X_sigma)
  p_nu <- ncol(X_nu)

  p <- p_mu + p_sigma + p_nu


  # ============================================================
  # 9. Índices dos parâmetros
  # ============================================================

  ind_mu <- seq_len(p_mu)

  ind_sigma <- p_mu + seq_len(p_sigma)

  ind_nu <- p_mu + p_sigma + seq_len(p_nu)


  # ============================================================
  # 10. Informação observada
  # ============================================================

  I_obs <- matrix(
    0,
    nrow = p,
    ncol = p
  )


  # ============================================================
  # 11. Hessiana individual
  #
  # Aqui usamos derivadas numéricas da log-verossimilhança
  # individual para obter também os termos cruzados.
  # ============================================================

  loglik_i <- function(theta, yi) {

    dLGG(
      y = yi,
      mu = theta[1],
      sigma = theta[2],
      nu = theta[3],
      log = TRUE
    )
  }


  # ============================================================
  # 12. Loop sobre as observações
  # ============================================================

  H_beta_list <- vector("list", n)

  for (i in seq_len(n)) {

    theta_i <- c(
      mu[i],
      sigma[i],
      nu[i]
    )


    # ----------------------------------------------------------
    # Hessiana da log-verossimilhança em relação a
    # (mu, sigma, nu)
    # ----------------------------------------------------------

    H_theta <- optimHess(
      par = theta_i,
      fn = loglik_i,
      yi = y[i]
    )


    # ----------------------------------------------------------
    # Matriz Jacobiana da transformação
    #
    # theta = (mu, sigma, nu)
    #
    # em relação ao vetor completo de coeficientes
    # ----------------------------------------------------------

    A_i <- matrix(
      0,
      nrow = 3,
      ncol = p
    )


    # mu
    A_i[1, ind_mu] <- X_mu[i, ]


    # sigma
    #
    # sigma = exp(eta_sigma)
    #
    # d sigma / d beta_sigma = sigma * X_sigma
    A_i[2, ind_sigma] <- sigma[i] * X_sigma[i, ]


    # nu
    A_i[3, ind_nu] <- X_nu[i, ]


    # ----------------------------------------------------------
    # Termo da regra da cadeia
    # ----------------------------------------------------------

    H_beta_i <- t(A_i) %*% H_theta %*% A_i


    # ----------------------------------------------------------
    # Correção da segunda derivada do link de sigma
    #
    # sigma = exp(eta_sigma)
    #
    # d2 sigma / d eta_sigma^2 = sigma
    # ----------------------------------------------------------

    score_sigma_i <- score_sigma(
      y = y[i],
      mu = mu[i],
      sigma = sigma[i],
      nu = nu[i]
    )

    H_beta_i[
      ind_sigma,
      ind_sigma
    ] <-
      H_beta_i[
        ind_sigma,
        ind_sigma
      ] +
      score_sigma_i *
      sigma[i] *
      tcrossprod(
        X_sigma[i, ],
        X_sigma[i, ]
      )


    # ----------------------------------------------------------
    # Informação observada = - Hessiana
    # ----------------------------------------------------------

    I_obs <- I_obs - H_beta_i

    H_beta_list[[i]] <- H_beta_i
  }


  # ============================================================
  # 13. Matriz de penalização
  # ============================================================

  P <- matrix(
    0,
    nrow = p,
    ncol = p
  )

  if (penalized && d > 0) {

    pos_smooth <- length(beta_mu_param) + 1

    for (j in seq_len(d)) {

      q_j <- smooth_info[[j]]$q

      lambda_j <- smooth_info[[j]]$lambda

      S_j <- smooth_info[[j]]$S

      # Pode existir mais de uma matriz de penalização
      # para um mesmo suavizador.
      if (is.list(S_j)) {

        P_j <- matrix(
          0,
          nrow = q_j,
          ncol = q_j
        )

        for (k in seq_along(S_j)) {

          P_j <- P_j +
            lambda_j[k] * S_j[[k]]
        }

      } else {

        P_j <- lambda_j * S_j
      }


      ind_j <- pos_smooth:(pos_smooth + q_j - 1)

      P[
        ind_j,
        ind_j
      ] <- P[
        ind_j,
        ind_j
      ] + P_j

      pos_smooth <- pos_smooth + q_j
    }
  }


  # ============================================================
  # 14. Informação observada penalizada
  # ============================================================

  I_obs_pen <- I_obs + P


  # ============================================================
  # 15. Função para obter inversa de forma segura
  #
  # Regularização somente se a matriz for singular ou
  # numericamente singular.
  # ============================================================

  safe_inverse <- function(M, tol = 1e-8) {

    M <- (M + t(M)) / 2

    ev <- eigen(
      M,
      symmetric = TRUE,
      only.values = TRUE
    )$values

    scale_M <- max(abs(ev))

    singular <- (
      any(!is.finite(ev)) ||
        scale_M == 0 ||
        min(abs(ev)) <= tol * scale_M
    )

    if (!singular) {

      return(list(
        inverse = solve(M),
        regularized = FALSE,
        lambda_reg = 0,
        eigenvalues = ev
      ))
    }


    # ----------------------------------------------------------
    # Regularização somente quando necessária
    # ----------------------------------------------------------

    lambda_reg <-
      max(tol * scale_M - min(ev), 0)

    lambda_reg <-
      max(lambda_reg, tol * scale_M)

    M_reg <-
      M + lambda_reg * diag(nrow(M))

    list(
      inverse = solve(M_reg),
      regularized = TRUE,
      lambda_reg = lambda_reg,
      eigenvalues = ev,
      matrix_regularized = M_reg
    )
  }


  # ============================================================
  # 16. Inversas
  # ============================================================

  inv_obs <- safe_inverse(
    I_obs,
    tol = tol
  )

  inv_obs_pen <- safe_inverse(
    I_obs_pen,
    tol = tol
  )


  # ============================================================
  # 17. Nomes dos parâmetros
  # ============================================================

  names_mu_param <- names(beta_mu_param)

  names_smooth <- character(0)

  if (d > 0) {

    for (j in seq_len(d)) {

      names_smooth <- c(
        names_smooth,
        colnames(smooth_info[[j]]$X)
      )
    }
  }

  names_sigma <- names(coef(fit, "sigma"))
  names_nu <- names(coef(fit, "nu"))

  parameter_names <- c(
    names_mu_param,
    names_smooth,
    paste0("sigma.", names_sigma),
    paste0("nu.", names_nu)
  )

  rownames(I_obs) <- parameter_names
  colnames(I_obs) <- parameter_names

  rownames(P) <- parameter_names
  colnames(P) <- parameter_names

  rownames(I_obs_pen) <- parameter_names
  colnames(I_obs_pen) <- parameter_names

  rownames(inv_obs$inverse) <- parameter_names
  colnames(inv_obs$inverse) <- parameter_names

  rownames(inv_obs_pen$inverse) <- parameter_names
  colnames(inv_obs_pen$inverse) <- parameter_names


  # ============================================================
  # 18. Retorno
  # ============================================================

  return(list(

    # Informação observada
    I_obs = I_obs,

    # Penalização
    P = P,

    # Informação observada penalizada
    I_obs_pen = I_obs_pen,

    # Inversa da informação observada
    I_obs_inv = inv_obs$inverse,

    # Inversa da informação penalizada
    I_obs_pen_inv = inv_obs_pen$inverse,

    # Informação sobre regularização
    regularized_obs = inv_obs$regularized,
    regularized_pen = inv_obs_pen$regularized,

    lambda_reg_obs = inv_obs$lambda_reg,
    lambda_reg_pen = inv_obs_pen$lambda_reg,

    eigenvalues_obs = inv_obs$eigenvalues,
    eigenvalues_pen = inv_obs_pen$eigenvalues,

    # Matrizes de desenho
    X_mu = X_mu,
    X_mu_param = X_mu_param,
    X_mu_smooth = X_mu_smooth,
    X_sigma = X_sigma,
    X_nu = X_nu,

    # Informação dos suavizadores
    smooth_info = smooth_info,

    # Hessianas individuais
    H_beta = H_beta_list,

    # Nomes
    parameter_names = parameter_names
  ))
}


