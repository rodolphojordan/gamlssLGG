.delta_responseLGG <- function(fit) {

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
  # 2. Resposta e valores ajustados
  # ============================================================

  y <- model.response(model.frame(fit))

  mu_hat <- fitted(fit, "mu")
  sigma_hat <- fitted(fit, "sigma")[1]
  nu_hat <- fitted(fit, "nu")[1]

  # Resíduo padronizado
  epsilon <- (y - mu_hat) / sigma_hat


  # ============================================================
  # 3. Matriz X dos efeitos paramétricos de mu
  # ============================================================

  X_mu_full <- fit$mu.x

  beta_mu_all <- coef(fit, "mu")

  ind_param_mu <- !is.na(beta_mu_all)

  X_beta <- X_mu_full[
    ,
    ind_param_mu,
    drop = FALSE
  ]


  # ============================================================
  # 4. Matriz N dos coeficientes dos termos suaves
  # ============================================================

  terSmo <- getSmo(fit)

  if (length(terSmo$smooth) > 0) {

    X_smooth_list <- lapply(
      terSmo$smooth,
      function(sm) {

        model.matrix(terSmo)[
          ,
          sm$first.para:sm$last.para,
          drop = FALSE
        ]
      }
    )

    N <- do.call(cbind, X_smooth_list)

  } else {

    N <- matrix(
      nrow = nrow(X_beta),
      ncol = 0
    )
  }


  # ============================================================
  # 5. Delta_beta
  # ============================================================

  Delta_beta <- sweep(
    X_beta,
    1,
    exp(epsilon) / (nu_hat * sigma_hat),
    "*"
  )


  # ============================================================
  # 6. Delta_gamma
  # ============================================================

  if (ncol(N) > 0) {

    Delta_gamma <- sweep(
      N,
      1,
      exp(nu_hat * epsilon) / sigma_hat^2,
      "*"
    )

  } else {

    Delta_gamma <- matrix(
      nrow = nrow(X_beta),
      ncol = 0
    )
  }


  # ============================================================
  # 7. Delta_sigma
  # ============================================================

  Delta_sigma <-
    (1 / (nu_hat * sigma_hat^2)) *
    (
      nu_hat * epsilon * exp(nu_hat * epsilon) +
        exp(nu_hat * epsilon) -
        1
    )


  # ============================================================
  # 8. Delta_nu
  # ============================================================

  Delta_nu <-
    -1 / (nu_hat^2 * sigma_hat) +
    exp(nu_hat * epsilon) /
    (nu_hat^2 * sigma_hat) -
    epsilon * exp(nu_hat * epsilon) /
    (nu_hat * sigma_hat)


  # ============================================================
  # 9. Matriz Delta
  # ============================================================

  Delta_y <- rbind(
    t(Delta_beta),
    t(Delta_gamma),
    t(Delta_sigma),
    t(Delta_nu)
  )


  # ============================================================
  # 10. Nomes das colunas
  # ============================================================

  names_beta <- colnames(X_beta)

  names_gamma <- colnames(N)

  names_sigma <- paste0(
    "sigma.",
    names(coef(fit, "sigma"))
  )

  names_nu <- paste0(
    "nu.",
    names(coef(fit, "nu"))
  )

  rownames(Delta_y) <- c(
    names_beta,
    names_gamma,
    names_sigma,
    names_nu
  )


  # ============================================================
  # 11. Retorno
  # ============================================================

  return(Delta_y)
}
