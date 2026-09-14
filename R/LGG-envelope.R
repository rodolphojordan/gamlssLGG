envelope <- function(modelo) {

  dados <- na.omit(eval(modelo$call$data))
  nsim <- 100
  n <- dim(dados)[1]

  r1 <- sort(residuals(modelo))
  m1 <- matrix(0, nrow = n, ncol = nsim)

  mu <- fitted(modelo, "mu")
  sigma <- fitted(modelo, "sigma")
  nu <- fitted(modelo, "nu")

  a2 <- replicate(
    nsim,
    rLGG(
      n = n,
      mu = mu,
      sigma = sigma,
      nu = nu
    )
  )

  for (i in 1:nsim) {
    dados$y <- a2[, i]

    aj <- update(
      modelo,
      y ~ .,
      data = dados
    )

    m1[, i] <- sort(residuals(aj))
  }

  li <- apply(m1, 1, quantile, 0.025)
  m  <- apply(m1, 1, quantile, 0.5)
  ls <- apply(m1, 1, quantile, 0.975)

  quantis <- qnorm((1:n - 0.5) / n)

  plot(
    rep(quantis, 2),
    c(li, ls),
    type = "n",
    ann = FALSE
  )

  lines(quantis, li)
  lines(quantis, m, lty = 2)
  lines(quantis, ls)
  points(quantis, r1, pch = 20, cex = 0.75)
}

