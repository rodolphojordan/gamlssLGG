# Tests for LGG package
context("LGG Distribution Tests")

test_that("dLGG integrates to 1", {
  # Test numeric integration
  result <- integrate(
    function(y) dLGG(y, mu = 0, sigma = 1, nu = 0.5),
    -Inf,
    Inf
  )
  expect_equal(result$value, 1, tolerance = 1e-6)
})

test_that("pLGG and qLGG are inverses", {
  p <- seq(0.01, 0.99, by = 0.01)
  q <- qLGG(p, mu = 0, sigma = 1, nu = 0.5)
  p2 <- pLGG(q, mu = 0, sigma = 1, nu = 0.5)
  expect_equal(p, p2, tolerance = 1e-6)
})

test_that("dLGG limit to normal when nu -> 0", {
  y <- c(-2, -1, 0, 1, 2)
  d1 <- dLGG(y, mu = 0, sigma = 1, nu = 1e-7)
  d2 <- dnorm(y, mean = 0, sd = 1)
  expect_equal(d1, d2, tolerance = 1e-6)
})

test_that("gamlss family works", {
  library(gamlss)
  set.seed(123)
  n <- 50
  x <- runif(n)
  mu <- 1 + 2*x
  y <- rLGG(n, mu = mu, sigma = 0.5, nu = 0.4)
  dados <- data.frame(y, x)

  fit <- gamlss(y ~ x, family = LGG(), data = dados,
                control = gamlss.control(n.cyc = 5, trace = FALSE))

  expect_true(inherits(fit, "gamlss"))
})
