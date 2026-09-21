#' Local influence diagnostics for the LGG model
#'
#' Computes local influence measures for a fitted LGG model
#' under different perturbation schemes.
#'
#' @param fit An object of class \code{"gamlss"} fitted with the
#'   \code{LGG} family.
#' @param perturbation Character string specifying the perturbation
#'   scheme. Currently available options are \code{"case"} and
#'   \code{"response"}.
#' @param penalized Logical. If \code{TRUE}, uses the penalized
#'   observed Fisher information matrix. Default is \code{TRUE}.
#' @param tol Numerical tolerance used when checking the
#'   singularity of the observed Fisher information matrix.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Delta}: the perturbation matrix;
#'   \item \code{F}: the curvature matrix;
#'   \item \code{eigenvalues}: eigenvalues of \code{F};
#'   \item \code{eigenvectors}: eigenvectors of \code{F};
#'   \item \code{Cdmax}: maximum normal curvature;
#'   \item \code{dmax}: direction associated with the maximum
#'         curvature;
#'   \item \code{Br}: conformal curvature for each observation;
#'   \item \code{limit}: cutoff defined as
#'         \code{mean(Br) + 3 * sd(Br)};
#'   \item \code{influential}: logical vector indicating observations
#'         whose conformal curvature exceeds the cutoff.
#' }
#'
#' @export
localInfluenceLGG <- function(
    fit,
    perturbation = c("case", "response"),
    penalized = TRUE,
    tol = 1e-8
) {

  # ============================================================
  # 1. Verificações
  # ============================================================

  if (!inherits(fit, "gamlss")) {
    stop("fit deve ser um objeto da classe 'gamlss'.")
  }

  if (fit$family[1] != "LGG") {
    stop("A função foi desenvolvida para a família LGG.")
  }

  perturbation <- match.arg(perturbation)


  # ============================================================
  # 2. Informação de Fisher
  # ============================================================

  fisher <- .observed_fisherLGG(
    fit = fit,
    penalized = penalized,
    tol = tol
  )


  # ============================================================
  # 3. Matriz Delta
  # ============================================================

  if (perturbation == "case") {

    Delta <- .delta_caseLGG(fit)

  } else if (perturbation == "response") {

    Delta <- .delta_responseLGG(fit)
  }


  # ============================================================
  # 4. Matriz de curvatura
  # ============================================================

  if (penalized) {

    InFisher <- fisher$I_obs_pen_inv

  } else {

    InFisher <- fisher$I_obs_inv
  }

  F <- t(Delta) %*% InFisher %*% Delta


  # ============================================================
  # 5. Decomposição espectral
  # ============================================================

  deF <- eigen(
    F,
    symmetric = TRUE
  )

  eigenvalues <- deF$values
  eigenvectors <- deF$vectors


  # Maior autovalor

  Cdmax <- eigenvalues[1]

  dmax <- eigenvectors[, 1]

  # Normalização

  if (Cdmax > 0) {
    dmax <- dmax / sqrt(Cdmax)
  }

  dmax <- abs(dmax)


  # ============================================================
  # 6. Curvatura normal conformal
  # ============================================================

  denominator <- sqrt(
    sum(eigenvalues^2)
  )

  if (denominator > 0) {

    Br <- abs(
      diag(F) / denominator
    )

  } else {

    Br <- rep(0, nrow(F))
  }


  # ============================================================
  # 7. Ponto de corte
  # ============================================================

  limit <- mean(Br) + 3 * sd(Br)

  influential <- Br > limit


  # ============================================================
  # 8. Retorno
  # ============================================================

  return(
    list(
      perturbation = perturbation,
      penalized = penalized,

      Delta = Delta,
      F = F,

      eigenvalues = eigenvalues,
      eigenvectors = eigenvectors,

      Cdmax = Cdmax,
      dmax = dmax,

      Br = Br,
      limit = limit,
      influential = influential,

      fisher = fisher
    )
  )
}
