
#' Hat Matrix for Penalized Thin Plate Splines
#'
#' Computes the hat matrix for the location component of a
#' GAMLSS model with parametric and multiple smooth terms.
#'
#' @param object A fitted GAMLSS model.
#'
#' @return A list containing the hat matrix and leverage values.
#'
#' @details
#' The hat matrix is defined as
#'
#' H = X (X'X + sigma^2 K)^(-1) X',
#'
#' where K is the block diagonal penalty matrix.
#' Multiple smooth terms are supported, provided that each
#' smooth term has one penalty matrix and one smoothing
#' parameter.
#'
#' @importFrom Matrix bdiag
#' @export

LGGhat <- function(object) {

  # ==========================================================
  # 1. Verify the fitted model
  # ==========================================================

  if (!inherits(object, "gamlss")) {
    stop("'object' must be a fitted gamlss model.")
  }


  # ==========================================================
  # 2. Parametric design matrix
  # ==========================================================

  X_mu_full <- object$mu.x

  beta_mu_all <- coef(object, "mu")

  if (ncol(X_mu_full) != length(beta_mu_all)) {
    stop("Incompatible dimensions between mu.x and mu coefficients.")
  }

  ind_param_mu <- !is.na(beta_mu_all)

  X_mu_param <- X_mu_full[
    ,
    ind_param_mu,
    drop = FALSE
  ]


  # ==========================================================
  # 3. Smooth components
  # ==========================================================

  terSmo <- gamlss::getSmo(object)

  r <- length(terSmo$smooth)

  if (r == 0) {
    stop("No smooth terms found in the fitted model.")
  }

  if (length(terSmo$sp) != r) {
    stop(
      "The number of smoothing parameters must equal ",
      "the number of smooth terms."
    )
  }


  # ==========================================================
  # 4. Smooth design matrix
  # ==========================================================

  X_smooth_full <- model.matrix(terSmo)

  if (ncol(X_smooth_full) < 2) {
    stop("The smooth design matrix has no smooth columns.")
  }

  # Remove the intercept
  X_smooth <- X_smooth_full[
    ,
    -1,
    drop = FALSE
  ]


  # ==========================================================
  # 5. Complete design matrix
  # ==========================================================

  X <- cbind(
    X_mu_param,
    X_smooth
  )


  # ==========================================================
  # 6. Penalty matrices
  # ==========================================================

  p <- ncol(X_mu_param)

  O <- matrix(0, p, p)

  K_smooth <- vector("list", r)

  S_list <- vector("list", r)

  for (j in seq_len(r)) {

    smooth_j <- terSmo$smooth[[j]]

    # One penalty matrix per smooth term
    if (length(smooth_j$S) != 1) {
      stop(
        "Smooth term ", j,
        " has more than one penalty matrix."
      )
    }

    S_j <- as.matrix(smooth_j$S[[1]])

    # Dimension of the penalty matrix
    if (nrow(S_j) != ncol(S_j)) {
      stop("Penalty matrix must be square.")
    }

    # Apply smoothing parameter
    K_smooth[[j]] <- terSmo$sp[j] * S_j

    S_list[[j]] <- S_j
  }


  # ==========================================================
  # 7. Complete penalty matrix
  # ==========================================================

  K <- as.matrix(
    do.call(
      Matrix::bdiag,
      c(list(O), K_smooth)
    )
  )

  if (ncol(X) != nrow(K)) {
    stop(
      "Incompatible dimensions between the design ",
      "matrix and penalty matrix. Check the smooth ",
      "design matrix and penalty matrices."
    )
  }


  # ==========================================================
  # 8. Estimate sigma
  # ==========================================================

  sigma_hat <- fitted(object, "sigma")[1]


  # ==========================================================
  # 9. Hat matrix
  # ==========================================================

  A <- crossprod(X) + sigma_hat^2 * K

  H <- X %*% solve(A, t(X))


  # ==========================================================
  # 10. Leverage measures
  # ==========================================================

  h <- diag(H)



  # ==========================================================
  # 11. Return results
  # ==========================================================

  return(
    list(
      H = H,
      leverage = h
    )
  )
}
