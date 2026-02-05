#' Gaussian Process–Controlled B-Spline Surface
#'
#' An emulator based on the Gaussian Process–controlled B-spline surface model
#' of Li et al. (2026). The response surface is represented using a
#' tensor-product B-spline basis, with a Gaussian process prior placed on the
#' spline coefficients. Model parameters are estimated via profile likelihood,
#' and the number of spline basis functions per dimension can be selected using
#' a sequential knot number selection (SKNS) procedure based on AIC.
#'
#' @param X A numeric matrix or data frame of predictors, scaled to lie in
#'   \eqn{(0,1)} in each dimension.
#' @param y A numeric response vector of length \eqn{n}.
#' @param m_min Integer or integer vector of length \eqn{p}. Initial (and minimum)
#'   number of B-spline basis functions per dimension. If scalar, the value is
#'   recycled across dimensions.
#' @param m_max Integer or integer vector of length \eqn{p}. Maximum number of
#'   B-spline basis functions per dimension allowed during SKNS.
#' @param degree Integer spline degree used in \code{splines::bs}. Default is 3
#'   (cubic splines).
#' @param mean_fn Mean function specification. Either \code{"linear"} (uses
#'   \code{Fmat = X} with no intercept) or \code{"none"} (no explicit mean function).
#'   Intercepts are handled implicitly through the spline basis.
#' @param tau Nonnegative scalar specifying the noise-to-signal ratio
#'   \eqn{\tau = \delta} in the model. The observation noise variance is
#'   \eqn{\sigma^2 \tau}. Currently treated as fixed.
#' @param psi_init Optional numeric vector of length \eqn{p} giving initial values
#'   for the log lengthscale parameters of the Gaussian process prior on spline
#'   coefficients. If \code{NULL}, a weak default is used.
#' @param warn_complexity Logical; if \code{TRUE}, prints a warning comparing the
#'   rough cubic-time cost proxy \eqn{n + M^3} (where \eqn{M = \prod_j m_j}) to the
#'   \eqn{n^3} cost of a standard Gaussian process.
#' @param method Optimization method passed to \code{\link[stats]{optim}} for
#'   profile likelihood maximization. Default is \code{"L-BFGS-B"}.
#' @param lower,upper Optional numeric vectors specifying lower and upper bounds
#'   for the optimization parameters \code{psi}. If \code{NULL}, reasonable
#'   defaults are used.
#' @param verbose Logical; should progress and diagnostic information be printed?
#' @param skns Logical; if \code{TRUE}, performs sequential knot number selection
#'   (SKNS) over dimensions to choose the number of spline basis functions.
#' @param m_grid Optional integer vector specifying candidate knot counts to try
#'   during SKNS. If provided, the same grid is used for all dimensions (subject
#'   to bounds \code{m_min} and \code{m_max}).
#' @param m_by Integer step size used to construct candidate knot grids when
#'   \code{m_grid} is \code{NULL}. Larger values result in coarser (faster)
#'   searches. Default is 3.
#' @param ... Additional arguments passed to \code{\link[stats]{optim}}.
#'
#' @details
#' The model takes the form
#' \deqn{y = F\beta + U\gamma + \varepsilon,}
#' where \eqn{U} is a tensor-product B-spline basis matrix and the spline
#' coefficients satisfy \eqn{\gamma \sim \mathcal{N}(0, \sigma^2 R(\psi))}.
#' The covariance \eqn{R(\psi)} is constructed as a separable Gaussian process
#' over the spline coefficient lattice, with one lengthscale parameter per input
#' dimension. The regression coefficients \eqn{\beta} and variance \eqn{\sigma^2}
#' are estimated in closed form, while the covariance hyperparameters \eqn{\psi}
#' are estimated by maximizing the profile likelihood.
#'
#' Sequential knot number selection (SKNS) proceeds dimension by dimension,
#' holding all other knot counts fixed while evaluating candidate values for a
#' single dimension, and selecting the value that maximizes the profile
#' likelihood. This greedy procedure provides a data-driven way to control
#' spline resolution while mitigating the curse of dimensionality.
#'
#' This implementation currently uses dense tensor-product basis matrices and
#' dense linear algebra, and is intended for low- to moderate-dimensional
#' problems. Future versions may exploit sparsity and Kronecker structure to
#' improve scalability.
#'
#' @references
#' Li, Y., Tian, Y., Mo, H., & Du, S. (2026).
#' Gaussian Process Controlled B-Spline Surface.
#' \emph{INFORMS Journal on Data Science}.
#'
#' @examples
#' \dontrun{
#' X <- lhs::maximinLHS(200, 2)
#' f <- function(x) x[1]^2 + x[1]*x[2]
#' y <- apply(X, 1, f) + rnorm(200, 0, 0.05)
#'
#' fit <- gpbss(X, y)
#' preds <- predict(fit, X)
#' }
#'
#' @export
gpbss <- function(X, y,
                      m_min = 4L,
                      m_max = 20L,
                      degree = 3L,
                      mean_fn = c("none", "linear"),
                      tau = 1e-6,
                      psi_init = NULL,
                      warn_complexity = TRUE,
                      method = "L-BFGS-B",
                      lower = NULL,
                      upper = NULL,
                      verbose = FALSE,
                      skns = TRUE,
                      m_grid = NULL,
                      m_by = 3,
                      ...) {

  mean_fn <- match.arg(mean_fn)

  ## ---- coerce / validate ----
  if (is.null(X)) stop("X must be provided.", call. = FALSE)
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  if (!is.matrix(X)) stop("X must be a matrix (or a vector coerced to matrix).", call. = FALSE)
  storage.mode(X) <- "double"

  n <- nrow(X)
  p <- ncol(X)

  if (is.null(y) || length(y) != n) stop("y must be a numeric vector of length nrow(X).", call. = FALSE)
  y <- as.numeric(y)

  m_min <- as.integer(m_min)
  if (length(m_min) == 1L) m_min <- rep(m_min, p)
  if (length(m_min) != p) stop("m_min must be length 1 or length ncol(X).", call. = FALSE)

  m_max <- as.integer(m_max)
  if (length(m_max) == 1L) m_max <- rep(m_max, p)
  if (length(m_max) != p) stop("m_max must be length 1 or length ncol(X).", call. = FALSE)
  if (any(m_max < m_min)) stop("All m_max must be >= m_min.", call. = FALSE)

  degree <- as.integer(degree)
  if (degree < 0L) stop("degree must be nonnegative.", call. = FALSE)

  if (!is.numeric(tau) || length(tau) != 1L || tau < 0) stop("tau must be a nonnegative scalar.", call. = FALSE)

  ## ---- mean/trend design Fmat (NO intercept) ----
  if (mean_fn == "linear") {
    Fmat <- as.matrix(X)
    colnames(Fmat) <- paste0("x", seq_len(p))
  } else {
    Fmat <- matrix(numeric(0), nrow = n, ncol = 0)
  }

  ## ---- default psi_init and bounds ----
  if (is.null(psi_init)) psi_init <- rep(log(0.5), p)  # log lengthscales
  psi_init <- as.numeric(psi_init)
  if (length(psi_init) != p) stop("psi_init must have length p.", call. = FALSE)

  if (is.null(lower)) lower <- rep(log(1e-3), p)
  if (is.null(upper)) upper <- rep(log(1e3),  p)

  ## ---- helpers: build 1D bases, tensor U, and R_fun ----
  build_B_list <- function(m_vec) {
    B_list <- vector("list", p)
    for (j in seq_len(p)) {
      Bj <- splines::bs(
        x = X[, j],
        df = as.integer(m_vec[j]),
        degree = degree,
        intercept = TRUE
      )
      Bj <- as.matrix(Bj)
      B_list[[j]] <- Bj
    }
    B_list
  }

  build_U <- function(B_list) {
    m_eff <- vapply(B_list, ncol, integer(1))
    M <- prod(m_eff)

    U <- matrix(0, nrow = n, ncol = M)
    for (i in seq_len(n)) {
      ui <- B_list[[1]][i, ]
      if (p >= 2) {
        for (j in 2:p) ui <- kronecker(ui, B_list[[j]][i, ])
      }
      U[i, ] <- ui
    }
    list(U = U, m_eff = m_eff, M = M)
  }

  make_R_fun <- function(m_eff) {
    m_eff_local <- m_eff
    p_local <- p
    function(psi) {
      psi <- as.numeric(psi)
      if (length(psi) != p_local) stop("psi length mismatch inside R_fun.", call. = FALSE)
      ell <- exp(psi)

      R <- NULL
      for (j in seq_len(p_local)) {
        idx <- seq_len(m_eff_local[j])
        D2 <- (outer(idx, idx, "-"))^2
        Rj <- exp(-0.5 * D2 / (ell[j]^2))
        R <- if (is.null(R)) Rj else kronecker(R, Rj)
      }
      R
    }
  }

  ## ---- SKNS candidate grid ----
  if (is.null(m_grid)) {
    # default m_by if user didn't provide one
    if (is.null(m_by)) m_by <- 2L

    m_by <- as.integer(m_by)
    if (length(m_by) == 1L) m_by <- rep(m_by, p)
    if (length(m_by) != p) stop("m_by must be length 1 or length p.", call. = FALSE)

    m_grid_list <- vector("list", p)
    for (j in seq_len(p)) {
      if (m_by[j] <= 0L) stop("m_by must be positive.", call. = FALSE)
      grid <- seq(m_min[j], m_max[j], by = m_by[j])
      if (tail(grid, 1L) != m_max[j]) grid <- c(grid, m_max[j])
      m_grid_list[[j]] <- unique(as.integer(grid))
    }

  } else {
    # If user supplies a single vector, reuse it for every dim but clamp to [m_min, m_max]
    if (is.numeric(m_grid) && is.vector(m_grid)) {
      m_grid_list <- lapply(seq_len(p), function(j) {
        g <- as.integer(m_grid)
        g <- g[g >= m_min[j] & g <= m_max[j]]
        if (length(g) == 0) g <- m_min[j]
        unique(g)
      })
    } else {
      stop("m_grid must be NULL or an integer vector.", call. = FALSE)
    }
  }

  ## ---- SKNS loop ----
  m_sel <- m_min
  skns_trace <- list()
  psi_warm <- psi_init

  if (isTRUE(skns) && p >= 1) {
    for (j in seq_len(p)) {

      cand <- m_grid_list[[j]]
      best_ll <- -Inf
      best_mj <- m_sel[j]
      best_fit <- NULL

      if (verbose) cat("SKNS dim", j, "candidates:", paste(cand, collapse = ", "), "\n")

      for (mj in cand) {
        m_try <- m_sel
        m_try[j] <- mj

        B_try <- build_B_list(m_try)
        tmpU <- build_U(B_try)
        U_try <- tmpU$U
        m_eff_try <- tmpU$m_eff
        M_try <- tmpU$M

        if (warn_complexity) {
          cost_this <- n + M_try^3
          cost_gp   <- n^3
          ratio <- cost_this / cost_gp
          if (verbose) {
            cat(sprintf("  mj=%d -> M=%d (%s) ; (n+M^3)/n^3 = %.3g\n",
                        mj, M_try, paste(m_eff_try, collapse="x"), ratio))
          }
        }

        R_fun_try <- make_R_fun(m_eff_try)

        fit_try <- fit_gpbss_profile(
          y = y, U = U_try, Fmat = Fmat, tau = tau, R_fun = R_fun_try,
          psi_init = psi_warm,
          method = method,
          lower = lower, upper = upper,
          jitter = 1e-8,
          sigma2_df = "n",
          verbose = FALSE,
          ...
        )

        ll_try <- fit_try$logLik
        if (is.finite(ll_try) && ll_try > best_ll) {
          best_ll <- ll_try
          best_mj <- mj
          best_fit <- fit_try
        }
      }

      # Commit best mj for this dimension
      m_sel[j] <- best_mj
      if (!is.null(best_fit)) psi_warm <- best_fit$psi_hat

      skns_trace[[j]] <- list(
        dim = j,
        chosen_m = best_mj,
        logLik = best_ll,
        psi_hat = if (!is.null(best_fit)) best_fit$psi_hat else NA
      )

      if (verbose) cat("  -> chose mj =", best_mj, " (logLik =", best_ll, ")\n")
    }
  }

  ## ---- final fit using selected m ----
  B_list <- build_B_list(m_sel)
  tmpU <- build_U(B_list)
  U <- tmpU$U
  m_eff <- tmpU$m_eff
  M <- tmpU$M

  cost_this <- n + M^3
  cost_gp   <- n^3
  ratio <- cost_this / cost_gp

  msg <- sprintf(
    paste0(
      "Dense tensor basis: M = prod(m_j) = %d (m_eff = %s). ",
      "Rough cubic proxies: n + M^3 = %.3e vs n^3 = %.3e (ratio = %.3g).\n"
    ),
    M, paste(m_eff, collapse = "x"), cost_this, cost_gp, ratio
  )

  if (verbose) {
    cat(msg)
  }

  if (ratio > 1 & warn_complexity) {
    warning(
      paste0(
        "Dense spline basis may be computationally inefficient relative to a full GP.\n",
        msg
      ),
      call. = FALSE
    )
  }

  R_fun <- make_R_fun(m_eff)

  fit <- fit_gpbss_profile(
    y = y, U = U, Fmat = Fmat, tau = tau, R_fun = R_fun,
    psi_init = psi_warm,
    method = method,
    lower = lower, upper = upper,
    jitter = 1e-8,
    sigma2_df = "n",
    verbose = verbose,
    ...
  )

  out <- list(
    X = X,
    y = y,
    n = n,
    p = p,
    degree = degree,
    mean_fn = mean_fn,
    tau = tau,
    m_min = m_min,
    m_max = m_max,
    m_sel = m_sel,
    m_eff = m_eff,
    M = M,
    B_list = B_list,
    Fmat = Fmat,
    U = U,
    R_fun = R_fun,
    fit = fit,
    skns_trace = skns_trace
  )
  class(out) <- "gpbspline"
  out
}



#' Profile log-likelihood for GP-on-spline-coefficients model
#'
#' Computes the Gaussian profile log-likelihood for the model
#' \deqn{y = F\beta + U\gamma + \varepsilon, \quad
#' \gamma \sim N(0, \sigma^2 R(\psi)), \quad \varepsilon \sim N(0, \sigma^2 \tau I).}
#'
#' For fixed hyperparameters \eqn{\psi} and \eqn{\tau}, the MLEs \eqn{\hat\beta} and
#' \eqn{\hat\sigma^2} are available in closed form via generalized least squares (GLS).
#' This function returns the *profile* log-likelihood after plugging these in.
#'
#' NOTE: This implementation forms \eqn{\Sigma_0 = U R(\psi) U^\top + \tau I} explicitly
#' and uses a Cholesky factorization. That is correct but not scalable for large n/M.
#'
#' @param psi Numeric vector of hyperparameters (passed to \code{R_fun}).
#' @param y Numeric response vector (length n).
#' @param U Numeric matrix (n x M) tensor-product spline basis (or any basis).
#' @param Fmat Numeric matrix (n x q) mean/trend design (e.g., intercept-only or linear).
#' @param tau Nonnegative scalar noise-to-signal ratio. Noise variance is \eqn{\sigma^2\tau}.
#' @param R_fun Function taking \code{psi} and returning either:
#'   (i) an MxM matrix \code{R}, or
#'   (ii) a list with element \code{R} containing the MxM matrix.
#' @param jitter Nonnegative scalar jitter added to the diagonal of \eqn{\Sigma_0}
#'   if needed for numerical stability.
#' @param sigma2_df Degrees-of-freedom divisor for \eqn{\hat\sigma^2}. Choose \code{"n"}
#'   (default) to divide by n (MLE) or \code{"n-q"} to divide by n-q.
#' @param return_details Logical; if TRUE, returns \eqn{\hat\beta}, \eqn{\hat\sigma^2},
#'   \eqn{\Sigma_0^{-1}} solves, etc.
#'
#' @return If \code{return_details=FALSE}: a single numeric value (profile logLik).
#'   If TRUE: a list with elements \code{logLik}, \code{beta_hat}, \code{sigma2_hat},
#'   \code{tau}, \code{psi}, and additional intermediates.
#'
#' @examples
#' # Minimal R_fun: independent coefficients (R = I)
#' R_fun <- function(psi) diag(1, 10)
#' n <- 50; M <- 10
#' U <- matrix(rnorm(n*M), n, M)
#' Fmat <- cbind(1)
#' y <- rnorm(n)
#' profile_loglik_gpbss(psi = numeric(0), y, U, Fmat, tau = 0.1, R_fun = R_fun)
#'
#' @export
profile_loglik_gpbss <- function(psi, y, U, Fmat, tau, R_fun,
                                 jitter = 1e-8,
                                 sigma2_df = c("n", "n-q"),
                                 return_details = FALSE) {

  sigma2_df <- match.arg(sigma2_df)

  ## ---- checks ----
  if (is.null(y) || is.null(U) || is.null(Fmat)) stop("y, U, and Fmat must be provided.", call. = FALSE)
  y <- as.numeric(y)
  if (!is.matrix(U)) U <- as.matrix(U)
  if (!is.matrix(Fmat)) Fmat <- as.matrix(Fmat)

  n <- length(y)
  if (nrow(U) != n) stop("nrow(U) must equal length(y).", call. = FALSE)
  if (nrow(Fmat) != n) stop("nrow(Fmat) must equal length(y).", call. = FALSE)
  if (!is.numeric(tau) || length(tau) != 1L || tau < 0) stop("tau must be a nonnegative scalar.", call. = FALSE)

  ## ---- build R(psi) ----
  Robj <- R_fun(psi)
  R <- if (is.list(Robj)) Robj$R else Robj
  if (!is.matrix(R)) R <- as.matrix(R)
  M <- ncol(U)
  if (!all(dim(R) == c(M, M))) stop("R(psi) must be an MxM matrix with M = ncol(U).", call. = FALSE)

  ## ---- Sigma0 = U R U' + tau I ----
  # Work with Sigma0 (no sigma^2 factor). Profile over sigma^2 afterwards.
  Sigma0 <- U %*% R %*% t(U)
  diag(Sigma0) <- diag(Sigma0) + tau

  ## ---- Cholesky with adaptive jitter ----
  # chol() requires symmetric PD; numerical issues happen.
  chol_try <- function(jit) {
    S <- Sigma0
    if (jit > 0) diag(S) <- diag(S) + jit
    chol(S)
  }

  L <- try(chol_try(0), silent = TRUE)
  jit_used <- 0
  if (inherits(L, "try-error")) {
    jit_used <- jitter
    ok <- FALSE
    for (k in 0:8) {
      L <- try(chol_try(jit_used), silent = TRUE)
      if (!inherits(L, "try-error")) { ok <- TRUE; break }
      jit_used <- jit_used * 10
    }
    if (!ok) stop("Cholesky failed for Sigma0 even after adding jitter.", call. = FALSE)
  }

  ## ---- helper: solve Sigma0^{-1} v via Cholesky ----
  solve_Sigma0 <- function(v) {
    # Solve Sigma0 x = v
    backsolve(L, forwardsolve(t(L), v))
  }

  ## ---- GLS beta_hat ----
  # beta_hat = (F' S^{-1} F)^{-1} F' S^{-1} y
  q <- ncol(Fmat)

  Sinv_y <- solve_Sigma0(y)

  if (q == 0L) {
    beta_hat <- numeric(0)
    r <- y
  } else {
    Sinv_F <- solve_Sigma0(Fmat)  # n x q

    Ft_Sinv_F <- crossprod(Fmat, Sinv_F)  # q x q
    Ft_Sinv_y <- crossprod(Fmat, Sinv_y)  # q x 1

    # Solve (Ft_Sinv_F) beta = Ft_Sinv_y
    ridge <- 1e-10
    beta_hat <- try(
      solve(Ft_Sinv_F + diag(ridge, nrow(Ft_Sinv_F)), Ft_Sinv_y),
      silent = TRUE
    )
    if (inherits(beta_hat, "try-error")) {
      stop("Failed to solve for beta_hat; Ft_Sinv_F may be singular.", call. = FALSE)
    }

    r <- y - drop(Fmat %*% beta_hat)
  }

  Sinv_r <- solve_Sigma0(r)
  quad <- drop(crossprod(r, Sinv_r)) # r' Sigma0^{-1} r

  denom <- if (sigma2_df == "n") n else max(1L, n - q)
  sigma2_hat <- quad / denom
  if (!is.finite(sigma2_hat) || sigma2_hat <= 0) {
    stop("Nonpositive or non-finite sigma2_hat encountered.", call. = FALSE)
  }

  ## ---- log|Sigma0| from Cholesky ----
  # Sigma0 = t(L) %*% L with chol() in R returns upper-tri by default
  # logdet = 2 * sum(log(diag(L)))
  logdet_Sigma0 <- 2 * sum(log(diag(L)))

  ## ---- profile log-likelihood ----
  # Up to additive constant:  -0.5 [ n log(sigma2_hat) + log|Sigma0| + denom ]
  # If using MLE denom=n, the constant term is n; for denom=n-q it becomes (n-q),
  # but since we're mostly doing optimization in psi, constant doesn't matter.
  ll <- -0.5 * (denom * log(sigma2_hat) + logdet_Sigma0 + denom)

  if (!return_details) return(ll)

  list(
    logLik = ll,
    beta_hat = drop(beta_hat),
    sigma2_hat = sigma2_hat,
    tau = tau,
    psi = psi,
    jitter_used = jit_used,
    logdet_Sigma0 = logdet_Sigma0,
    quad = quad,
    denom = denom
  )
}


#' Fit GP-on-spline-coefficients model by profile likelihood
#'
#' Optimizes hyperparameters \eqn{\psi} by maximizing the profile log-likelihood
#' returned by \code{\link{profile_loglik_gpbss}}, and returns the resulting
#' \eqn{\hat\psi}, \eqn{\hat\beta}, and \eqn{\hat\sigma^2}.
#'
#' @param y Numeric response vector (length n).
#' @param U Numeric matrix (n x M) tensor-product spline basis.
#' @param Fmat Numeric matrix (n x q) mean/trend design.
#' @param tau Nonnegative scalar noise-to-signal ratio. (For now treated as fixed.)
#' @param R_fun Function mapping \code{psi} to an MxM prior covariance \code{R(psi)}.
#' @param psi_init Numeric vector initial guess for \code{psi}.
#' @param method Optimization method passed to \code{optim}. Default \code{"L-BFGS-B"}.
#' @param lower,upper Optional bounds for \code{psi} (used if method supports bounds).
#' @param jitter Jitter passed to \code{profile_loglik_gpbss}.
#' @param sigma2_df Degrees-of-freedom divisor for \eqn{\hat\sigma^2}, passed through.
#' @param verbose Logical; print optimizer progress and final value.
#' @param ... Additional arguments passed to \code{optim}.
#'
#' @return A list with elements \code{psi_hat}, \code{beta_hat}, \code{sigma2_hat},
#'   \code{logLik}, \code{optim}, and the fixed \code{tau}.
#'
#' @export
fit_gpbss_profile <- function(y, U, Fmat, tau, R_fun,
                              psi_init,
                              method = "L-BFGS-B",
                              lower = NULL, upper = NULL,
                              jitter = 1e-8,
                              sigma2_df = c("n", "n-q"),
                              verbose = FALSE,
                              ...) {

  #sigma2_df <- match.arg(sigma2_df)

  if (is.null(psi_init)) stop("psi_init must be provided.", call. = FALSE)
  psi_init <- as.numeric(psi_init)

  obj_fn <- function(psi) {
    # optim minimizes; we maximize logLik
    ll <- profile_loglik_gpbss(
      psi = psi, y = y, U = U, Fmat = Fmat, tau = tau, R_fun = R_fun,
      jitter = jitter, sigma2_df = sigma2_df, return_details = FALSE
    )
    if (!is.finite(ll)) return(1e30)
    -ll
  }

  if (is.null(lower)) lower <- rep(-Inf, length(psi_init))
  if (is.null(upper)) upper <- rep( Inf, length(psi_init))

  opt <- stats::optim(
    par = psi_init,
    fn = obj_fn,
    method = method,
    lower = lower,
    upper = upper,
    ...
  )

  if (verbose) {
    cat("optim convergence:", opt$convergence, "\n")
    cat("optim value (neg logLik):", opt$value, "\n")
  }

  # Recompute details at optimum
  dets <- profile_loglik_gpbss(
    psi = opt$par, y = y, U = U, Fmat = Fmat, tau = tau, R_fun = R_fun,
    jitter = jitter, sigma2_df = sigma2_df, return_details = TRUE
  )

  list(
    psi_hat = opt$par,
    beta_hat = dets$beta_hat,
    sigma2_hat = dets$sigma2_hat,
    logLik = dets$logLik,
    tau = tau,
    optim = opt
  )
}


#' Posterior Predictive Samples for GP-B-spline Emulator
#'
#' Generates posterior predictive samples for a fitted \code{gpbspline} model.
#' Returns a \code{samples x n_test} matrix, consistent with other \pkg{duqling}
#' emulator \code{predict} methods.
#'
#' @param object An object of class \code{gpbspline} returned by \code{gpbspline()}.
#' @param newdata New input matrix \code{(n_test x p)}. If \code{NULL}, uses training inputs.
#' @param samples Integer number of posterior predictive samples to draw. Default is 1000.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric matrix of dimension \code{samples x n_test}, where each row is
#'   one posterior predictive draw.
#'
#' @export
predict.gpbspline <- function(object, newdata = NULL, samples = 1000, ...) {
  if (!inherits(object, "gpbspline")) stop("object must be of class 'gpbspline'.", call. = FALSE)

  if (is.null(newdata)) newdata <- object$X
  if (is.vector(newdata)) newdata <- matrix(newdata, ncol = object$p)
  if (!is.matrix(newdata)) stop("newdata must be a matrix (or vector coerced to matrix).", call. = FALSE)
  if (ncol(newdata) != object$p) stop("newdata must have ncol equal to training X.", call. = FALSE)

  samples <- as.integer(samples)
  if (length(samples) != 1L || samples < 1L) stop("samples must be a positive integer.", call. = FALSE)

  # --- pull fitted pieces ---
  fit <- object$fit
  if (is.null(fit$beta_hat) || is.null(fit$sigma2_hat) || is.null(fit$psi_hat)) {
    stop("object$fit must contain beta_hat, sigma2_hat, and psi_hat.", call. = FALSE)
  }

  y <- object$y
  X <- object$X
  U_train <- object$U
  F_train <- object$Fmat
  tau <- object$tau
  sigma2_hat <- as.numeric(fit$sigma2_hat)
  beta_hat <- as.numeric(fit$beta_hat)
  psi_hat <- as.numeric(fit$psi_hat)
  R_fun <- object$R_fun
  B_list_train <- object$B_list

  n <- nrow(X)
  n_test <- nrow(newdata)
  p <- object$p

  # --- build F* consistent with training (your convention: no intercept) ---
  if (ncol(F_train) == 0L) {
    F_star <- matrix(numeric(0), nrow = n_test, ncol = 0)
  } else {
    F_star <- as.matrix(newdata)
    if (ncol(F_star) != ncol(F_train)) stop("Fmat dimension mismatch.", call. = FALSE)
  }

  # --- evaluate 1D spline bases at newdata using training bs() attributes ---
  eval_bs_from_attr <- function(B_train, x_new) {
    knots <- attr(B_train, "knots")
    bdry  <- attr(B_train, "Boundary.knots")
    deg   <- attr(B_train, "degree")
    if (is.null(bdry) || length(bdry) != 2L) {
      stop("Stored spline basis is missing Boundary.knots.", call. = FALSE)
    }

    tvec <- c(rep(bdry[1], deg + 1L),
              if (!is.null(knots) && length(knots) > 0) knots else numeric(0),
              rep(bdry[2], deg + 1L))

    Bnew <- splines::splineDesign(knots = tvec, x = x_new, ord = deg + 1L, outer.ok = TRUE)
    Bnew <- as.matrix(Bnew)

    # Safety: splineDesign should match training ncol
    if (ncol(Bnew) != ncol(B_train)) {
      stop("Basis evaluation mismatch: ncol(Bnew) != ncol(B_train).", call. = FALSE)
    }
    Bnew
  }

  B_list_star <- vector("list", p)
  for (j in seq_len(p)) {
    B_list_star[[j]] <- eval_bs_from_attr(B_list_train[[j]], newdata[, j])
  }

  # --- build dense tensor basis U* (same column ordering as training) ---
  m_eff <- vapply(B_list_star, ncol, integer(1))
  M <- prod(m_eff)
  U_star <- matrix(0, nrow = n_test, ncol = M)
  for (i in seq_len(n_test)) {
    ui <- B_list_star[[1]][i, ]
    if (p >= 2) {
      for (j in 2:p) ui <- kronecker(ui, B_list_star[[j]][i, ])
    }
    U_star[i, ] <- ui
  }

  # --- compute posterior for gamma ---
  # A = (R^{-1} + (1/tau) U'U)^{-1}
  R <- R_fun(psi_hat)
  if (!is.matrix(R) || nrow(R) != M || ncol(R) != M) {
    stop("R_fun(psi_hat) must return an M x M matrix matching U*.", call. = FALSE)
  }

  # R^{-1} (jitter if needed)
  Rinv <- try(solve(R), silent = TRUE)
  if (inherits(Rinv, "try-error")) {
    jit <- max(1e-10, 1e-10 * mean(diag(R)))
    Rinv <- solve(R + diag(jit, M))
  }

  UU <- crossprod(U_train)  # M x M
  Kmat <- Rinv + (1 / tau) * UU

  A <- try(solve(Kmat), silent = TRUE)
  if (inherits(A, "try-error")) {
    jit <- max(1e-10, 1e-10 * mean(diag(Kmat)))
    A <- solve(Kmat + diag(jit, M))
  }

  # r = y - F beta
  if (ncol(F_train) == 0L) {
    r <- y
  } else {
    r <- as.numeric(y - drop(F_train %*% beta_hat))
  }

  Ur <- crossprod(U_train, r)             # M x 1
  gamma_hat <- (1 / tau) * (A %*% Ur)     # M x 1

  # predictive mean
  mu <- numeric(n_test)
  if (ncol(F_star) > 0L) mu <- mu + as.numeric(F_star %*% beta_hat)
  mu <- mu + as.numeric(U_star %*% gamma_hat)

  # predictive covariance: sigma2 * (U* A U*' + tau I)
  # (this is the expensive part; OK for now)
  Cov <- sigma2_hat * (U_star %*% A %*% t(U_star) + diag(tau, n_test))

  # draw samples jointly (keeps cross-point dependence)
  preds <- matrix(NA_real_, nrow = samples, ncol = n_test)
  jitC <- max(1e-10, 1e-10 * mean(diag(Cov)))
  L <- try(chol(Cov + diag(jitC, n_test)), silent = TRUE)

  if (!inherits(L, "try-error")) {
    Z <- matrix(rnorm(samples * n_test), nrow = n_test)
    # chol is upper-tri U so U' %*% Z gives covariance
    preds <- t(mu + t(L) %*% Z)
  } else {
    # fallback: eigen (slower but robust)
    ev <- eigen(Cov, symmetric = TRUE)
    vals <- pmax(ev$values, 0)
    sqrtC <- ev$vectors %*% diag(sqrt(vals), n_test) %*% t(ev$vectors)
    Z <- matrix(rnorm(samples * n_test), nrow = n_test)
    preds <- t(mu + sqrtC %*% Z)
  }

  preds
}

#' Plot Diagnostics for GP-B-spline Objects
#'
#' Produces simple diagnostic plots for a fitted \code{gpbspline} model:
#' (1) posterior predictive mean vs observed \code{y} with ~95\% intervals;
#' (2) histogram of residuals with a normal curve overlay.
#'
#' @param x An object of class \code{gpbspline} as returned by \code{gpbspline()}.
#' @param ... Additional plotting arguments (currently ignored).
#'
#' @return Called for its side effect (plots). Invisibly returns \code{NULL}.
#'
#' @seealso \code{\link{predict.gpbspline}}
#'
#' @export
plot.gpbspline <- function(x, ...) {
  if (!inherits(x, "gpbspline")) stop("x must be of class 'gpbspline'.", call. = FALSE)

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))

  preds <- stats::predict(x)  # samples x n
  yhat <- apply(preds, 2, mean)
  yy <- x$y

  ci <- 2 * apply(preds, 2, stats::sd)

  graphics::plot(yy, yhat, pch = 16, xlab = "y", ylab = "E[y | data]")
  graphics::segments(x0 = yy, y0 = yhat - ci, y1 = yhat + ci, col = "orange")
  graphics::points(yy, yhat, pch = 16)
  graphics::abline(0, 1, col = "dodgerblue")

  rr <- yy - yhat
  graphics::hist(rr,
                 breaks = ceiling(length(rr)^0.33 * diff(range(rr)) / (3.5 * stats::sd(rr))),
                 freq = FALSE, main = "", xlab = "residuals")
  graphics::curve(stats::dnorm(x, mean(rr), stats::sd(rr)),
                  add = TRUE, col = "orange", lwd = 2)

  invisible(NULL)
}
