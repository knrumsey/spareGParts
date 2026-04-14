#' Relevant Vector Machine
#'
#' Bayesian relevance vector machine regression with Gaussian radial basis
#' functions and optional Bayesian model averaging over a discrete set of
#' kernel lengthscales.
#'
#' This implementation is included because the faster \pkg{kernlab}
#' implementation does not provide posterior predictive samples in the form
#' needed for \pkg{duqling}.
#'
#' @param X A data frame or matrix of predictors scaled to lie between 0 and 1.
#' @param y A response vector of length \code{nrow(X)}.
#' @param max_basis Maximum number of basis functions retained after optional
#'   LASSO screening when the full basis is too wide.
#' @param qlscale Discrete lengthscale set on the quantile scale, formed from
#'   quantiles of pairwise distances among sampled rows of \code{X}.
#' @param lscale Optional discrete set of kernel lengthscales. If supplied,
#'   this overrides \code{qlscale}.
#' @param lscale_probs Prior probabilities corresponding to the candidate
#'   lengthscales.
#' @param prune_thresh Basis functions with \code{alpha_i > prune_thresh} are
#'   effectively removed from the model to speed computation.
#' @param drop_models Should candidate lengthscale models with negligible
#'   posterior probability be dropped?
#' @param tol Relative tolerance for early stopping of hyperparameter
#'   optimization.
#' @param maxiter Maximum number of EM iterations.
#' @param mc_cores Number of cores used to parallelize across candidate
#'   lengthscales.
#' @param verbose Logical; print progress?
#'
#' @details
#' For each candidate kernel lengthscale, a radial-basis-function design matrix
#' is built using training inputs as candidate centers. If the basis is too
#' large, an initial LASSO screen is used to reduce its width. Conditional on
#' that basis, the RVM hyperparameters are estimated by iterative type-II
#' maximum likelihood updates for the coefficient precisions and observation
#' noise variance. When multiple candidate lengthscales are supplied, posterior
#' model probabilities are computed from the marginal likelihood values and used
#' for Bayesian model averaging in prediction.
#'
#' Predictive draws are for the observed response, so posterior predictive
#' variance includes the fitted residual noise variance.
#'
#' @references
#' Tipping, Michael. "The relevance vector machine." Advances in neural
#' information processing systems 12 (1999).
#'
#' Tipping, Michael E. "Sparse Bayesian learning and the relevance vector
#' machine." Journal of Machine Learning Research 1 (2001): 211-244.
#'
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391 * ((x[1] - 0.4) * (x[2] - 0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- rvm(X, y)
#' @export
rvm <- function(X, y,
                max_basis = 1000,
                qlscale = c(0.2, 0.5),
                lscale = NULL,
                lscale_probs = NULL,
                prune_thresh = 1e6,
                drop_models = TRUE,
                tol = 5e-3,
                maxiter = 2000,
                mc_cores = 1,
                verbose = TRUE) {
  X <- as.matrix(X)

  if (stats::sd(y) == 0) {
    y <- y + stats::rnorm(length(y), 0, 1e-6)
  }

  if (is.null(lscale)) {
    ndist <- min(1000, nrow(X))
    ind_dist <- sample(nrow(X), ndist, replace = FALSE)
    dists <- as.matrix(stats::dist(X[ind_dist, , drop = FALSE]))
    pairwise_dists <- dists[lower.tri(dists)]
    lscale <- as.numeric(stats::quantile(pairwise_dists, probs = qlscale))
  }

  nl <- length(lscale)
  if (is.null(lscale_probs)) {
    lscale_probs <- rep(1 / nl, nl)
  } else {
    lscale_probs <- as.numeric(lscale_probs)
    lscale_probs <- lscale_probs / sum(lscale_probs)
  }

  mc_cores <- min(mc_cores, nl)

  mu_y <- mean(y)
  sigma_y <- stats::sd(y)
  y <- (y - mu_y) / sigma_y

  fit_one_lscale <- function(i) {
    Phi_full <- make_Phi(X, lscale = lscale[i])

    if (ncol(Phi_full) > max_basis + 1) {
      lasso_fit <- glmnet::glmnet(
        Phi_full, y,
        alpha = 1,
        intercept = FALSE,
        lambda.min.ratio = 1e-6,
        nlambda = 100
      )

      coef_path <- as.matrix(stats::coef(lasso_fit))[-1, , drop = FALSE]
      nnz <- apply(coef_path, 2, function(x) sum(x != 0))
      ok <- which(nnz <= max_basis)

      if (length(ok) == 0) {
        coefs <- as.numeric(stats::coef(lasso_fit, s = min(lasso_fit$lambda))[-1])
        selected <- order(abs(coefs), decreasing = TRUE)[seq_len(max_basis)]
      } else {
        lambda_idx <- ok[which.max(nnz[ok])]
        coefs <- as.numeric(stats::coef(lasso_fit, s = lasso_fit$lambda[lambda_idx])[-1])
        selected <- which(coefs != 0)
      }

      selected <- sort(unique(selected))
      Phi_curr <- Phi_full[, selected, drop = FALSE]
    } else {
      Phi_curr <- Phi_full
      selected <- seq_len(ncol(Phi_full))
    }

    fit <- optimize_hyperpars(
      Phi_curr, y,
      prune_thresh = prune_thresh,
      tol = tol,
      maxiter = maxiter,
      verbose = verbose
    )
    fit$selected <- selected
    list(Phi = Phi_curr, Fit = fit)
  }

  if (mc_cores == 1) {
    out_list <- lapply(seq_along(lscale), fit_one_lscale)
  } else {
    out_list <- parallel::mclapply(seq_along(lscale), fit_one_lscale, mc.cores = mc_cores)
  }

  Phi <- lapply(out_list, `[[`, "Phi")
  Fits <- lapply(out_list, `[[`, "Fit")

  log_evidence <- vapply(Fits, function(fit) fit$log_marginal_lik, numeric(1))
  log_prior <- log(lscale_probs)
  max_logev <- max(log_evidence + log_prior)
  post_probs <- exp((log_evidence + log_prior) - max_logev)
  post_probs <- post_probs / sum(post_probs)

  if (drop_models) {
    keep <- which(post_probs > 1e-6)
    Fits <- Fits[keep]
    Phi <- Phi[keep]
    lscale <- lscale[keep]
    lscale_probs <- lscale_probs[keep]
    post_probs <- post_probs[keep]
    post_probs <- post_probs / sum(post_probs)
  }

  out <- list(
    Fits = Fits,
    Phi = Phi,
    lscale = lscale,
    lscale_probs = lscale_probs,
    post_probs = post_probs,
    X_train = X,
    y_train = y,
    y_center = mu_y,
    y_scale = sigma_y,
    call = match.call()
  )
  class(out) <- "rvm"
  out
}


rbf_kernel <- function(x, y, l) {
  exp(-sum((x - y)^2) / (2 * l^2))
}


make_Phi <- function(X, centers = NULL, kernel = rbf_kernel, lscale = 0.1) {
  if (is.null(centers)) {
    centers <- X
  }

  N <- nrow(X)
  M <- nrow(centers)
  Phi <- matrix(NA_real_, nrow = N, ncol = M)

  for (i in seq_len(N)) {
    for (j in seq_len(M)) {
      Phi[i, j] <- kernel(X[i, ], centers[j, ], lscale)
    }
  }

  Phi <- cbind(1, Phi)
  return(Phi)
}


chol_solve_with_jitter <- function(A, b = NULL,
                                   jitter_seq = c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) {
  n <- nrow(A)

  for (jit in jitter_seq) {
    A_jit <- A
    if (jit > 0) {
      diag(A_jit) <- diag(A_jit) + jit
    }

    R <- try(chol(A_jit), silent = TRUE)
    if (!inherits(R, "try-error")) {
      if (is.null(b)) {
        return(list(chol = R, jitter = jit))
      } else {
        x <- backsolve(R, forwardsolve(t(R), b))
        return(list(chol = R, x = x, jitter = jit))
      }
    }
  }

  stop("Cholesky failed even after jitter escalation.")
}


optimize_hyperpars <- function(Phi, y,
                               prune_thresh = 1e6,
                               tol = 1e-4,
                               maxiter = 500,
                               verbose = TRUE) {
  N <- nrow(Phi)
  M <- ncol(Phi)
  pruned <- integer(0)
  alpha <- rep(1, M)
  sigma2 <- 1
  Phi_full <- Phi

  for (iter in seq_len(maxiter)) {
    keep_set <- setdiff(seq_len(M), pruned)
    if (!(1 %in% keep_set)) {
      warning("Intercept was pruned; restoring it.")
      keep_set <- sort(unique(c(1, keep_set)))
    }

    Phi <- Phi_full[, keep_set, drop = FALSE]
    Sigma_inv <- diag(alpha[keep_set], nrow = length(keep_set)) + crossprod(Phi) / sigma2

    chol_out <- chol_solve_with_jitter(Sigma_inv)
    R <- chol_out$chol
    Sigma <- chol2inv(R)
    mu <- Sigma %*% (crossprod(Phi, y) / sigma2)

    gamma <- 1 - alpha[keep_set] * diag(Sigma)

    alpha_new <- rep(1e9, M)
    alpha_new[keep_set] <- as.numeric(gamma / pmax(mu^2, .Machine$double.eps))
    alpha_new[1] <- 1e-9

    denom <- max(N - sum(gamma), 1e-8)
    sigma2_new <- max(1e-7, sum((y - Phi %*% mu)^2) / denom)

    ind_to_prune <- which(alpha_new > prune_thresh)
    pruned <- union(pruned, ind_to_prune)
    pruned <- setdiff(pruned, 1)

    rel_change_alpha <- max(abs(alpha_new - alpha) / (abs(alpha) + .Machine$double.eps))
    rel_change_sigma <- abs(sigma2_new - sigma2) / (abs(sigma2) + .Machine$double.eps)

    if (verbose && (iter %% 100 == 0)) {
      cat(
        "iteration ", iter,
        "\n\tmax relative change = ",
        max(rel_change_alpha, rel_change_sigma),
        "\n",
        sep = ""
      )
    }

    alpha <- alpha_new
    sigma2 <- sigma2_new

    if (max(rel_change_alpha, rel_change_sigma) < tol) {
      break
    }
  }

  if (iter == maxiter) {
    warning("maxiter reached before convergence was obtained.")
  }

  keep_set <- setdiff(seq_len(M), pruned)
  if (!(1 %in% keep_set)) {
    keep_set <- sort(unique(c(1, keep_set)))
  }

  Phi <- Phi_full[, keep_set, drop = FALSE]
  inv_alpha <- 1 / alpha[keep_set]
  C <- sigma2 * diag(N) + tcrossprod(Phi, Phi * rep(inv_alpha, each = nrow(Phi)))

  cholC_out <- chol_solve_with_jitter(C)
  cholC <- cholC_out$chol
  sCy <- backsolve(cholC, forwardsolve(t(cholC), y))

  logdetC <- 2 * sum(log(diag(cholC)))
  log_marginal_lik <- -0.5 * (N * log(2 * pi) + logdetC + crossprod(y, sCy))

  list(
    mu = mu,
    Sigma = Sigma,
    alpha = alpha,
    sigma2 = sigma2,
    keep_set = keep_set,
    log_marginal_lik = as.numeric(log_marginal_lik),
    iter = iter
  )
}


#' Posterior Predictive Sampling for RVM Objects
#'
#' Generates posterior predictive samples for a fitted \code{rvm} object,
#' allowing for Bayesian model averaging over kernel lengthscales.
#'
#' @param object An object of class \code{rvm} as returned by \code{rvm()}.
#' @param newdata A matrix or data frame of new input locations. If \code{NULL},
#'   predictions are generated for the training data.
#' @param samples Number of posterior predictive samples.
#' @param nugget Logical. Should predictive draws include the fitted residual
#'   noise variance? Defaults to \code{TRUE}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A numeric matrix of dimension \code{samples x nrow(newdata)}
#'   containing posterior predictive draws.
#' @export
predict.rvm <- function(object, newdata = NULL, samples = 1000, nugget = TRUE, ...) {
  if (is.null(newdata)) {
    newdata <- object$X_train
  }
  newdata <- as.matrix(newdata)

  n_models <- length(object$Fits)
  n_test <- nrow(newdata)
  preds <- matrix(NA_real_, nrow = samples, ncol = n_test)

  wts <- object$post_probs
  wts <- wts / sum(wts)

  mu_post <- s2_post <- matrix(NA_real_, nrow = n_models, ncol = n_test)

  for (j in seq_len(n_models)) {
    fit <- object$Fits[[j]]
    lscale <- object$lscale[j]

    mu <- fit$mu
    Sigma <- fit$Sigma
    sigma2 <- fit$sigma2
    keep_set <- fit$keep_set

    K <- make_Phi(newdata, object$X_train, lscale = lscale)
    K <- K[, fit$selected, drop = FALSE]
    K <- K[, keep_set, drop = FALSE]

    mu_post[j, ] <- as.vector(K %*% mu)
    s2_post[j, ] <- rowSums((K %*% Sigma) * K)

    if (nugget) {
      s2_post[j, ] <- s2_post[j, ] + sigma2
    }
  }

  if (n_models == 1) {
    for (i in seq_len(samples)) {
      preds[i, ] <- stats::rnorm(
        n_test,
        mean = mu_post[1, ],
        sd = sqrt(pmax(s2_post[1, ], 0))
      )
    }
  } else {
    mix_ids <- sample(seq_len(n_models), samples, replace = TRUE, prob = wts)
    for (i in seq_len(samples)) {
      j <- mix_ids[i]
      preds[i, ] <- stats::rnorm(
        n_test,
        mean = mu_post[j, ],
        sd = sqrt(pmax(s2_post[j, ], 0))
      )
    }
  }

  if (!is.null(object$y_scale) && !is.null(object$y_center)) {
    preds <- preds * object$y_scale + object$y_center
  }

  preds
}


#' Plot Diagnostics for RVM Objects
#'
#' Plots diagnostic summaries for a fitted \code{rvm} model.
#'
#' @param x An object of class \code{rvm} as returned by \code{rvm()}.
#' @param ... Additional plotting arguments.
#' @export
plot.rvm <- function(x, ...) {
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar), add = TRUE)

  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))

  preds <- stats::predict(x)
  yhat <- colMeans(preds)
  yy <- x$y_center + x$y_scale * x$y_train
  ci <- 2 * apply(preds, 2, stats::sd)

  plot(yy, yhat, pch = 16, xlab = "y", ...)
  graphics::segments(x0 = yy, y0 = yhat - ci, y1 = yhat + ci, col = "orange")
  graphics::points(yy, yhat, pch = 16)
  graphics::abline(0, 1, col = "dodgerblue")

  rr <- yy - yhat
  graphics::hist(
    rr,
    breaks = ceiling(length(rr)^0.33 * diff(range(rr)) / (3.5 * stats::sd(rr))),
    freq = FALSE
  )
  graphics::curve(
    stats::dnorm(x, mean(rr), stats::sd(rr)),
    add = TRUE,
    col = "orange",
    lwd = 2
  )
}
