#' Relevant Vector Machine
#'
#' A Bayesian RVM implementation (Tipping 1999). This implementation exists only because the (much more efficient) kernlab implementation doesn't provide probabilistic predictions (or the capability to do so).
#'
#' @param X A dataframe or matrix of predictors scaled to be between 0 and 1
#' @param y a reponse vector of length n
#' @param max_basis Bounds complexity by performing LASSO on initial basis functions when ncol(X) exceeds max_basis.
#' @param qlscale Discrete lengthscale set on quantile scale (quantiles of the distribution of pairwise X distances).
#' @param lscale Discrete lengthscale set. Overrides qlscale when specified.
#' @param lscale_probs Prior probabilities corresponding to the lengthscale set.
#' @param prune_thresh When alpha_i exceed prune_thresh, we effectively set it to infinity speeding up future matrix inverse solves.
#' @param drop_models Should models (corresponding to lscale) be dropped if they have sufficiently small posterior probability?
#' @param tol Tolerance for early stopping of hyperparameter optimization.
#' @param maxiter Number of iterations of EM algorithm before stopping.
#' @param mc_cores How many cores to use (for parallelizing over various lengthscales)
#' @details Algorithm has complexity O(nm^2). Candidate points are greedily selected to maximize a scoring criterion (Eq. 8 in Keerthi & Chu), conditional on kernel parameters. Once a subset is obtained, the GPfit package is used to estimate kernel parameters and the process repeats.
#' @references
#' Tipping, Michael. "The relevance vector machine." Advances in neural information processing systems 12 (1999).
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- rvm(X, y)
#' @export
rvm <- function(X, y, max_basis=1000, qlscale=c(0.2, 0.5), lscale=NULL, lscale_probs=NULL, prune_thresh=1e6, drop_models=TRUE, tol=5e-3, maxiter=2000, mc_cores=1){
  if(sd(y) == 0) y <- y + rnorm(y, 0, 1e-6)

  if(is.null(lscale)){
    ndist <- min(1000, nrow(X))
    ind_dist <- sample(nrow(X), ndist, replace=FALSE)
    dists <- as.matrix(dist(X[ind_dist,,drop=FALSE]))
    pairwise_dists <- dists[lower.tri(dists)]
    lscale <- quantile(pairwise_dists, probs = qlscale)
  }

  nl <- length(lscale)
  if(is.null(lscale_probs)){
    lscale_probs <- rep(1/nl, nl)
  }
  mc_cores <- min(mc_cores, length(lscale))

  # Get initial hyperparameter estimates
  mu_y <- mean(y)
  sigma_y <- sd(y)
  y <- (y-mu_y) / sigma_y
  sigma0 <- 1
  #alpha0 should happen separately for each lengthscale

  selected <- "should dissapear"

  # Send phi matrices to helper function
  Phi <- Fits <- list()
  if(mc_cores == 1){
    # Just loop over lscale's
    for(i in seq_along(lscale)){
      Phi_full <- make_Phi(X, lscale=lscale[i])
      # LASSO screening (optional: only if Phi is wide)
      if(ncol(Phi_full) > max_basis + 1){
        # Fit LASSO path
        lasso_fit <- glmnet(Phi_full, y, alpha=1, intercept=FALSE, lambda.min.ratio = 1e-6, nlambda=100)
        # Find lambda that yields at most max_basis nonzero coefficients
        nnz <- apply(coef(lasso_fit)[-1, , drop=FALSE], 2, function(x) sum(x != 0))
        lambda_idx <- max(which(nnz <= max_basis))
        coefs <- coef(lasso_fit, s=lasso_fit$lambda[lambda_idx])[-1]
        selected <- which(as.numeric(coefs) != 0)
        Phi_curr <- Phi_full[, selected, drop=FALSE]
      } else {
        Phi_curr <- Phi_full
        selected <- seq_len(ncol(Phi_full))
      }
      Phi[[i]] <- Phi_curr

      #Phi[[i]] <- make_Phi(X, lscale=lscale[i])
      # Do lasso here?
      fit <- optimize_hyperpars(Phi[[i]], y, prune_thresh, tol, maxiter)
      fit$selected <- selected
      Fits[[i]] <- fit
    }
  }else{
    # Use mclapply
    out <- mclapply(seq_along(lscale), function(i) {
      Phi_full <- make_Phi(X, lscale = lscale[i])
      if (ncol(Phi_full) > max_basis + 1) {
        # LASSO screening
        lasso_fit <- glmnet(Phi_full, y, alpha = 1, intercept = FALSE,
                            lambda.min.ratio = 1e-6, nlambda = 100)
        nnz <- apply(coef(lasso_fit)[-1, , drop = FALSE], 2, function(x) sum(x != 0))
        ok <- which(nnz <= max_basis)
        if (length(ok) == 0) {
          # Fallback: keep max_basis largest by absolute coefficient at minimal lambda
          coefs <- coef(lasso_fit, s = min(lasso_fit$lambda))[-1]
          selected <- order(abs(coefs), decreasing = TRUE)[seq_len(max_basis)]
        } else {
          lambda_idx <- ok[which.max(nnz[ok])]
          coefs <- coef(lasso_fit, s = lasso_fit$lambda[lambda_idx])[-1]
          selected <- which(as.numeric(coefs) != 0)
        }
        Phi_curr <- Phi_full[, selected, drop = FALSE]
      } else {
        Phi_curr <- Phi_full
        selected <- seq_len(ncol(Phi_full))
      }
      # Fit RVM on the (potentially reduced) Phi
      fit <- optimize_hyperpars(Phi_curr, y, prune_thresh, tol, maxiter)
      fit$selected <- selected
      list(Phi = Phi_curr, Fit = fit)
    }, mc.cores = mc_cores)

    Phi <- lapply(out, `[[`, "Phi")
    Fits <- lapply(out, `[[`, "Fit")
  }

  # Evaluate posterior probabilities for each lscale
  log_evidence <- sapply(Fits, function(fit) fit$log_marginal_lik)
  log_prior <- log(lscale_probs)
  max_logev <- max(log_evidence + log_prior)  # for numerical stability
  post_probs <- exp((log_evidence + log_prior) - max_logev)
  post_probs <- post_probs / sum(post_probs)

  # Drop negligible models
  if(drop_models){
    keep <- which(post_probs > 1e-6)
    Fits <- Fits[keep]
    lscale <- lscale[keep]
    post_probs <- post_probs[keep]
    # Renormalize, just in case
    post_probs <- post_probs / sum(post_probs)
  }

  out <- list(
    Fits = Fits,
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
  return(out)
}


# Default RBF kernel (user can override)
rbf_kernel <- function(x, y, l) exp(-sum((x-y)^2) / (2*l^2))

make_Phi <- function(X, centers = NULL, kernel = rbf_kernel, lscale = 0.1) {
  if(is.null(centers)) centers <- X
  N <- nrow(X)
  M <- nrow(centers)
  Phi <- matrix(NA, nrow = N, ncol = M)
  for (i in 1:N) {
    for (j in 1:M) {
      Phi[i, j] <- kernel(X[i, ], centers[j, ], lscale)
    }
  }
  Phi <- cbind(1, Phi)
}


optimize_hyperpars <- function(Phi, y, prune_thresh=1e6, tol=1e-4, maxiter=500) {
  N <- nrow(Phi)
  M <- ncol(Phi)
  pruned <- NULL
  alpha <- rep(1, M)
  sigma2 <- 1
  Phi_full <- Phi
  # EM loop
  for(iter in 1:maxiter){
    # Compute Sigma and mu (posterior covariance and mean of weights)
    keep_set <- setdiff(1:M, pruned)
    Phi <- Phi_full[, keep_set, drop=FALSE]
    Sigma_inv <- diag(alpha[keep_set]) + crossprod(Phi) / sigma2

    # Compute Sigma with try-catch and increasing jitter if needed
    for (jit in c(0, 1e-8, 1e-6, 1e-4)) {
      Sigma_inv_jit <- Sigma_inv + diag(jit, ncol(Sigma_inv))
      Sigma <- try(solve(Sigma_inv_jit), silent = TRUE)
      if (!inherits(Sigma, "try-error")) {
        break
      }
    }
    mu <- Sigma %*% (t(Phi) %*% y) / sigma2
    # Compute gamma values (sparsity relevance)
    gamma <- 1 - alpha[keep_set] * diag(Sigma)
    # Update alpha and sigma2
    alpha_new <- rep(1e9, M)
    alpha_new[keep_set] <- gamma / (mu^2)
    sigma2_new <- sum((y - Phi %*% mu)^2) / (N - sum(gamma))

    # Prune step
    ind_to_prune <- which(alpha_new > prune_thresh)
    pruned <- c(pruned, ind_to_prune)

    # Convergence check
    rel_change_alpha <- max(abs(alpha_new - alpha) / (abs(alpha) + .Machine$double.eps))
    rel_change_sigma <- abs(sigma2_new - sigma2) / (abs(sigma2) + .Machine$double.eps)
    if(rel_change_alpha < tol && rel_change_sigma < tol) break
    print(max(rel_change_alpha, rel_change_sigma))
    alpha <- alpha_new
    sigma2 <- sigma2_new
  }
  if(iter == maxiter){
    warning("maxiter reached before convergence was obtained.")
  }
  # Compute log marginal likelihood (see Tipping 2001, Eqn 7)
  Phi <- Phi_full
  C <- sigma2 * diag(N) + tcrossprod(Phi, Phi * rep(1/alpha, each = nrow(Phi)))
  cholC <- chol(C)
  logdetC <- 2*sum(log(diag(cholC)))
  log_marginal_lik <- -0.5 * (N*log(2*pi) + logdetC + t(y) %*% solve(C, y))

  list(mu = mu, Sigma = Sigma, alpha = alpha, sigma2 = sigma2, keep_set=keep_set,
       log_marginal_lik = as.numeric(log_marginal_lik), iter=iter)
}


#' Posterior Predictive Sampling for RVM Objects
#'
#' Generates posterior predictive samples for a fitted \code{rvm} object, allowing for Bayesian model averaging over kernel lengthscales.
#'
#' @param object An object of class \code{rvm} as returned by \code{rvm()}.
#' @param newdata A matrix or data frame of new input locations at which to generate predictions. If \code{NULL}, predictions are generated for the training data.
#' @param samples Number of posterior predictive samples to draw (default is 1000).
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' For each posterior sample, a model (lengthscale) is sampled according to its posterior probability, and a Gaussian predictive draw is made using the corresponding RVM fit. If only a single model is present, all samples are drawn from that model.
#'
#' @return
#' A numeric matrix of dimension \code{samples x nrow(newdata)} containing posterior predictive draws.
#'
#' @seealso \code{\link{rvm}}
#'
#' @examples
#' # Assume fit is a fitted rvm object
#' Xtest <- matrix(runif(50), ncol = ncol(fit$X_train))
#' pred_draws <- predict(fit, Xtest, samples = 500)
#' # Posterior predictive mean
#' yhat <- colMeans(pred_draws)
#' @export
predict.rvm <- function(object, newdata = NULL, samples = 1000, ...) {
  if (is.null(newdata)){
    newdata <- object$X_train
  }
  n_models <- length(object$Fits)
  n_test <- nrow(newdata)
  preds <- matrix(NA, nrow = samples, ncol = n_test)

  # Posterior weights (normalize for safety)
  wts <- object$post_probs

  # Precompute all model predictive means and variances
  mu_post <- s2_post <- matrix(NA, nrow=n_models, ncol=n_test)
  for (j in seq_len(n_models)) {
    fit <- object$Fits[[j]]
    lscale <- object$lscale[j]

    # Define X's
    Xj <- newdata
    Xj_train <- object$X_train

    # Extract model information
    mu <- fit$mu
    Sigma <- fit$Sigma
    sigma2 <- fit$sigma2
    keep_set <- fit$keep_set

    # Make the K vectors
    K <- make_Phi(Xj, Xj_train, lscale=lscale)
    K <- K[,fit$selected] # Drop LASSO columns
    K <- K[,keep_set]     # Induce rvm sparsity

    # Get mu and sigma for each testing point
    mu_post[j,] <- as.vector(K %*% mu)
    s2_post[j,] <- sigma2 + rowSums((K %*% Sigma) * K)
  }

  # For single model, just use that one repeatedly
  if (n_models == 1){
    for (i in 1:samples) {
      preds[i, ] <- rnorm(n_test, mean = mu_post[1, ], sd = sqrt(s2_post[1, ]))
    }
  }else{
    # For multiple models, sample which model to use for each posterior sample
    mix_ids <- sample(seq_len(n_models), samples, replace = TRUE, prob = wts)
    for (i in 1:samples){
      j <- mix_ids[i]
      preds[i, ] <- rnorm(n_test, mean = mu_post[j, ], sd = sqrt(s2_post[j, ]))
    }
  }

  # Transform back to original scale if needed
  if (!is.null(object$y_scale) && !is.null(object$y_center)) {
    preds <- preds * object$y_scale + object$y_center
  }
  return(preds)
}

#' Plot Diagnostics for RVM Objects
#'
#' Plots diagnostic summaries for a fitted \code{rvm} model, including predicted vs observed and residual histogram.
#'
#' @param x An object of class \code{rvm} as returned by \code{rvm()}.
#' @param ... Additional plotting arguments (currently ignored).
#'
#' @details
#' By default, produces a two-panel plot: (1) a scatterplot of posterior predictive means vs observed \code{y}, with 95\% predictive intervals; (2) a histogram of residuals with a normal curve overlay.
#'
#' @return
#' Called for its side effect (plots). Invisibly returns \code{NULL}.
#'
#' @seealso \code{\link{predict.rvm}}
#'
#' @examples
#' # Assume fit is a fitted rvm object
#' plot(fit)
#' @export
plot.rvm <- function(x, ...){
  opar <- graphics::par(no.readonly = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), oma = c(0,
                                                              0, 0, 0))
  preds <- stats::predict(x)
  yhat <- apply(preds, 2, mean)
  yy <- x$y_center + x$y_scale * x$y_train
  ci <- 2 * apply(preds, 2, stats::sd)
  plot(yy, yhat, pch = 16, xlab = "y")
  graphics::segments(x0 = yy, y0 = yhat - ci, y1 = yhat +
                       ci, col = "orange")
  graphics::points(yy, yhat, pch = 16)
  graphics::abline(0, 1, col = "dodgerblue")
  rr <- yy - yhat
  graphics::hist(rr, breaks = ceiling(length(rr)^0.33 * diff(range(rr))/(3.5 *
                                                                           stats::sd(rr))), freq = F)
  #xx = seq(range(rr)[1], range(rr)[2], length.out = 100)
  graphics::curve(stats::dnorm(x, mean(rr), stats::sd(rr)),
                  add = TRUE, col = "orange", lwd=2)
  graphics::par(opar)
}

