#' Matching Pursuit GP
#'
#' A subset of data (SoD) approximate Gaussian process from Keerthi & Chu (2005)
#'
#' @param X A dataframe or matrix of predictors scaled to be between 0 and 1
#' @param y a reponse vector of length n
#' @param m Subset size for first loop. Successive iterations scale the subset size by \code{m_scale}.
#' @param cache_size Size of candidate cache (between m and n) for first loop. Successive iterations scaled by \code{m_scale}.
#' @param refresh_rate Fraction of the candidate set that is new during each iteration. Corresponds roughly to kappa from Keerthi & Chu paper.
#' @param sig2 Initial sigma squared value for Gaussian kernel
#' @param ell A vector (or scalar for isotropic kernel) of initial lengthscales.
#' @param loops Number of times do we loop between subset and fitting step (see details).
#' @param p_vech Cutoff for switching to Vecchia approx. E.g., for \code{p = ncol(X) <= p_vech} we use \code{GP_fit} package for GP fitting, and we use \code{GpGp} otherwise.
#' @param m_scale Multiplicative scaling factor for the subset size at each new iteration
#' @param verbose Logical.
#' @param visualize Generates a "movie" showing the selection process on the first two columns. Default of 0 indicates no visualization. Non-zero values indicates the frame rate in seconds (recommend 0.1). This is not sophisticated.
#' @details Algorithm has complexity O(nm^2). Candidate points are greedily selected to maximize a scoring criterion (Eq. 8 in Keerthi & Chu), conditional on kernel parameters. Once a subset is obtained, the GPfit package is used to estimate kernel parameters and the process repeats.
#' @references
#' Keerthi, Sathiya, and Wei Chu. "A matching pursuit approach to sparse gaussian process regression." Advances in neural information processing systems 18 (2005).
#'
#' Liu, Haitao, et al. "When Gaussian process meets big data: A review of scalable GPs." IEEE transactions on neural networks and learning systems 31.11 (2020): 4405-4423.
#'
#' Chalupka, Krzysztof, Christopher KI Williams, and Iain Murray. "A framework for evaluating approximation methods for Gaussian process regression." The Journal of Machine Learning Research 14.1 (2013): 333-350.
#'
#' Ranjan, P., Haynes, R., and Karsten, R. (2011). A Computationally Stable Approach to Gaussian Process Interpolation of Deterministic Computer Simulation Data, Technometrics, 53(4), 366 - 378.
#'
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- mpgp(X, y)
#' @export
mpgp <- function(X, y, m=NULL,
                      cache_size=100,
                      refresh_rate=0.59,
                      sig2=NULL, ell=NULL,
                      loops=2, p_vech=20, m_scale=1.5,
                      verbose=TRUE,
                      visualize=0, ...){
  if(var(y) == 0){
    if(verbose){
      warning("response is constant? Adding small noise for stability.")
    }
    y <- y + rnorm(y, 0, 1e-6)
  }
  if(is.null(sig2)){
    sig2 <- var(y)
  }
  if(is.null(ell)){
    ell <- cheap_lengthscale(X)
  }

  n <- length(y)
  if(is.null(m)){
    m <- max(50, floor(sqrt(n)))
  }
  m <- min(m, 1200)

  cache_size <- min(cache_size, m/2)
  kappa <- ceiling(cache_size * refresh_rate)

  if(verbose){
    cat("Loops remaining: ", loops, "\n",
        "\t Cache size: ", cache_size, "\n",
        "\t Subset size: ", m, "\n")
  }

  # Sample initial cache
  SoD_set <- NULL
  tmp_set <- cluster::pam(X,
                          k=cache_size,
                          pamonce=6,
                          keep.diss=FALSE, keep.data=FALSE)$id.med

  # Calculate the full kernel row K_i for each of the chosen i's (in tmp_set)
  K_cache <- matrix(NA, nrow=cache_size, ncol=n)
  for(j in seq_along(tmp_set)){
    K_cache[j,] <- kernel_function(X[tmp_set[j],,drop=FALSE], X)
  }

  # Store SoD K matrix
  K_SoD <- matrix(NA, nrow=m, ncol=n)

  cnt <- 0
  alpha_I <- NULL
  while(cnt < m){
    if(verbose){
      tmp <- rep(cnt, 3) == round(m * c(0.25, 0.5, 0.75))
      if(any(tmp)){
        cat("\t\t Progress:", c("25%", "50%", "75%")[tmp], "\n")
      }
    }

    # Evaluate criterion for each c
    Delta_crit <- alpha_star_i <- rep(NA, cache_size)
    for(j in seq_along(tmp_set)) {
      i <- tmp_set[j]
      K_i <- K_cache[j, ]
      if(cnt == 0) {
        K_I <- NULL
        tilde_ki <- numeric(0)
      } else {
        K_I <- K_SoD[1:cnt, , drop=FALSE]
        tilde_ki <- K_SoD[1:cnt, i, drop=TRUE]
      }
      K_ii <- 1 # K_i[i] # Or set as appropriate
      res <- delta_i(K_i, y, K_I, alpha_I, sig2, tilde_ki, K_ii)
      Delta_crit[j] <- res$Delta
      alpha_star_i[j] <- res$alpha_star
    }

    # Put the best one in the set
    cnt <- cnt + 1
    best_index <- which.max(Delta_crit)
    SoD_set <- c(SoD_set, tmp_set[best_index])
    K_SoD[cnt,] <- K_cache[best_index,,drop=FALSE]
    alpha_I <- c(alpha_I, alpha_star_i[best_index])

    # Update cache
    Delta_crit[best_index] <- -Inf
    ord <- order(Delta_crit)            # Lowest values first
    tmp_set <- tmp_set[ord]             # Re-order tmp_set
    K_cache <- K_cache[ord,,drop=FALSE] # Re-order rows of the cache

    # Sample new indices
    candidates <- setdiff(1:n, union(tmp_set, SoD_set))
    new_ind <- sample(candidates, kappa)
    tmp_set[1:kappa] <- new_ind

    # Calculate new K_i values
    for(j in seq_len(kappa)){
      K_cache[j,] <- kernel_function(X[tmp_set[j],,drop=FALSE], X)
    }

    if(visualize > 0){
      pts <- rep(1, n)
      col <- rep("grey", n)
      pts[tmp_set] <- 16
      col[tmp_set] <- "black"
      pts[SoD_set] <- 15
      col[SoD_set] <- "orange"
      plot(X[,1:2], pch=pts, col=col, main=paste0("Loops left:", loops))
      Sys.sleep(0.2)
    }
  }

  # Fit GP
  if(verbose){
    cat("\t\t Progress: 100%", "\n",
        "\t\t Fitting GP...", "\n\n")
  }

  X_sub <- X[SoD_set,,drop=FALSE]
  y_sub <- y[SoD_set]

  if(ncol(X) <= p_vech){
    fit <- GPfit::GP_fit(X_sub, y_sub,
                         corr=list(type="exponential", power=2),
                         ...)
    fit$ell <- 10^(-fit$beta/2) * 0.70711
    fit$Xfull <- X
    fit$yfull <- y
    fit$type <- "GPfit"
  }else{
    fit <- GpGp::fit_model(y_sub, X_sub, X=matrix(rep(1, nrow(X_sub)), ncol=1), silent=TRUE)
    fit$ell <- fit$covparms[2]
    fit$Xfull <- X
    fit$yfull <- y
    fit$type <- "GpGp"
  }


  # Either iterate or return
  if(loops == 1){
    class(fit) <- c("mpgp")
    return(fit)
  }else{
    return(mpgp(X, y, round(m * m_scale), round(m * cache_size), refresh_rate,
                               fit$sig2, fit$ell, loops-1, p_vech, m_scale,# Recursive stuff
                               verbose, visualize, ...))
  }
}

kernel_function <- function(x, X_mat, l = 0.1){
  if (is.null(dim(x)))                   # allow x as a plain vector
    x <- matrix(x, nrow = 1)

  d <- ncol(X_mat)
  if (length(l) == 1L)                   # isotropic ⇒ recycle scalar
    l <- rep(l, d)
  if (length(l) != d)
    stop("length(l) must be either 1 or ncol(X_mat)")

  scaled_diffs <- sweep(X_mat, 2, x)
  scaled_diffs <- sweep(scaled_diffs, 2, l, FUN = "/")
  distsq <- rowSums(scaled_diffs^2)

  return(exp(-0.5 * distsq))
}


delta_i <- function(K_i, y, K_I = NULL, alpha_I = NULL, sigma2 = 1, tilde_ki = NULL, K_ii = 1) {
  n <- length(y)

  if(length(alpha_I) == 0) {
    pred <- rep(0, n)
    tilde_term <- 0
  } else {
    pred <- as.numeric(crossprod(alpha_I, K_I))  # length-n prediction vector
    tilde_term <- sum(tilde_ki * alpha_I)
  }

  resid <- y - pred  # length-n

  numer <- sum(K_i * resid) - sigma2 * tilde_term
  denom <- sigma2 * K_ii + sum(K_i^2)

  alpha_star <- numer / denom
  Delta <- 0.5 * alpha_star^2 * denom
  list(Delta = Delta, alpha_star = alpha_star)
}

cheap_lengthscale <- function(X, frac = 0.05){
  n  <- nrow(X)
  d  <- ncol(X)
  r  <- min(500, max(30, round(frac * n)))
  idx <- sample(n, r)

  Xsub <- X[idx, , drop = FALSE]
  ell <- numeric(d)
  for (j in seq_len(d)) {
    v  <- Xsub[, j]
    diffs <- abs(outer(v, v, "-"))
    ell[j] <- median( diffs[upper.tri(diffs, diag = FALSE)])
  }
  ell[ell == 0] <- 1e-6
  return(ell)
}


#' Predict Method for class mpgp
#'
#' See \code{mpgp()} for details.
#'
#' @param object An object returned by the \code{mpgp} function.
#' @param newdata A dataframe of the same dimension as the training data.
#' @param samples How many posterior samples should be taken at each test point? If 0 or FALSE, then the MAP estimate is returned.
#' @param ... Additional arguments to predict
#' @details Predict function for mpgp. Just a wrapper for predict.GP but matches the output expected by \code{duqling} package.
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f)
#' fit <- mpgp(X, y)
#' predict(fit)
#'
#' @importFrom stats predict
#' @export
predict.mpgp <- function(object, newdata=NULL, samples=1000){
  if(is.null(newdata)){
    newdata <- object$Xfull
  }

  if(fit$type == "GPfit"){
    class(object) <- c("GP")
    out <- GPfit::predict.GP(object, newdata)
    preds <- matrix(NA, nrow=samples, ncol=nrow(newdata))
    for(i in 1:nrow(newdata)){
      preds[,i] <- out$Y_hat[i] + sqrt(out$MSE[i]) * rnorm(samples, 0, 1)
    }
  }else{
    preds <- GpGp::cond_sim(object, newdata, rep(1, nrow(newdata)), nsims=samples)
    preds <- t(preds)
  }
  return(preds)
}

#' Plot Method for class mpgp
#'
#' See \code{mpgp} for details.
#'
#' @param x An object returned by the \code{mpgp} function.
#' @param ... additional arguments passed to \code{plot}
#' @details Plot function for mpgp.
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- mpgp(X, y)
#' plot(fit)
#' @export
plot.mpgp <- function(x, ...){
  pred <- predict(x, x$Xfull, samples=1000)
  yhat <- colMeans(pred)
  plot(x$yfull, yhat, ...)
  graphics::abline(0, 1, lwd=2, col='orange')

  ci <- apply(pred, 2, function(yy) stats::quantile(yy, c(0.025, 0.975)))
  for(i in 1:ncol(ci)){
    graphics::segments(x$yfull[i], ci[1,i], x$yfull[i], ci[2,i])
  }
}

