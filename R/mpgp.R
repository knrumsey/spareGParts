#' Matching Pursuit GP (new)
#'
#' A subset-of-data approximate Gaussian process using matching pursuit style
#' subset selection, with hetGP::mleHomGP as the primary fitting engine.
#'
#' @param X A dataframe or matrix of predictors scaled to be between 0 and 1
#' @param y A response vector of length n
#' @param m Subset size for first loop. Successive iterations scale the subset
#'   size by \code{m_scale}.
#' @param cache_size Size of candidate cache for first loop. Successive
#'   iterations are scaled by \code{m_scale}.
#' @param refresh_rate Fraction of the candidate set refreshed during each
#'   iteration.
#' @param sig2 Initial process variance value for subset scoring.
#' @param ell A vector (or scalar for isotropic kernel) of initial lengthscales.
#' @param loops Number of matching-pursuit / GP-refit loops.
#' @param p_vech Cutoff for switching to GpGp. For \code{p <= p_vech},
#'   \code{hetGP::mleHomGP} is used. Otherwise \code{GpGp} is used.
#' @param m_scale Multiplicative scaling factor for subset and cache size at
#'   each new iteration.
#' @param verbose Logical.
#' @param visualize Generates a crude visualization on the first two columns.
#' @param ... Additional arguments passed to the GP fitting backend.
#' @export
mpgp <- function(X, y, m = NULL,
                     cache_size = 100,
                     refresh_rate = 0.59,
                     sig2 = NULL, ell = NULL,
                     loops = 2, p_vech = 100, m_scale = 1.5,
                     verbose = TRUE,
                     visualize = 0, ...) {

  X <- as.matrix(X)
  y <- as.numeric(y)

  if (stats::var(y) < 1e-12) {
    if (verbose) {
      warning("response is nearly constant; adding small noise for stability.")
    }
    y <- y + stats::rnorm(length(y), 0, 1e-6)
  }

  if (is.null(sig2)) {
    sig2 <- stats::var(y)
  }
  if (is.null(ell)) {
    ell <- cheap_lengthscale(X)
  }

  n <- length(y)
  if (is.null(m)) {
    m <- max(50, floor(sqrt(n)))
  }
  m <- min(n - 1, m, 1200)

  cache_size <- max(m, min(n - 1, cache_size))
  kappa <- ceiling(cache_size * refresh_rate)

  if (verbose) {
    cat("Loops remaining: ", loops, "\n",
        "\t Cache size: ", cache_size, "\n",
        "\t Subset size: ", m, "\n")
  }

  # Initial candidate cache via PAM medoids
  SoD_set <- integer(0)
  tmp_set <- cluster::pam(
    X,
    k = cache_size,
    pamonce = 6,
    keep.diss = FALSE,
    keep.data = FALSE
  )$id.med

  # Cache full kernel rows for candidate points
  K_cache <- matrix(NA_real_, nrow = cache_size, ncol = n)
  for (j in seq_along(tmp_set)) {
    K_cache[j, ] <- kernel_function(X[tmp_set[j], , drop = FALSE], X, l = ell)
  }

  # Store selected kernel rows
  K_SoD <- matrix(NA_real_, nrow = m, ncol = n)

  cnt <- 0
  alpha_I <- numeric(0)

  while (cnt < m) {
    if (verbose) {
      tmp <- rep(cnt, 3) == round(m * c(0.25, 0.5, 0.75))
      if (any(tmp)) {
        cat("\t\t Progress:", c("25%", "50%", "75%")[tmp], "\n")
      }
    }

    Delta_crit <- alpha_star_i <- rep(NA_real_, cache_size)

    for (j in seq_along(tmp_set)) {
      i <- tmp_set[j]
      K_i <- K_cache[j, ]

      if (cnt == 0) {
        K_I <- NULL
        tilde_ki <- numeric(0)
      } else {
        K_I <- K_SoD[1:cnt, , drop = FALSE]
        tilde_ki <- K_SoD[1:cnt, i, drop = TRUE]
      }

      K_ii <- 1
      res <- delta_i(K_i, y, K_I, alpha_I, sig2, tilde_ki, K_ii)
      Delta_crit[j] <- res$Delta
      alpha_star_i[j] <- res$alpha_star
    }

    # Add best point to subset
    cnt <- cnt + 1
    best_index <- which.max(Delta_crit)
    SoD_set <- c(SoD_set, tmp_set[best_index])
    K_SoD[cnt, ] <- K_cache[best_index, , drop = FALSE]
    alpha_I <- c(alpha_I, alpha_star_i[best_index])

    # Reorder cache so weakest candidates come first
    Delta_crit[best_index] <- -Inf
    ord <- order(Delta_crit)
    tmp_set <- tmp_set[ord]
    K_cache <- K_cache[ord, , drop = FALSE]

    # Refresh up to kappa cache points, but near the end there may be fewer
    # unused candidates available than kappa
    candidates <- setdiff(seq_len(n), union(tmp_set, SoD_set))
    kappa2 <- min(kappa, length(candidates))

    if (kappa2 > 0) {
      new_ind <- sample(candidates, kappa2)
      tmp_set[seq_len(kappa2)] <- new_ind

      for (j in seq_len(kappa2)) {
        K_cache[j, ] <- kernel_function(
          X[tmp_set[j], , drop = FALSE],
          X,
          l = ell
        )
      }
    }

    if (visualize > 0) {
      pts <- rep(1, n)
      col <- rep("grey", n)
      pts[tmp_set] <- 16
      col[tmp_set] <- "black"
      pts[SoD_set] <- 15
      col[SoD_set] <- "orange"
      plot(X[, 1:2], pch = pts, col = col, main = paste0("Loops left: ", loops))
      Sys.sleep(visualize)
    }
  }

  if (verbose) {
    cat("\t\t Progress: 100%\n",
        "\t\t Fitting GP...\n\n")
  }

  X_sub <- X[SoD_set, , drop = FALSE]
  y_sub <- y[SoD_set]

  if (ncol(X) <= p_vech) {
    fit <- hetGP::mleHomGP(
      X = X_sub,
      Z = y_sub,
      covtype = "Gaussian",
      ...
    )

    fit$backend <- "homGP"

    # IMPORTANT:
    # This conversion assumes hetGP Gaussian correlation is of the form
    # exp( - sum((x - x')^2 / theta) ), while kernel_function uses
    # exp( -0.5 * sum((x - x')^2 / ell^2) ).
    # Under that convention, ell = sqrt(theta / 2).
    fit$ell_mp <- sqrt(fit$theta / 2)

    # Process variance and nugget variance on the response scale
    fit$sig2_mp <- fit$nu_hat
    fit$nugget_mp <- fit$nu_hat * fit$g

    fit$Xfull <- X
    fit$yfull <- y

  } else {
    fit <- GpGp::fit_model(
      y_sub,
      X_sub,
      X = matrix(1, nrow = nrow(X_sub), ncol = 1),
      covfun_name = "exponential_isotropic",
      silent = TRUE
    )

    fit$backend <- "GpGp"
    fit$ell_mp <- sqrt(fit$covparms[2] / 2)
    fit$sig2_mp <- fit$covparms[1]
    fit$nugget_mp <- if (length(fit$covparms) >= 3) {
      fit$covparms[1] * fit$covparms[3]
    } else {
      0
    }

    fit$Xfull <- X
    fit$yfull <- y
  }

  if (loops == 1) {
    class(fit) <- c("mpgp", class(fit))
    return(fit)
  } else {
    return(
      mpgp(
        X = X,
        y = y,
        m = round(m * m_scale),
        cache_size = round(cache_size * m_scale),
        refresh_rate = refresh_rate,
        sig2 = fit$sig2_mp,
        ell = fit$ell_mp,
        loops = loops - 1,
        p_vech = p_vech,
        m_scale = m_scale,
        verbose = verbose,
        visualize = visualize,
        ...
      )
    )
  }
}


kernel_function <- function(x, X_mat, l = 0.1) {
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1)
  }

  d <- ncol(X_mat)
  if (length(l) == 1L) {
    l <- rep(l, d)
  }
  if (length(l) != d) {
    stop("length(l) must be either 1 or ncol(X_mat)")
  }

  scaled_diffs <- sweep(X_mat, 2, x)
  scaled_diffs <- sweep(scaled_diffs, 2, l, FUN = "/")
  distsq <- rowSums(scaled_diffs^2)

  exp(-0.5 * distsq)
}


delta_i <- function(K_i, y, K_I = NULL, alpha_I = NULL,
                    sigma2 = 1, tilde_ki = NULL, K_ii = 1) {

  if (length(alpha_I) == 0) {
    pred <- rep(0, length(y))
    tilde_term <- 0
  } else {
    pred <- as.numeric(crossprod(alpha_I, K_I))
    tilde_term <- sum(tilde_ki * alpha_I)
  }

  resid <- y - pred
  numer <- sum(K_i * resid) - sigma2 * tilde_term
  denom <- sigma2 * K_ii + sum(K_i^2)

  alpha_star <- numer / denom
  Delta <- 0.5 * alpha_star^2 * denom

  list(Delta = Delta, alpha_star = alpha_star)
}


cheap_lengthscale <- function(X, frac = 0.05) {
  n <- nrow(X)
  d <- ncol(X)
  r <- min(500, max(30, round(frac * n)))
  idx <- sample(n, r)

  Xsub <- X[idx, , drop = FALSE]
  ell <- numeric(d)

  for (j in seq_len(d)) {
    v <- Xsub[, j]
    diffs <- abs(outer(v, v, "-"))
    ell[j] <- median(diffs[upper.tri(diffs, diag = FALSE)])
  }

  ell[ell == 0] <- 1e-6
  ell
}


#' Predict method for class mpgp
#'
#' @param object An object returned by \code{mpgp()}.
#' @param newdata A matrix/dataframe of predictors.
#' @param samples Number of predictive draws to return. If 0 or FALSE, returns
#'   a 1 x n matrix of predictive means.
#' @param nugget Logical. If TRUE, include nugget / observation noise in the
#'   predictive variance when available.
#' @param ... Ignored.
#' @export
predict.mpgp <- function(object, newdata = NULL, samples = 1000,
                             nugget = TRUE, ...) {

  if (is.null(newdata)) {
    newdata <- object$Xfull
  }
  newdata <- as.matrix(newdata)

  return_mean_only <- identical(samples, 0) || identical(samples, FALSE)

  if (object$backend == "homGP") {
    obj <- object
    class(obj) <- "homGP"
    out <- stats::predict(obj, x = newdata)

    mu <- as.numeric(out$mean)
    var_pred <- as.numeric(out$sd2)

    if (nugget) {
      if (!is.null(out$nugs)) {
        var_pred <- var_pred + as.numeric(out$nugs)
      } else if (!is.null(object$nugget_mp)) {
        var_pred <- var_pred + object$nugget_mp
      }
    }

    var_pred <- pmax(var_pred, 0)

    if (return_mean_only) {
      return(matrix(mu, nrow = 1))
    }

    preds <- matrix(NA_real_, nrow = samples, ncol = nrow(newdata))
    for (i in seq_len(nrow(newdata))) {
      preds[, i] <- mu[i] + sqrt(var_pred[i]) * stats::rnorm(samples)
    }
    return(preds)
  }

  if (object$backend == "GpGp") {
    if (return_mean_only) {
      mu <- tryCatch(
        {
          as.numeric(
            GpGp::predictions(
              object,
              newdata,
              matrix(1, nrow = nrow(newdata), ncol = 1)
            )$mean
          )
        },
        error = function(e) {
          colMeans(
            t(
              GpGp::cond_sim(
                object,
                newdata,
                matrix(1, nrow = nrow(newdata), ncol = 1),
                nsims = 200
              )
            )
          )
        }
      )
      return(matrix(mu, nrow = 1))
    }

    preds <- GpGp::cond_sim(
      object,
      newdata,
      matrix(1, nrow = nrow(newdata), ncol = 1),
      nsims = samples
    )
    preds <- t(preds)

    if (nugget && !is.null(object$nugget_mp) &&
        is.finite(object$nugget_mp) && object$nugget_mp > 0) {
      preds <- preds + matrix(
        stats::rnorm(length(preds), sd = sqrt(object$nugget_mp)),
        nrow = nrow(preds),
        ncol = ncol(preds)
      )
    }

    return(preds)
  }

  stop("Unknown backend in predict.mpgp().")
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
