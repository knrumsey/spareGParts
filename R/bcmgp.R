#' Bayesian Committee Machine (GP)
#'
#' Many subset of data GPs are aggregated via the Robust Bayesian Committee algorithm of
#'
#' @param X A dataframe or matrix of predictors scaled to be between 0 and 1
#' @param y a reponse vector of length n
#' @param M Number of data partitions (defaults to \code{floor(sqrt(n))}).
#' @param max_size The maximum partition size. Ignored when \code{partition="random"}.
#' @param partition Data partition type; Either "random" or "cluster"; the latter uses a fast PAM algorithm to form the data partitions.
#' @param expert_weight How should the experts be weighted? Setting to \code{"uniform"} gives beta_i = 1 (as in Tresp 2000) while the default \code{"varying"} sets beta_i as formulated in Deisenroth and Ng 2015.
#' @param isotropic Logical. Should each input have the same lengthscale? (Defaults to \code{FALSE} when \code{ncol(X) <= 5} and \code{TRUE} otherwise).
#' @details Algorithm has complexity O(n(n/M)^2).
#' @references
#' Tresp, Volker. "A Bayesian committee machine." Neural computation 12.11 (2000): 2719-2741.
#'
#' Deisenroth, Marc, and Jun Wei Ng. "Distributed gaussian processes." International conference on machine learning. PMLR, 2015.
#'
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- bcmgp(X, y)
#' @export
bcmgp <- function(X, y,
                  M=NULL, max_size=NULL,
                  partition="cluster",
                  expert_weight="varying",
                  isotropic=FALSE){
  if(max(X) > 1 | min(X) < 0) warning("Expecting inputs on (0, 1) scale")
  mu_y <- mean(y)
  sigma_y <- sd(y)
  y_std <- (y - mean(y)) / sd(y)

  n <- nrow(X)
  p <- ncol(X)
  if(is.null(M)) M <- floor(sqrt(n)/2)
  if(is.null(max_size)) max_size <- ceiling(n/M)
  if(max_size * M  < n){
    warning("max_size is too small.")
    max_size <- ceiling(n/M)
  }
  if(partition == "cluster"){
    clust <- modified_pam(X, M, max_size)
  }else{
    clust <- sample(rep(1:M, ceiling(n/M)))[1:n]
  }
  if(is.null(isotropic)){
    if(p <= 5){
      isotropic <- FALSE
    }else{
      isotropic <- TRUE
    }
  }

  X_list <- y_list <- list()
  for(m in seq_len(M)){
    ind <- which(clust == m)
    X_list[[m]] <- X[ind,,drop=FALSE]
    y_list[[m]] <- y_std[ind]
  }
  params <- fit_gp_shared_hypers(X_list, y_list,
                                 isotropic=isotropic, upper_ell=2)

  fit <- list(X_list=X_list, y_list=y_list, params=params, clust=clust, weight=expert_weight, mu_y=mu_y, sigma_y=sigma_y)
  class(fit) <- "bcmgp"
  return(fit)
}

#' Predict Method for class bcmgp
#'
#' See \code{bcmgp()} for details.
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
#' fit <- bcmgp(X, y)
#' predict(fit)
#'
#' @importFrom stats predict
#' @export
predict.bcmgp <- function(object, newdata=NULL, samples=1000){
  if(is.null(newdata)){
    if(length(object$X_list) > 1){
      newdata <- do.call(rbind, object$X_list)
    }else{
      newdata <- object$X_list[[1]]
    }

  }

  M <- length(object$X_list)
  sigma_prior <- sqrt(object$params$sigma^2 + object$params$tau^2)
  preds <- matrix(NA, nrow=samples, ncol=nrow(newdata))
  for(i in 1:nrow(newdata)){
    xx <- newdata[i,,drop=FALSE]
    mu <- sigma <- rep(NA, M)
    beta <- rep(1, M)
    for(m in 1:M){
      gp_params <- gp_predict(xx,
                              object$X_list[[m]], object$y_list[[m]],
                              object$params$ell,
                              object$params$sigma^2,
                              object$params$tau^2)

      mu[m] <- gp_params$mean
      sigma[m] <- sqrt(gp_params$variance)
      if(object$weight == "varying"){
        beta[m] <- log(sigma_prior) - log(sigma[m])
      }
    }

    sigma2_post <- (sum(beta * sigma^(-2)) + (1 - sum(beta)) * sigma_prior^(-2))^(-1)
    mean_post <- sigma2_post * sum(beta * sigma^(-2) * mu)
    preds[,i] <- rnorm(samples, mean_post, sqrt(abs(sigma2_post)))
  }

  preds <- object$mu_y + object$sigma_y * preds
  return(preds)
}

#' Plot Method for class bcmgp
#'
#' See \code{bcmgp} for details.
#'
#' @param x An object returned by the \code{bcmgp} function.
#' @param ... additional arguments passed to \code{plot}
#' @details Plot function for bcmgp.
#' @examples
#' X <- lhs::maximinLHS(100, 2)
#' f <- function(x) 10.391*((x[1]-0.4)*(x[2]-0.6) + 0.36)
#' y <- apply(X, 1, f) + stats::rnorm(100, 0, 0.1)
#' fit <- bcmgp(X, y)
#' plot(fit)
#' @export
plot.bcmgp <- function(x, ...){
  yfull <- x$mu_y + x$sigma_y * unlist(x$y_list)
  pred <- predict(x, samples=1000)
  yhat <- colMeans(pred)
  plot(yfull, yhat, ...)
  graphics::abline(0, 1, lwd=2, col='orange')

  ci <- apply(pred, 2, function(yy) stats::quantile(yy, c(0.025, 0.975)))
  for(i in 1:ncol(ci)){
    graphics::segments(yfull[i], ci[1,i], yfull[i], ci[2,i])
  }
}


modified_pam <- function(X, M, max_size){
  max_size <- max(max_size, ceiling(nrow(X) / M))

  clust <- cluster::pam(X, k=M, pamonce=6)
  medioids <- clust$medoids
  tab <- table(clust$clustering)
  flag <- max(tab) > max_size
  while(flag){
    clust_indx <- which.max(tab)[1]
    samp_indx <- tab[clust_indx]
    point_indx <- which(clust$clustering == clust_indx)[samp_indx]
    point <- X[point_indx,]

    invalid_clusts <- which(tab >= max_size)
    medioids[invalid_clusts,] <- NA #rep(1e9, length(invalid_clusts)*ncol(X))
    new_clust <- which.min(apply(medioids, 1, function(xx) mean((xx-point)^2)))

    # Edit cluster here
    clust$clustering[point_indx] <- new_clust

    tab <- table(clust$clustering)
    flag <- max(tab) > max_size
  }
  return(clust$clustering)
}


fit_gp_shared_hypers <- function(X_list, y_list,
                                 isotropic = FALSE,
                                 lower_ell = 1e-2,
                                 upper_ell = 3,
                                 lower_sigma = 1e-6,
                                 lower_tau   = 1e-4,
                                 jitter = 1e-8,
                                 maxit = 200,
                                 verbose = FALSE)
{
  if (!is.list(X_list)) X_list <- list(X_list)
  if (!is.list(y_list)) y_list <- list(y_list)

  M <- length(X_list);  d <- ncol(X_list[[1L]])
  stopifnot(all(vapply(X_list, ncol, 1L) == d))

  ## ---- initial values --------------------------------------------------- ##
  y_all <- unlist(y_list, use.names = FALSE)
  theta0 <- c(rep(0.5, if (isotropic) 1 else d),  # ell(·)
              sd(y_all),                          # sigma
              0.1 * sd(y_all) + 1e-6)             # tau

  ## ---- bounds for L-BFGS-B --------------------------------------------- ##
  lower <- log(c( rep(lower_ell, if (isotropic) 1 else d),
                  lower_sigma,
                  lower_tau))
  upper <- log(c( rep(upper_ell, if (isotropic) 1 else d),
                  Inf,
                  Inf))


  ## ---- NLL for one expert ---------------------------------------------- ##
  nll_one <- function(X, y, ell, sigma2, tau2)
  {
    K <- kfun(X, X, ell, sigma2);  diag(K) <- diag(K) + tau2 + jitter
    L <- tryCatch(chol(K), error = function(e) NULL)
    if (is.null(L)) return(1e14)                       # large finite penalty
    alpha <- backsolve(L, forwardsolve(L, y, trans = TRUE))
    0.5 * sum(y * alpha) + sum(log(diag(L))) + 0.5 * length(y) * log(2*pi)
  }

  ## ---- total NLL -------------------------------------------------------- ##
  nll_total <- function(theta_log)
  {
    ell    <- exp(theta_log[ seq_len(if (isotropic) 1 else d) ])
    sigma2 <- exp(2*theta_log[ length(ell) + 1L ])
    tau2   <- exp(2*theta_log[ length(ell) + 2L ])

    val <- 0
    for (m in seq_len(M))
      val <- val + nll_one(X_list[[m]], y_list[[m]], ell, sigma2, tau2)
    if (verbose) cat("NLL:", format(val, digits = 6), " ell:", ell,
                     " sig:", sqrt(sigma2), " tau:", sqrt(tau2), "\n")
    val
  }

  ## ---- optimisation ----------------------------------------------------- ##
  opt <- optim(par = log(theta0), fn = nll_total,
               method = "L-BFGS-B", lower = lower, upper = upper,
               control = list(maxit = maxit))

  ## ---- unpack (repeat ell if isotropic) --------------------------------- ##
  ell_est <- if (isotropic) rep( exp(opt$par[1L]), d )
  else            exp( opt$par[ seq_len(d) ] )

  list(ell   = ell_est,
       sigma = exp(opt$par[length(ell_est) + 1L]),
       tau   = exp(opt$par[length(ell_est) + 2L]),
       value = opt$value,
       convergence = opt$convergence,
       counts = opt$counts,
       message = opt$message)
}


kfun <- function(X1, X2, ell, sigma2){
  if (length(ell) == 1L) ell <- rep(ell, ncol(X1))          # isotropic
  Z1 <- sweep(X1, 2L, ell, FUN = "/")
  Z2 <- sweep(X2, 2L, ell, FUN = "/")
  r2 <- outer(rowSums(Z1^2), rowSums(Z2^2), "+") - 2 * tcrossprod(Z1, Z2)
  sigma2 * exp(-0.5 * pmax(r2, 0))
}

gp_predict <- function(xstar, X_train, y_train, ell, sigma2, tau2){
  # Use your kernel
  K <- kfun(X_train, X_train, ell, sigma2) + diag(tau2, nrow(X_train))
  Ks <- kfun(X_train, matrix(xstar, nrow=1), ell, sigma2) # n x 1
  Kss <- sigma2 # rbf kernel at same point is sigma2

  # Numerically stable solve
  # Attempt to solve with increasing jitter values
  jitter_values <- c(0, 1e-9, 1e-6, 1e-3)
  alpha <- NULL
  success <- FALSE

  for (j in jitter_values) {
    Kj <- K
    if (j > 0) {
      Kj <- Kj + diag(j, nrow(Kj))
    }

    attempt <- try(solve(Kj, b), silent = TRUE)  # Replace 'b' with your RHS
    if (!inherits(attempt, "try-error") && all(is.finite(attempt))) {
      alpha <- attempt
      success <- TRUE
      break
    }
  }
  mu <- as.numeric(t(Ks) %*% alpha)
  v <- solve(K, Ks)     # (n x 1)
  var <- as.numeric(Kss - t(Ks) %*% v + tau2)
  var <- max(var, 1e-10) # prevent negative variance
  return(list(mean = mu, variance = var))
}
