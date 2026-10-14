# Undirected binary network with additive latent-group effects and binary similarity.
softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

validate_network <- function(Y, x) {
  if (!is.matrix(Y) || nrow(Y) != ncol(Y) || nrow(Y) < 4L ||
      anyNA(Y) || any(!Y %in% c(0, 1)) || !isTRUE(all.equal(Y, t(Y))) ||
      any(diag(Y) != 0)) stop("Y must be a symmetric binary matrix with zero diagonal.")
  if (length(x) != nrow(Y) || anyNA(x) || any(!x %in% c(0, 1)))
    stop("x must contain one binary covariate per node.")
}

make_dyads <- function(Y, x) {
  validate_network(Y, x)
  ind <- which(upper.tri(Y), arr.ind = TRUE)
  list(n = nrow(Y), i = ind[, 1L], j = ind[, 2L], y = Y[ind],
       s = as.numeric(x[ind[, 1L]] == x[ind[, 2L]]))
}

simulate_network <- function(setting, seed) {
  if (length(setting$theta) != 2L || length(setting$gamma) != 2L ||
      any(setting$gamma <= 0) || abs(sum(setting$gamma) - 1) > 1e-10)
    stop("This implementation requires two groups and positive proportions summing to one.")
  n <- as.integer(setting$n)
  if (n < 4L || !all(is.finite(c(setting$theta, setting$beta)))) stop("Invalid setting.")
  set.seed(seed)
  # iid memberships follow the manuscript's mixture model. Group counts are random.
  z <- sample.int(2L, n, replace = TRUE, prob = setting$gamma)
  # Exactly balanced binary categories for even n; randomized independently of z.
  x <- sample(rep(c(0L, 1L), length.out = n))
  ind <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  eta <- setting$theta[z[ind[, 1L]]] + setting$theta[z[ind[, 2L]]] +
    setting$beta * as.numeric(x[ind[, 1L]] == x[ind[, 2L]])
  Y <- matrix(0L, n, n)
  Y[ind] <- rbinom(nrow(ind), 1L, plogis(eta))
  Y <- Y + t(Y)
  list(Y = Y, x = x, z = z, seed = seed, setting = setting)
}

# Network component of the ELBO and its FULL gradient/Hessian, with alpha fixed.
network_objective <- function(par, alpha, dyads, use_covariate = TRUE,
                              derivatives = TRUE) {
  p <- if (use_covariate) 3L else 2L
  if (length(par) != p) stop("Parameter dimension does not match the model.")
  value <- 0; gradient <- numeric(p); hessian <- matrix(0, p, p)
  for (k in 1:2) for (l in 1:2) {
    design <- cbind(rep(as.numeric(k == 1L) + as.numeric(l == 1L), length(dyads$y)),
                    rep(as.numeric(k == 2L) + as.numeric(l == 2L), length(dyads$y)))
    if (use_covariate) design <- cbind(design, dyads$s)
    eta <- drop(design %*% par)
    w <- alpha[dyads$i, k] * alpha[dyads$j, l]
    value <- value + sum(w * (dyads$y * eta - softplus(eta)))
    if (derivatives) {
      prob <- plogis(eta)
      gradient <- gradient + drop(crossprod(design, w * (dyads$y - prob)))
      hessian <- hessian - crossprod(design, design * as.numeric(w * prob * (1 - prob)))
    }
  }
  list(value = value, gradient = gradient, hessian = hessian)
}

elbo <- function(par, alpha, gamma, dyads, use_covariate = TRUE) {
  if (any(alpha <= 0) || any(gamma <= 0)) stop("alpha and gamma must be strictly positive.")
  network_objective(par, alpha, dyads, use_covariate, FALSE)$value +
    sum(alpha * (matrix(log(gamma), nrow(alpha), 2L, byrow = TRUE) - log(alpha)))
}

# Separable MM E-step. For K=2 the simplex-constrained diagonal QP has an
# exact closed form, so quadprog and the original n*K by n*K QP are unnecessary.
# It maximizes the same surrogate, on alpha_ik >= eps and sum_k alpha_ik = 1.
e_step_mm <- function(par, alpha, gamma, dyads, use_covariate = TRUE, eps = 1e-8) {
  if (eps <= 0 || eps >= 0.5 || any(alpha < eps * (1 - 1e-6))) stop("Invalid positivity floor.")
  n <- dyads$n
  coefficient <- matrix(0, n, 2L)
  for (k in 1:2) for (l in 1:2) {
    eta <- par[k] + par[l] + if (use_covariate) par[3L] * dyads$s else 0
    lp <- dyads$y * eta - softplus(eta)
    log_matrix <- matrix(0, n, n)
    log_matrix[cbind(dyads$i, dyads$j)] <- lp
    log_matrix[cbind(dyads$j, dyads$i)] <- lp
    coefficient[, k] <- coefficient[, k] + drop(log_matrix %*% alpha[, l]) / 2
  }
  # S = sum_ik c_ik * alpha_ik^2 + d_ik * alpha_ik.
  quadratic <- (coefficient - 1) / alpha
  D <- -2 * quadratic
  d <- matrix(log(gamma), n, 2L, byrow = TRUE) - log(alpha) + 1
  first <- (d[, 1L] - d[, 2L] + D[, 2L]) / (D[, 1L] + D[, 2L])
  first <- pmin(pmax(first, eps), 1 - eps)
  cbind(first, 1 - first)
}

m_step_newton <- function(par, alpha, dyads, use_covariate = TRUE,
                          max_iter = 500L, tol = 1e-4) {
  backtracks <- 0L; regularizations <- 0L
  for (iteration in seq_len(max_iter)) {
    old <- network_objective(par, alpha, dyads, use_covariate)
    if (max(abs(old$gradient)) < 1e-7)
      return(list(par = par, converged = TRUE, iterations = iteration - 1L,
                  backtracks = backtracks, regularizations = regularizations))
    information <- -old$hessian
    step <- NULL
    for (ridge in c(0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2)) {
      candidate_step <- tryCatch(solve(information + diag(ridge, length(par)), old$gradient),
                                 error = function(e) NULL)
      if (!is.null(candidate_step) && all(is.finite(candidate_step)) &&
          sum(candidate_step * old$gradient) > 0) {
        step <- candidate_step
        regularizations <- regularizations + as.integer(ridge > 0)
        break
      }
    }
    if (is.null(step)) return(list(par = par, converged = FALSE, iterations = iteration,
                                   backtracks = backtracks, regularizations = regularizations))
    accepted <- FALSE
    scale <- 1
    for (attempt in 0:30) {
      proposed <- par + scale * step
      if (all(is.finite(proposed))) {
        value <- network_objective(proposed, alpha, dyads, use_covariate, FALSE)$value
        if (is.finite(value) && value >= old$value) { accepted <- TRUE; break }
      }
      scale <- scale / 2; backtracks <- backtracks + 1L
    }
    if (!accepted) return(list(par = par, converged = FALSE, iterations = iteration,
                               backtracks = backtracks, regularizations = regularizations))
    change <- sqrt(sum((proposed - par)^2)); par <- proposed
    if (change < tol) return(list(par = par, converged = TRUE, iterations = iteration,
                                  backtracks = backtracks, regularizations = regularizations))
  }
  list(par = par, converged = FALSE, iterations = max_iter,
       backtracks = backtracks, regularizations = regularizations)
}

initialize_fit <- function(Y, seed, method = "degree", eps = 1e-8) {
  set.seed(seed); n <- nrow(Y)
  if (method == "random") {
    a <- runif(n, 0.1, 0.9)
    alpha <- cbind(a, 1 - a)
  } else if (method == "degree") {
    degree <- rowSums(Y) + runif(n, -0.1, 0.1)
    labels <- kmeans(degree, centers = 2L, nstart = 1L)$cluster
    if (mean(degree[labels == 1L]) > mean(degree[labels == 2L])) labels <- 3L - labels
    alpha <- matrix(0.05, n, 2L); alpha[cbind(seq_len(n), labels)] <- 0.95
    # Random perturbations provide distinct initializations without using true labels.
    a <- pmin(pmax(alpha[, 1L] + runif(n, -0.025, 0.025), eps), 1 - eps)
    alpha <- cbind(a, 1 - a)
  } else if (method == "mmsb") {
    if (!requireNamespace("lda", quietly = TRUE)) stop("Install lda to use initialization='mmsb'.")
    mmsb <- lda::mmsb.collapsed.gibbs.sampler(Y, 2L, num.iterations = 1000L,
                                            alpha = 0.5, burnin = 20L,
                                            beta.prior = list(1, diag(7, 2L) + 1))
    counts <- t(mmsb$document_sums)
    a <- counts[, 1L] / rowSums(counts)
    a[!is.finite(a)] <- 0.5
    a <- pmin(pmax(a, eps), 1 - eps); alpha <- cbind(a, 1 - a)
  } else stop("Unknown initialization method.")
  density <- sum(Y) / (n * (n - 1))
  center <- qlogis(pmin(pmax(density, 0.01), 0.99)) / 2
  list(alpha = alpha, theta = center + c(-0.3, 0.3), beta = 0)
}

fit_vem <- function(Y, x, use_covariate = TRUE, initial,
                    max_outer = 500L, max_inner = 500L,
                    tol_outer = 1e-4, tol_inner = 1e-4, alpha_floor = 1e-8) {
  started <- proc.time()[[3L]]
  dyads <- make_dyads(Y, x)
  alpha <- initial$alpha
  if (!all(dim(alpha) == c(nrow(Y), 2L)) || any(!is.finite(alpha)) ||
      any(alpha < alpha_floor) || max(abs(rowSums(alpha) - 1)) > 1e-8)
    stop("Invalid initial membership probabilities.")
  gamma <- colMeans(alpha)
  par <- c(initial$theta, if (use_covariate) initial$beta else NULL)
  old <- elbo(par, alpha, gamma, dyads, use_covariate)
  trace <- data.frame(iteration = 0L, elbo = old)
  converged <- FALSE; reason <- "outer_iteration_limit"
  inner_iterations <- 0L; all_inner_converged <- TRUE
  backtracks <- 0L; regularizations <- 0L
  for (iteration in seq_len(max_outer)) {
    proposed_alpha <- e_step_mm(par, alpha, gamma, dyads, use_covariate, alpha_floor)
    proposed_gamma <- colMeans(proposed_alpha)
    m <- m_step_newton(par, proposed_alpha, dyads, use_covariate, max_inner, tol_inner)
    value <- elbo(m$par, proposed_alpha, proposed_gamma, dyads, use_covariate)
    if (!is.finite(value)) stop("Nonfinite ELBO.")
    if (value < old - 1e-8 * (1 + abs(old))) stop("ELBO decreased beyond numerical tolerance.")
    alpha <- proposed_alpha; gamma <- proposed_gamma; par <- m$par
    inner_iterations <- inner_iterations + m$iterations
    all_inner_converged <- all_inner_converged && m$converged
    backtracks <- backtracks + m$backtracks
    regularizations <- regularizations + m$regularizations
    trace <- rbind(trace, data.frame(iteration = iteration, elbo = value))
    relative_change <- abs(value - old) / max(1, abs(old))
    old <- value
    if (relative_change < tol_outer) {
      converged <- m$converged
      reason <- if (converged) "elbo_tolerance" else "inner_optimizer_not_converged"
      break
    }
  }
  list(theta = par[1:2], beta = if (use_covariate) par[3L] else NA_real_,
       gamma = gamma, alpha = alpha, labels = max.col(alpha, ties.method = "first"),
       elbo = old, converged = converged, reason = reason,
       all_inner_converged = all_inner_converged,
       outer_iterations = iteration, inner_iterations = inner_iterations,
       backtracks = backtracks, regularizations = regularizations, trace = trace,
       runtime_s = unname(proc.time()[[3L]] - started))
}
