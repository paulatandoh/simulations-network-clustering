source("config/settings.R"); source("R/model.R"); source("R/metrics.R"); source("R/experiment.R")
close_to <- function(a, b, tol = 1e-6) stopifnot(max(abs(a - b)) < tol)
settings <- simulation_settings()
stopifnot(identical(settings$S2$gamma, c(0.25, 0.75)), settings$S3$beta == 0.5,
          identical(settings$S4$theta, c(-0.35, 0.35)))
for (setting in settings) {
  setting$n <- 20L
  d <- simulate_network(setting, 180L)
  stopifnot(identical(d$Y, simulate_network(setting, 180L)$Y),
            all(d$Y == t(d$Y)), all(diag(d$Y) == 0), sum(d$x) == 10L)
}
d <- simulate_network(modifyList(settings$S1, list(n = 20L)), 122L)
a <- initialize_fit(d$Y, 123L)$alpha; dyads <- make_dyads(d$Y, d$x)
p <- c(-0.4, 0.35, 0.6); h <- 1e-5
for (use_cov in c(FALSE, TRUE)) {
  par <- if (use_cov) p else p[1:2]
  analytical <- network_objective(par, a, dyads, use_cov)
  g <- numeric(length(par)); H <- matrix(0, length(par), length(par))
  for (r in seq_along(par)) {
    inc <- rep(0, length(par)); inc[r] <- h
    plus <- network_objective(par + inc, a, dyads, use_cov)
    minus <- network_objective(par - inc, a, dyads, use_cov)
    g[r] <- (plus$value - minus$value) / (2 * h)
    H[, r] <- (plus$gradient - minus$gradient) / (2 * h)
  }
  close_to(g, analytical$gradient); close_to(H, analytical$hessian)
  gamma <- colMeans(a)
  # Independent exhaustive dyad calculation checks the complete ELBO/entropy sum.
  manual <- 0
  for (i in 1:(nrow(d$Y) - 1L)) for (j in (i + 1L):nrow(d$Y))
    for (k in 1:2) for (l in 1:2) {
      eta <- par[k] + par[l] + if (use_cov) par[3] * as.numeric(d$x[i] == d$x[j]) else 0
      manual <- manual + a[i, k] * a[j, l] * (d$Y[i, j] * eta - softplus(eta))
    }
  manual <- manual + sum(a * (matrix(log(gamma), nrow(a), 2, byrow = TRUE) - log(a)))
  close_to(manual, elbo(par, a, gamma, dyads, use_cov))
  updated <- e_step_mm(par, a, gamma, dyads, use_cov)
  stopifnot(all(updated >= 1e-8 * (1 - 1e-6)), max(abs(rowSums(updated) - 1)) < 1e-12,
            elbo(par, updated, gamma, dyads, use_cov) >= elbo(par, a, gamma, dyads, use_cov) - 1e-7)
  # Independently reconstruct the full surrogate: check tangency/minorization
  # and compare one closed-form node update against numerical optimization.
  surrogate <- function(candidate) {
    value <- 0
    for (u in seq_along(dyads$y)) {
      i <- dyads$i[u]; j <- dyads$j[u]
      for (k in 1:2) for (l in 1:2) {
        eta <- par[k] + par[l] + if (use_cov) par[3] * dyads$s[u] else 0
        lp <- dyads$y[u] * eta - softplus(eta)
        weight <- candidate[i,k]^2 * a[j,l] / (2*a[i,k]) +
                  candidate[j,l]^2 * a[i,k] / (2*a[j,l])
        value <- value + weight * lp
      }
    }
    value + sum(candidate * (matrix(log(gamma), nrow(a), 2, byrow=TRUE) -
                              log(a) - candidate/a + 1))
  }
  close_to(surrogate(a), elbo(par, a, gamma, dyads, use_cov))
  stopifnot(surrogate(updated) <= elbo(par, updated, gamma, dyads, use_cov) + 1e-7,
            surrogate(updated) >= surrogate(a) - 1e-7)
  objective <- function(first) {
    candidate <- a; candidate[1,] <- c(first,1-first); surrogate(candidate)
  }
  numerical <- optimize(objective, interval=c(1e-8,1-1e-8), maximum=TRUE, tol=1e-9)
  close_to(updated[1,1], numerical$maximum, tol=1e-5)
  fit <- fit_vem(d$Y, d$x, use_cov, initialize_fit(d$Y, 12L), max_outer = 100L)
  stopifnot(all(is.finite(fit$theta)), all(diff(fit$trace$elbo) >= -1e-7),
            all(fit$alpha > 0), abs(sum(fit$gamma) - 1) < 1e-10)
}
close_to(partition_scores(c(1,1,2,2), c(2,2,1,1)), c(RI=1, ARI=1))
close_to(partition_scores(c(1,1,2,2), c(1,2,1,2)), c(RI=1/3, ARI=-0.5))
fixture <- list(theta=c(0.5,-0.5), gamma=c(0.75,0.25),
                alpha=matrix(c(0.8,0.2,0.1,0.9), 2, byrow=TRUE), labels=c(1L,2L))
aligned <- align_fit(fixture, c(-0.5,0.5))
close_to(aligned$theta,c(-0.5,0.5)); close_to(aligned$gamma,c(0.25,0.75))
stopifnot(identical(aligned$labels, c(2L,1L)))
# All-failed groups must remain loggable, without causing a summary crash.
stopifnot(nrow(summarize_replicates(data.frame())) == 0L)
message("All checks passed: generation, derivatives, full ELBO, surrogate tangency/minorization, exact QP, MM ascent, positivity, label alignment, RI/ARI.")
