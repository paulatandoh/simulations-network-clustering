# Label-invariant partition scores without fossil/mclust dependencies.
choose2 <- function(x) x * (x - 1) / 2
partition_scores <- function(truth, estimate) {
  if (length(truth) != length(estimate) || length(truth) < 2L || anyNA(truth) || anyNA(estimate))
    stop("Partitions must have equal lengths, at least two nodes, and no missing values.")
  tab <- table(truth, estimate)
  both <- sum(choose2(tab)); rows <- sum(choose2(rowSums(tab)))
  cols <- sum(choose2(colSums(tab))); total <- choose2(length(truth))
  RI <- (total - rows - cols + 2 * both) / total
  expected <- rows * cols / total; denominator <- (rows + cols) / 2 - expected
  ARI <- if (abs(denominator) < 1e-12) {
    as.numeric(abs(both - (rows + cols) / 2) < 1e-12)
  } else (both - expected) / denominator
  c(RI = RI, ARI = ARI)
}

# Match the TWO fitted theta values to generating theta by minimum squared error.
# Apply the SAME permutation to gamma, alpha and labels. Never sort gamma separately.
# Ground truth is used only after fitting to evaluate parameter recovery.
align_fit <- function(fit, theta_truth) {
  permutations <- list(c(1L, 2L), c(2L, 1L))
  errors <- vapply(permutations, function(p) sum((fit$theta[p] - theta_truth)^2), numeric(1))
  p <- permutations[[which.min(errors)]]
  fit$theta <- fit$theta[p]; fit$gamma <- fit$gamma[p]; fit$alpha <- fit$alpha[, p, drop = FALSE]
  fit$labels <- match(fit$labels, p)
  fit$permutation <- p
  fit
}

fit_metrics <- function(fit, data, model) {
  aligned <- align_fit(fit, data$setting$theta)
  scores <- partition_scores(data$z, aligned$labels)
  row <- data.frame(
    setting = data$setting$id, model = model, n = nrow(data$Y), network_seed = data$seed,
    theta1 = aligned$theta[1L], theta2 = aligned$theta[2L], beta = aligned$beta,
    gamma1 = aligned$gamma[1L], gamma2 = aligned$gamma[2L],
    RSE_theta = sqrt(sum((aligned$theta - data$setting$theta)^2)),
    RSE_gamma = sqrt(sum((aligned$gamma - data$setting$gamma)^2)),
    RSE_beta = if (model == "covariate") abs(aligned$beta - data$setting$beta) else NA_real_,
    RI = unname(scores["RI"]), ARI = unname(scores["ARI"]),
    ELBO = fit$elbo, converged = fit$converged, reason = fit$reason,
    all_inner_converged = fit$all_inner_converged,
    outer_iterations = fit$outer_iterations, inner_iterations = fit$inner_iterations,
    selected_start_runtime_s = fit$runtime_s,
    backtracks = fit$backtracks, regularizations = fit$regularizations,
    edges = sum(data$Y) / 2, stringsAsFactors = FALSE
  )
  list(row = row, aligned_fit = aligned)
}

# Manuscript averages RSE within each replication. This differs from an overall RMSE.
summarize_replicates <- function(rows) {
  if (!nrow(rows)) return(data.frame())
  keys <- unique(rows[, c("setting", "model")])
  metrics <- c("theta1", "theta2", "beta", "gamma1", "gamma2", "RSE_theta", "RSE_gamma", "RSE_beta", "RI", "ARI", "ELBO", "runtime_s", "between_model_ARI")
  result <- list(); index <- 0L
  for (r in seq_len(nrow(keys))) {
    group <- rows[rows$setting == keys$setting[r] & rows$model == keys$model[r], ]
    for (subset in c("all_successful", "converged_only")) {
      sub <- if (subset == "converged_only") group[group$converged %in% TRUE, ] else group
      for (metric in metrics) {
        if (!metric %in% names(sub)) next
        values <- sub[[metric]][is.finite(sub[[metric]])]
        index <- index + 1L
        result[[index]] <- data.frame(setting = keys$setting[r], model = keys$model[r], subset = subset,
           metric = metric, n_successful_fits = nrow(group), n_converged_fits = sum(group$converged %in% TRUE),
           n_values = length(values), mean = if (length(values)) mean(values) else NA_real_,
           sd = if (length(values) > 1L) sd(values) else NA_real_, stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(result)) return(data.frame())
  do.call(rbind, result)
}
