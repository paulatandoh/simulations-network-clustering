# Each replicate generates ONE network; both models fit that same network.
# Independent data-generation and initialization seeds are saved explicitly.
run_replication <- function(setting, replication, setting_number, control, folder) {
  seed <- as.integer(control$base_seed + setting_number * 100000L + replication * 1000L)
  data <- simulate_network(setting, seed)
  data_id <- sprintf("%s_rep%03d", setting$id, replication)
  if (control$save_networks) saveRDS(data, file.path(folder, paste0(data_id, "_data.rds")))
  fits <- list(network_only = list(), covariate = list()); starts <- list(); index <- 0L
  totals <- c(network_only = 0, covariate = 0)
  for (start in seq_len(control$starts)) {
    start_seed <- seed + start
    initial <- initialize_fit(data$Y, start_seed, control$initialization, control$alpha_floor)
    for (model in names(fits)) {
      index <- index + 1L; began <- proc.time()[[3L]]
      fit <- tryCatch(fit_vem(data$Y, data$x, model == "covariate", initial,
          control$max_outer, control$max_inner, control$tol_outer, control$tol_inner,
          control$alpha_floor), error = function(e) e)
      elapsed <- unname(proc.time()[[3L]] - began); totals[model] <- totals[model] + elapsed
      failed <- inherits(fit, "error")
      starts[[index]] <- data.frame(setting = setting$id, replication = replication,
        model = model, start = start, network_seed = seed, initialization_seed = start_seed,
        numerical_failure = failed, error = if (failed) conditionMessage(fit) else "",
        converged = if (failed) FALSE else fit$converged,
        reason = if (failed) "error" else fit$reason,
        ELBO = if (failed) NA_real_ else fit$elbo, runtime_s = elapsed,
        outer_iterations = if (failed) NA_integer_ else fit$outer_iterations,
        inner_iterations = if (failed) NA_integer_ else fit$inner_iterations)
      if (!failed) { fit$start <- start; fit$initialization_seed <- start_seed; fits[[model]][[length(fits[[model]]) + 1L]] <- fit }
    }
  }
  selected <- list(); rows <- list()
  for (model in names(fits)) {
    candidates <- fits[[model]]
    if (!length(candidates)) next
    # Prefer converged starts, then select the largest ELBO (fixed model dimension).
    ok <- vapply(candidates, function(f) f$converged, logical(1))
    pool <- if (any(ok)) which(ok) else seq_along(candidates)
    best <- candidates[[pool[which.max(vapply(candidates[pool], function(f) f$elbo, numeric(1)))]]]
    selected[[model]] <- best
    metrics <- fit_metrics(best, data, model)
    row <- metrics$row; row$replication <- replication; row$selected_start <- best$start
    row$initialization_seed <- best$initialization_seed; row$runtime_s <- unname(totals[model])
    row$failed_starts <- sum(vapply(starts, function(x) x$model == model && x$numerical_failure, logical(1)))
    rows[[model]] <- row
    saveRDS(metrics$aligned_fit, file.path(folder, paste0(data_id, "_", model, "_fit.rds")))
  }
  if (length(rows)) {
    agreement <- if (length(selected) == 2L) unname(partition_scores(selected$network_only$labels, selected$covariate$labels)["ARI"]) else NA_real_
    for (model in names(rows)) rows[[model]]$between_model_ARI <- agreement
  }
  list(rows = if (length(rows)) do.call(rbind, rows) else data.frame(),
       starts = do.call(rbind, starts))
}

run_experiment <- function(settings, control, output = "results", resume = TRUE) {
  if (control$reps < 1L || control$starts < 1L) stop("Use positive reps and starts.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  manifest <- list(settings = settings, control = control, R_version = R.version.string,
                   RNGkind = RNGkind(), format_version = 1L)
  manifest_path <- file.path(output, "manifest.rds")
  if (file.exists(manifest_path)) {
    existing <- readRDS(manifest_path)
    if (!identical(existing, manifest)) stop("Settings, control, or R version changed. Choose a new output directory.")
    if (!resume) stop("Existing output directory. Use resume=TRUE or a new directory.")
  } else saveRDS(manifest, manifest_path)
  dput(list(settings = settings, control = control), file = file.path(output, "configuration.R"))
  capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
  all_rows <- list(); all_starts <- list(); task <- 0L
  for (setting_number in seq_along(settings)) {
    setting <- settings[[setting_number]]
    folder <- file.path(output, setting$id)
    dir.create(folder, showWarnings = FALSE)
    for (replication in seq_len(control$reps)) {
      task <- task + 1L
      path <- file.path(folder, sprintf("rep_%03d.rds", replication))
      if (resume && file.exists(path)) {
        result <- readRDS(path)
        message(setting$id, " replicate ", replication, ": resumed")
      } else {
        message(setting$id, " replicate ", replication, "/", control$reps)
        seed_index <- if (is.null(setting$seed_index)) match(setting$id, c("S1", "S2", "S3", "S4")) else setting$seed_index
        result <- run_replication(setting, replication, seed_index, control, folder)
        temp <- paste0(path, ".tmp")
        saveRDS(result, temp)
        if (!file.rename(temp, path)) stop("Could not save replicate checkpoint.")
      }
      if (nrow(result$rows)) all_rows[[length(all_rows) + 1L]] <- result$rows
      all_starts[[task]] <- result$starts
      # Refresh cumulative exports after EVERY replication.
      rows <- if (length(all_rows)) do.call(rbind, all_rows) else data.frame()
      starts <- do.call(rbind, all_starts)
      write.csv(rows, file.path(output, "replicates.csv"), row.names = FALSE)
      write.csv(starts, file.path(output, "starts.csv"), row.names = FALSE)
      write.csv(summarize_replicates(rows), file.path(output, "summary.csv"), row.names = FALSE)
      failure_summary <- aggregate(cbind(numerical_failure, converged) ~ setting + model,
                                   data = starts, FUN = sum)
      write.csv(failure_summary, file.path(output, "fit_status.csv"), row.names = FALSE)
    }
  }
  invisible(list(replicates = rows, starts = starts, summary = summarize_replicates(rows)))
}
