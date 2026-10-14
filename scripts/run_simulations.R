# Run from the repository root: Rscript scripts/run_simulations.R --smoke
source("config/settings.R")
source("R/model.R")
source("R/metrics.R")
source("R/experiment.R")
args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  idx <- match(paste0("--", name), args)
  if (is.na(idx)) default else {
    if (idx == length(args)) stop("Missing value after --", name)
    args[idx + 1L]
  }
}
settings <- simulation_settings(); control <- simulation_control()
smoke <- "--smoke" %in% args
if (smoke) {
  for (k in seq_along(settings)) settings[[k]]$n <- 30L
  control$reps <- 2L; control$starts <- 2L
  control$max_outer <- 150L; control$max_inner <- 100L
}
control$reps <- as.integer(option("reps", control$reps))
control$starts <- as.integer(option("starts", control$starts))
control$base_seed <- as.integer(option("seed", control$base_seed))
control$initialization <- option("init", control$initialization)
control$save_networks <- "--save-networks" %in% args
selected <- strsplit(option("settings", "S1,S2,S3,S4"), ",", fixed = TRUE)[[1L]]
if (any(!selected %in% names(settings))) stop("Unknown setting; choose S1,S2,S3,S4.")
# Keep master setting positions for stable seeds, even when running one setting alone.
full_settings <- settings
settings <- full_settings[selected]
# run_experiment uses the master index recorded in each setting.
for (id in names(settings)) settings[[id]]$seed_index <- match(id, names(full_settings))
output <- option("output", if (smoke) "results/smoke" else "results/full")
run_experiment(settings, control, output)
message("Finished. Inspect ", output, "/summary.csv and starts.csv before using results.")
