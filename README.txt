Model-Based Clustering of Network Data with Binary Covariate Information
Authors: Paul H. Atandoh and Kevin H. Lee.
This repository contains the R code for the four simulation settings reported in the manuscript. Each setting uses 100 independently generated networks with 100 nodes and two latent groups. A network-only model and a covariate-informed model are fitted to the same network within each replication.
The simulations examine group imbalance, weaker covariate effects and reduced structural separation. Incorporating covariate information substantially improves recovery of the generating network parameters, while clustering accuracy is broadly comparable between the two fitted models. Reduced structural separation presents the greatest clustering challenge among the settings examined.
Data-generating model
For an undirected binary network, edges are independent conditional on latent group memberships and the observed binary covariate:
```text
logit P(Y_ij = 1 | z, x) = theta[z_i] + theta[z_j] + beta * I(x_i == x_j)
```
Latent group memberships are sampled independently with probabilities `gamma`. The binary covariate is generated independently of latent membership and has 50 nodes in each category for n=100. The network-only fit fixes beta=0; the covariate-informed fit estimates beta. Both fits analyze the same covariate-generated network.
Setting	Description	n	theta1	theta2	beta	gamma1	gamma2
S1	Balanced groups; base setting	100	-0.50	0.50	1.0	0.50	0.50
S2	Unbalanced latent groups	100	-0.50	0.50	1.0	0.25	0.75
S3	Weaker covariate effect	100	-0.50	0.50	0.5	0.50	0.50
S4	Reduced structural separation	100	-0.35	0.35	1.0	0.50	0.50
Software and files
The manuscript simulation configuration uses base R and degree-based initialization. No additional R packages are required for this configuration; R 4.1 or later is recommended.
Path	Purpose
`config/settings.R`	Four generating settings and simulation controls
`R/model.R`	Network generation, initialization and model fitting
`R/metrics.R`	Label alignment, parameter error and clustering metrics
`R/experiment.R`	Replication driver, summaries and checkpoints
`scripts/run_simulations.R`	Main simulation entry point
`tests/run_checks.R`	Mathematical and implementation checks
`docs/implementation.md`	Algorithm, metrics and output conventions
`docs/verification.md`	Software verification and completed-run evidence
Open `network-clustering-simulations.Rproj` in RStudio or Posit Cloud and run commands from the project root, where `README.md`, `R/` and `config/` are located.
Run the simulation study
*Make sure to have all files in the same directory before running the codes
First run the checks in the R Console:
```r
source("tests/run_checks.R")
```
Then run the complete experiment:
```r
source("scripts/run_simulations.R")
```
The defaults in `config/settings.R` are:
100 replications per setting, with n=100 and K=2;
five shared initialization states per replication for the paired fits;
base seed 20261012;
degree-based initialization;
maximum 500 outer and 500 inner iterations;
inner and outer tolerances of 1e-4;
membership probability floor of 1e-8.
Results are written to `results/full/`. The complete study produces 400 generated networks and 800 selected model fits.
The equivalent Terminal command is:
```sh
Rscript scripts/run_simulations.R
```
A small pipeline check can be run separately:
```sh
Rscript scripts/run_simulations.R --smoke
```
The smoke configuration uses n=30, two replications per setting and two starts. Its results are for checking execution and are separate from the manuscript's n=100 experiment.
Resume after an interruption
Reopen the same project and rerun the same command with the same output folder and configuration. Keep all files under `results/full/`.
Completed replication checkpoints are reused. An unfinished replication is rerun. Configuration or R-version changes require a new output folder. For example, a separate output folder can be selected from the Terminal:
```sh
Rscript scripts/run_simulations.R --output results/new-run
```
Results and reporting
Each output folder contains:
File	Contents
`replicates.csv`	Selected-fit parameters, errors, RI/ARI, convergence, seeds and runtime
`starts.csv`	Diagnostics for every initialization
`summary.csv`	Parameter/error/clustering means and sample standard deviations by setting and model
`fit_status.csv`	Counts of converged and failed starts
`S1/` through `S4/`	Replication checkpoints and selected-fit RDS objects
`manifest.rds`, `configuration.R`, `sessionInfo.txt`	Run settings and reproducibility metadata
For each replication, parameter errors are defined as:
```text
RSE_theta = sqrt(sum((theta_hat - theta)^2))
RSE_gamma = sqrt(sum((gamma_hat - gamma)^2))
RSE_beta  = abs(beta_hat - beta)
```
Errors for theta and gamma are calculated jointly across their components. The reported means and standard deviations summarize these per-replication errors. RI and ARI compare the estimated partition with the generating latent-group partition.
Labels are aligned consistently for evaluation: the same permutation is applied to theta, gamma, membership probabilities and hard labels. True generating labels/parameters are not used to initialize or fit either model.
The highest-ELBO converged start is selected. If no start converges, the highest-ELBO finite start is retained and explicitly flagged. `summary.csv` contains both `all_successful` and `converged_only` subsets. In the completed manuscript run, all 800 selected fits converged, so these subsets describe the same selected fits.
`runtime_s` includes all fitting starts for a model; `selected_start_runtime_s` is the selected fit's runtime. Initialization and data-generation time are excluded from model-fitting runtime.
Generated files are excluded from version control by `.gitignore`. Retain the complete output directory locally for restart and reproducibility. The `--save-networks` option additionally saves each generated network, covariate vector and generating labels.
Verification
The automated checks cover derivatives, the ELBO, surrogate minorization/tangency, E-step ascent, positivity, simplex constraints, label alignment and RI/ARI. GitHub Actions runs these checks and a small end-to-end simulation after a push or pull request. The completed 100-replication results are documented separately in `docs/verification.md`.
See `docs/GITHUB_UPLOAD.md` for upload instructions. The real-data applications are organized in a separate analysis project.