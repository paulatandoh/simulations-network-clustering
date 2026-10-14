# Simulation implementation

## Paired experiment

Every replication generates one undirected binary network using the setting's theta, beta and gamma. The observed binary covariate is balanced independently of latent labels. Both fitted models use that same network and the same initialization states. Network-only fitting omits the covariate coefficient; it does not generate a second network without covariates.

## Fitting

The E-step maximizes a separable MM surrogate. With two groups, each node's simplex-constrained quadratic maximization has a closed-form solution. Updated probabilities satisfy the explicit lower bound alpha_ik >= 1e-8. The M-step jointly updates network and covariate coefficients with full Newton derivatives and backtracking; all cross-Hessian terms are included. Stable logistic and softplus calculations avoid avoidable overflow. Any information-matrix regularization and backtracking are logged.

The complete ELBO sums the dyad term and mixing/entropy term across all applicable nodes/groups. Outer convergence requires the implementation's tolerance condition and a converged inner solve. Failed starts and fallback finite nonconverged selections remain visible in output diagnostics.

## Initialization and evaluation

Degree-based initialization is the configuration used for the manuscript simulation run. An optional MMSB path (`--init mmsb`) requires `lda`, and a random-initialization path is also available. These are alternative configurations; use separate output folders when changing initialization.

No generating labels or parameter values enter fitting or initialization. For evaluation, a theta-based permutation aligns the estimated groups and is applied jointly to theta, gamma, alpha and hard labels.

Theta/gamma error is the Euclidean distance across their components for each replication; beta error is absolute error. Means and sample SDs across replications are reported, rather than a pooled RMSE. RI and ARI evaluate latent-partition recovery; between-model ARI measures agreement between the fitted partitions.

## Files and restart behavior

A configuration manifest guards completed-replication reuse. Completed replication records are reused when restarting with the same settings/control and R version; an unfinished replication is rerun. Changes require a separate output directory. Selected-fit RDS objects retain parameter estimates and ELBO traces. Runtime conventions are given in the root README.
