# Verification and completed simulation results

## Software checks

The default degree-initialization path was validated in R 4.3.3 on Linux using base R. The checks in `tests/run_checks.R` cover:

- stable network generation;
- analytic gradient/full Hessian against central finite differences;
- the complete ELBO against exhaustive dyad summation;
- surrogate tangency and minorization;
- the exact two-group quadratic solution against numerical optimization;
- MM/ELBO ascent, positivity and simplex constraints;
- label alignment and RI/ARI fixtures.

End-to-end pipeline checks covered all four settings, subset execution with stable setting-specific seeds, and reuse of completed checkpoints. These small checks are distinct from the complete manuscript experiment.

## Completed manuscript experiment

The supplied completed simulation results contain 100 replications for each of S1–S4 and two selected model fits per replication: 400 generated networks and 800 selected fits. Each setting/model combination contains 100 records; all selected fits are marked converged. The following rounded summaries agree with the updated simulation tables:

| Setting | Mean RSE_theta: network-only | Mean RSE_theta: covariate | Mean ARI: network-only | Mean ARI: covariate |
|---|---:|---:|---:|---:|
| S1 | 0.3317 | 0.0415 | 0.9509 | 0.9513 |
| S2 | 0.3611 | 0.0598 | 0.9445 | 0.9445 |
| S3 | 0.1749 | 0.0455 | 0.9651 | 0.9643 |
| S4 | 0.3331 | 0.0476 | 0.7906 | 0.7906 |

These results support improved recovery of the generating network parameters and broadly comparable clustering accuracy between models. They do not support a claim of consistently improved clustering accuracy.

Run settings and session information should be retained with generated results. The optional `lda`/MMSB path is not the default manuscript configuration and was not part of the local degree-initialization checks.
