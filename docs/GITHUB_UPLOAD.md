# Upload the simulation project to GitHub

1. Download and extract the updated ZIP on your computer.
2. Open https://github.com/new and create a repository, for example `covariate-network-clustering-simulations`. Suggested description: `R code for the four simulation settings in Model-Based Clustering of Network Data with Binary Covariate Information.` Choose Public for access by readers/reviewers, or Private while preparing it. Leave the automatic README and gitignore options unselected because both files are supplied.
3. Open the new repository. Use the upload-files link on an empty repository, or `Add file` > `Upload files` on an existing repository.
4. Open the extracted `network-clustering-simulations` folder on your computer. Upload its CONTENTS to the repository root: README.md, R/, config/, scripts/, tests/, docs/, results/README.md, the Rproj file, .gitignore and .github/. Include the hidden dotfiles. Uploading the ZIP alone does not expose the runnable project structure.
5. Use commit message `Add manuscript simulation code and documentation`, then confirm the upload/commit.
6. Check that README.md is at the top level, with R/, config/, scripts/ and tests/ beside it. In particular, scripts/run_simulations.R and .github/workflows/r-checks.yml should exist at those paths.
7. Open Actions and inspect the R verification workflow. It runs the automated checks and a smoke experiment, not the complete 100-replication study.
8. To publish the simulation output tables too, upload the completed full-run summary.csv, replicates.csv, starts.csv, fit_status.csv, configuration.R and sessionInfo.txt together under a clearly identified results subfolder or attach them to a versioned release. Keep smoke outputs separate. The supplied .gitignore excludes generated results during Git-based uploads; add intended publication results explicitly if using Git.

This guide describes the upload steps; the repository has not been created or published by this package.

Official GitHub instructions:
- https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository
- https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository
