# Third-party components and their licences

MR Studio's own code (the R engine in `prototype/engine`, the interface, the analysis contracts and the desktop shell) is released
under the MIT licence (`LICENSE`). The installer also bundles R, R packages and PLINK. Each of these is a separate program under its
own licence, unchanged except where noted. They are listed below with the versions or commits that the build pipeline pins
(`.github/workflows/prototype.yml`, `prototype/ci/runtime.yml`). The licence text of each package is installed with it, in its
folder inside the bundled R library. The exact versions in each installer, including the packages that these depend on, are
recorded by the pipeline in `runtime-lock-<system>.txt`.

Licences are given as each package's `DESCRIPTION` file states them. "GPL-2 | GPL-3" means the user may choose either version.

## Runtime and tools

| Component | Version | Licence | Source |
|---|---|---|---|
| R | 4.5 (conda-forge) | GPL-2 \| GPL-3 | r-project.org |
| PLINK | 1.9, official build of 2026-09-27 (checked by SHA-256) | GPL-3 | cog-genomics.org/plink |
| Tauri (desktop shell) and its Rust crates | see `Cargo.lock` | MIT or Apache-2.0 (Tauri); each crate states its own | tauri.app |

## R packages that MR Studio calls for its methods

| Package | Version or commit | Licence | Source |
|---|---|---|---|
| TwoSampleMR | commit 5e93ee7 | MIT | GitHub MRCIEU/TwoSampleMR |
| MendelianRandomization | 0.10.0 | GPL-2 \| GPL-3 | CRAN |
| MRMix | commit 56afdb2 | GPL-2 | GitHub gqi/MRMix |
| MRPRESSO | commit 3e3c92d | GPL-3 | GitHub rondolab/MR-PRESSO |
| RadialMR | commit 7ca35c0 | GPL-3 | GitHub WSpiller/RadialMR |
| mr.raps | commit dd79b5b (its rsnps import removed before installing) | GPL-3 | GitHub qingyuanzhao/mr.raps |
| simex | 1.8 (CRAN source, checked by SHA-256) | GPL-3 | CRAN |
| GenomicSEM | commit da95d43 | GPL-3 | GitHub GenomicSEM/GenomicSEM |
| simsalapar | 1.0-13 (CRAN source, checked by SHA-256) | GPL-2 \| GPL-3 | CRAN |
| MRcML | commit 58a3481 | MIT | GitHub xue-hr/MRcML |
| MRCD | commit 11ca095 | MIT | GitHub xue-hr/MRCD |
| BiDirectCausal | commit b463ff1 | MIT | GitHub xue-hr/BiDirectCausal |
| lhcMR | commit aee7f34 | MIT | GitHub LizaDarrous/lhcMR |
| mixtools | 2.0.0.1 | GPL (>= 2) | CRAN |
| MRlap | commit 660f026 | GPL-2 | GitHub n-mounier/MRlap |
| MRAPSS | commit 85c18ea | MIT (the repository's `LICENSE` file; its `DESCRIPTION` licence field is left unfilled) | GitHub YangLabHKUST/MR-APSS |
| cause | commit 6257e6c (compiled; where a system cannot build it, the app says CAUSE is not installed) | GPL (>= 2) | GitHub jean997/cause |
| loo | 2.8.0 (pinned for cause, see `prototype/ci/runtime.yml`) | GPL (>= 3) | CRAN |
| boot | 1.3-32 | Unlimited | CRAN |
| ieugwasr | conda-forge | MIT | CRAN |

## Other R packages in the bundled runtime

| Licence | Packages |
|---|---|
| MIT | cli, dplyr, mixsqp, purrr, ggplot2, gtable, jsonlite, later, magrittr, mgsub (1.7.3), nanoparquet, plotly, readr, scales, stringr, tidyr |
| GPL-2 | cowplot, doParallel, glmnet, kernlab, numDeriv |
| GPL-3 | gdata, ggnewscale, rmarkdown, rslurm, shiny, splitstackshape, ggrepel |
| GPL (>= 2) | digest, expm, gridExtra, httpuv, lattice, lavaan, nortest, pbapply, psych, Rcpp, rootSolve, sfsmisc |
| GPL (>= 3) | ashr, pracma, RcppParallel |
| GPL (any version) | gridBase, knitr, segmented |
| GPL-2 \| GPL-3 | e1071, MASS, png |
| LGPL (>= 2) or (>= 2.1) | DBI, R.utils, RSQLite, survival |
| MPL-2.0 | data.table |
| Apache-2.0 | foreach, iterators |
| BSD-3-Clause | colorspace |
| Artistic-2.0 | intervals, matrixStats |

## Not bundled

- **BWMR** (Zhao et al. 2020): its repository states no licence, so it may not be redistributed. MR Studio runs its own implementation
  of the published algorithm instead (MIT, part of MR Studio), validated against the package in `validation/zhao2020_bwmr`.
- **Reference data** (1000 Genomes European LD panel, the LDSC authors' LD scores, OpenGWAS summary statistics) are not in the
  installer. MR Studio downloads them on request from their original sources, under those sources' terms, and checks each file.
