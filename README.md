# Multivariate Erlang mixture inference

This repository contains the manuscript and the complete numerical materials
for **“Approximation by mixtures of multivariate Erlang distributions”** by
Hien Duy Nguyen.

The numerical demonstration is intentionally small.  It illustrates the
inferential devices derived in the manuscript for the Erlang-smoothed Bayesian
bootstrap: the local density bias–variance decomposition, the different
resolution choices for density estimation and bounded actuarial functionals,
the finite occupied-cell representation, and exact posterior means and joint
covariances for an aggregate exceedance probability and a capped layer.

## Contents

- `main.tex` — self-contained manuscript source.
- `R/numerical_demonstration.R` — seeded base-R simulation and plotting code.
- `results/density_resolution.csv` — density bias, variance, posterior risk,
  and occupied-component summaries for the resolution sweep.
- `results/rate_summary.csv` — density and actuarial-functional risk summaries
  across sample sizes.
- `results/run_metadata.txt` — seed, replication counts, integration grid,
  empirical slopes, and R session information.
- `figures/numerical_demonstration.pdf` — generated figure used in the manuscript.
- `figures/numerical_demonstration.png` — generated screen-resolution copy.

No external data are used.

## Reproduce the numerical results

From the repository root, run:

```sh
Rscript R/numerical_demonstration.R
```

The script uses only base R and writes the CSV summaries, run metadata, and
both figure formats to their repository paths.  The fixed seed is `20260926`.
On the system recorded in `results/run_metadata.txt`, the full run took less
than one minute.

A short installation check is available as:

```sh
Rscript R/numerical_demonstration.R --quick
```

The quick mode overwrites the result files with low-replication checks; rerun
the full command before using or committing the reported results.

## Numerical design

The target is a dependent bivariate lognormal distribution,

```text
mu    = (0, 0.25)
Sigma = ((0.36, 0.264), (0.264, 0.64))
```

so the log-scale correlation is `0.55`.  The experiment uses sample sizes
`N = 200, 500, 1000, 2000`.

- Local density loss is evaluated on `[0,6]^2` with a `100 x 100` midpoint grid
  and 150 replications per sample size.  The manuscript's rule
  `n_D = ceiling(sqrt(N))` gives `n_D = 15, 23, 32, 45`, respectively.
  A separate 120-replication sweep at `N = 1000` uses
  `n = 20, 32, 50, 100, 200` to display the bias–variance trade-off; `32` is
  the prescribed resolution at that sample size.
- The functional experiment uses 1000 replications and the sufficient
  undersmoothing choice `n_Theta = N`.  It studies
  `Pr(X1 + X2 > 4)` and
  `E[min((X1 + X2 - 2.5)_+, 2.5)]`.
- Reference functional values are obtained by deterministic
  conditional-lognormal quadrature inside the script.

These calculations illustrate the estimator analysed in the manuscript.  They
do not benchmark an EM or adaptive mixture fit, study an uncapped stop loss,
or establish frequentist coverage of Bayesian-bootstrap credible sets.

## Build the manuscript

A standard TeX installation with `latexmk` and pdfLaTeX is sufficient:

```sh
latexmk -pdf main.tex
```

The generated PDF and PNG figures are included in the repository.  Rerun the
numerical script before compiling if you change the plotting code or results.

The bibliography is included directly in `main.tex`; BibTeX and Biber are not
required.

## Rights

The manuscript and code are provided here for scholarly review and
reproducibility.  No separate reuse licence is asserted in this repository.
