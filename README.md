# POISE

<!-- badges: start -->
<!-- badges: end -->

**POISE** provides fast, non-iterative estimation of the rate ratio
between two groups of negative binomial counts. A one-step estimator
anchored at the Poisson group rates reproduces the negative binomial
maximum likelihood fit and returns a matching closed-form standard error
and Wald inference, without the iterative joint optimization of
`MASS::glm.nb`. The package also reconstructs the maximum likelihood
analysis from published group summaries, sizes trials following the
method of Tang (2015), and compares operating characteristics by
simulation.

## Installation

You can install the development version of POISE from
[GitHub](https://github.com/gosukehommaEX/POISE) with:

``` r
# install.packages("pak")
pak::pak("gosukehommaEX/POISE")
```

## Example

Estimate the rate ratio from subject-level data:

``` r
library(POISE)

set.seed(1)
n     <- 200
fu    <- runif(n, 0.5, 1)
grp   <- rep(c("T", "C"), each = n / 2)
theta <- ifelse(grp == "T", 0.5, 1)
y     <- rnbinom(n, mu = fu * theta, size = 2)

EstRROneStep(y, fu, grp, control = "C")
```

Reconstruct the same analysis from group-level summaries only, without
individual patient data:

``` r
EstRRSummary(rate_T = 0.9, rate_C = 1.5,
             followup_T = 1, followup_C = 1,
             n_T = 200, n_C = 200, nu = 1.5)
```

Size a trial following Tang (2015):

``` r
SampleSizeTang(rate_T = 0.2, rate_C = 0.4, tau = 2, nu = 1,
               dropout_T = 0.2, dropout_C = 0.2, beta = 0.10)$n_C
```

## Main functions

- `EstRROneStep()` one-step closed-form rate ratio estimator.
- `EstRRSummary()` reconstruct the maximum likelihood fit from summaries.
- `print.rrfit()` formatted display of a fitted rate ratio object.
- `GenSimData()` simulate two-group negative binomial data.
- `CompareOC()` compare operating characteristics of the estimators.
- `ReconError()` quantify the summary reconstruction error.
- `SampleSizeTang()` sample size for a negative binomial rate comparison.

## Reproducing the manuscript

The scripts in `inst/scripts` regenerate the manuscript figures and
tables. Run `data_generate_main.R` first to cache the simulation results
under `inst/scripts/data`, then `table_and_figure_manuscript.R` to write
the figures (PDF and EPS) and tables (LaTeX) to `inst/scripts/results`.

## Reference

Tang Y (2015). Sample size estimation for negative binomial regression
comparing rates of recurrent events with unequal follow-up time.
*Journal of Biopharmaceutical Statistics*.
