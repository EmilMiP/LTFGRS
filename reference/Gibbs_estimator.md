# Wrapper around the Gibbs Sampler that returns formatted liability estimates for the proband

Wrapper around the Gibbs Sampler that returns formatted liability
estimates for the proband

## Usage

``` r
Gibbs_estimator(cov, tbl, out, tol = 0.01, burn_in = 1000, phen_names = NULL)
```

## Arguments

- cov:

  Covariance (kinship matrix times heritability with corrected diagonal)
  matrix

- tbl:

  Tibble with lower and upper bounds for the Gibbs sampler

- out:

  Vector indicating if genetic ans/or full liabilities should be
  estimated

- tol:

  Convergence criteria, tolerance

- burn_in:

  Number of burn-in iterations

- phen_names:

  Names of the phenotypes being analyzed

## Value

Formatted liability estimate(s) and standard error(s) of the mean for
the proband.

## Examples

``` r
# uninformative sampling:
Gibbs_estimator(cov = diag(3), tbl = tibble::tibble(lower = rep(-Inf, 3),
upper = rep(Inf, 3)), out = 1:2, tol = 0.01, burn_in = 1000)
#> # A tibble: 1 × 4
#>   genetic_est genetic_se full_est full_se
#>         <dbl>      <dbl>    <dbl>   <dbl>
#> 1   -0.000781    0.00318  0.00394 0.00326
```
