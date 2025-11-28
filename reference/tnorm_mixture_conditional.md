# Title: Calculates mean and variance of mixture of two truncated normal distributions

Title: Calculates mean and variance of mixture of two truncated normal
distributions

## Usage

``` r
tnorm_mixture_conditional(mu, var, lower, upper, K_i, K_pop)
```

## Arguments

- mu:

  Mean value of normal distribution.

- var:

  Variance of normal distribution.

- lower:

  Lower threshold (can be -Inf).

- upper:

  Upper threshold (can be Inf).

- K_i:

  (Stratified) cumulative incidence proportion for the individual.

- K_pop:

  Population prevalence (cumulative incidence proportion).

## Value

mean and variance of mixture distribution between two truncated normal
distributions

## Examples

``` r
tnorm_mixture_conditional(mu = 0, var = 1, lower = -Inf, upper = Inf, K_i = 0, K_pop = 0.01)
#> $mean
#> [1] 0
#> 
#> $var
#> [1] 1
#> 
tnorm_mixture_conditional(mu = 0, var = 1, lower = -Inf, upper = 2, K_i = .01, K_pop = 0.05)
#> $mean
#> [1] 0.04287187
#> 
#> $var
#> [1] 1.083906
#> 
```
