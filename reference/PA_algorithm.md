# Title Pearson-Aitken algorithm to calculate mean values in truncated multivariate normal distributions

Title Pearson-Aitken algorithm to calculate mean values in truncated
multivariate normal distributions

## Usage

``` r
PA_algorithm(mu, covmat, target_id, lower, upper, K_i = NA, K_pop = NA)
```

## Arguments

- mu:

  vector of means

- covmat:

  covariance matrix, contaning kinship coefficient and heritability on
  each entry (except diagnoal, which is 1 for full liabilities and h2
  for genetic liabilities)

- target_id:

  ID of target individual (or genetic liability), i.e. rowname in covmat
  to return expected genetic liability for

- lower:

  vector of lower thresholds

- upper:

  vector of upper thresholds

- K_i:

  vector of stratified CIPs for each individual. Only used for
  estimating genetic liability under the mixture model.

- K_pop:

  vector of population CIPs. Only used for estimating genetic liability
  under the mixture model.

## Value

A list with two elements: est (expected genetic liability, given input
data) and var (variance of genetic liability, given input data).
