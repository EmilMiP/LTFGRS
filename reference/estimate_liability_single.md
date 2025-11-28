# Estimating the genetic or full liability

`estimate_liability_single` estimates the genetic component of the full
liability and/or the full liability for a number of individuals based on
their family history.

## Usage

``` r
estimate_liability_single(
  .tbl = NULL,
  family_graphs = NULL,
  h2 = 0.5,
  pid = "pid",
  fid = "fid",
  family_graphs_col = "fam_graph",
  role = NULL,
  out = c(1),
  tol = 0.01,
  useMixture = FALSE,
  method = "PA"
)
```

## Arguments

- .tbl:

  A matrix, list or data frame that can be converted into a tibble. Must
  have at least five columns that hold the family identifier, the
  personal identifier, the role and the lower and upper thresholds. Note
  that the role must be one of the following abbreviations - `g`
  (Genetic component of full liability) - `o` (Full liability) - `m`
  (Mother) - `f` (Father) - `c[0-9]*.[0-9]*` (Children) - `mgm`
  (Maternal grandmother) - `mgf` (Maternal grandfather) - `pgm`
  (Paternal grandmother) - `pgf` (Paternal grandfather) - `s[0-9]*`
  (Full siblings) - `mhs[0-9]*` (Half-siblings - maternal side) -
  `phs[0-9]*` (Half-siblings - paternal side) - `mau[0-9]*`
  (Aunts/Uncles - maternal side) - `pau[0-9]*` (Aunts/Uncles - paternal
  side). Defaults to `NULL`.

- family_graphs:

  A tibble with columns pid and family_graph_col. See prepare_graph for
  construction of the graphs. The family graphs Defaults to NULL.

- h2:

  A number representing the heritability on liability scale for a single
  phenotype. Must be non-negative. Note that under the liability
  threshold model, the heritability must also be at most 1. Defaults to
  0.5.

- pid:

  A string holding the name of the column in `.tbl` (or `family` and
  `threshs`) that hold the personal identifier(s). Defaults to "PID".

- fid:

  A string holding the name of the column in `.tbl` or `family` that
  holds the family identifier. Defaults to "fid".

- family_graphs_col:

  Name of column with family graphs in family_graphs. Defaults to
  "fam_graph".

- role:

  A string holding the name of the column in `.tbl` that holds the role.
  Each role must be chosen from the following list of abbreviations -
  `g` (Genetic component of full liability) - `o` (Full liability) - `m`
  (Mother) - `f` (Father) - `c[0-9]*.[0-9]*` (Children) - `mgm`
  (Maternal grandmother) - `mgf` (Maternal grandfather) - `pgm`
  (Paternal grandmother) - `pgf` (Paternal grandfather) - `s[0-9]*`
  (Full siblings) - `mhs[0-9]*` (Half-siblings - maternal side) -
  `phs[0-9]*` (Half-siblings - paternal side) - `mau[0-9]*`
  (Aunts/Uncles - maternal side) - `pau[0-9]*` (Aunts/Uncles - paternal
  side). Defaults to "role".

- out:

  A character or numeric vector indicating whether the genetic component
  of the full liability, the full liability or both should be returned.
  If `out = c(1)` or `out = c("genetic")`, the genetic liability is
  estimated and returned. If `out = c(2)` or `out = c("full")`, the full
  liability is estimated and returned. If `out = c(1,2)` or
  `out = c("genetic", "full")`, both components are estimated and
  returned. Defaults to `c(1)`.

- tol:

  A number that is used as the convergence criterion for the Gibbs
  sampler. Equals the standard error of the mean. That is, a tolerance
  of 0.2 means that the standard error of the mean is below 0.2.
  Defaults to 0.01.

- useMixture:

  Logical indicating whether the mixture model should be used to
  calculate the genetic liability. Requires K_i and K_pop columns as
  well as lower and upper. Defaults to FALSE.

- method:

  Estimation method used to estimate the (genetic) liability. Defaults
  to "PA". Current implementation of PA only supports estimates of
  genetic liability. For full or both genetic and full liability
  estimates use "Gibbs".

## Value

If `family` and `threshs` are two matrices, lists or data frames that
can be converted into tibbles, if `family` has two columns named like
the strings represented in `pid` and `fid`, if `threshs` has a column
named like the string given in `pid` as well as a column named "lower"
and a column named "upper" and if the liability-scale heritability `h2`,
`out`, `tol` and `always_add` are of the required form, then the
function returns a tibble with either four or six columns (depending on
the length of out). The first two columns correspond to the columns
`fid` and `pid` ' present in `family`. If `out` is equal to `c(1)` or
`c("genetic")`, the third and fourth column hold the estimated genetic
liability as well as the corresponding standard error, respectively. If
`out` equals `c(2)` or `c("full")`, the third and fourth column hold the
estimated full liability as well as the corresponding standard error,
respectively. If `out` is equal to `c(1,2)` or `c("genetic","full")`,
the third and fourth column hold the estimated genetic liability as well
as the corresponding standard error, respectively, while the fifth and
sixth column hold the estimated full liability as well as the
corresponding standard error, respectively.

## Details

This function can be used to estimate either the genetic component of
the full liability, the full liability or both. It is possible to input
either

## See also

[`future_apply`](https://future.apply.futureverse.org/reference/future_apply.html),
[`estimate_liability_multi`](https://emilmip.github.io/LTFGRS/reference/estimate_liability_multi.md),
[`estimate_liability`](https://emilmip.github.io/LTFGRS/reference/estimate_liability.md)

## Examples

``` r
sims <- simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL,
add_ind = TRUE, h2 = 0.5, n_sim=10, pop_prev = .05)
#
estimate_liability_single(.tbl = sims$thresholds,
h2 = 0.5, pid = "indiv_ID", fid = "fid", role = "role", out = c(1),
tol = 0.01)
#> The number of workers is 1
#> # A tibble: 10 × 4
#>    fid    indiv_ID     est   var
#>    <chr>  <chr>      <dbl> <dbl>
#>  1 fid_1  fid_1    -0.0309 0.481
#>  2 fid_2  fid_2    -0.0169 0.490
#>  3 fid_3  fid_3    -0.0320 0.482
#>  4 fid_4  fid_4    -0.0138 0.491
#>  5 fid_5  fid_5     1.28   0.246
#>  6 fid_6  fid_6     0.476  0.423
#>  7 fid_7  fid_7     0.414  0.420
#>  8 fid_8  fid_8    -0.0131 0.491
#>  9 fid_9  fid_9     0.830  0.417
#> 10 fid_10 fid_10   -0.0369 0.478
#
sims <- simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE,
h2 = 0.5, n_sim=10, pop_prev = .05)
#> Warning: Neither fam_vec nor n_fam is specified...
#
estimate_liability_single(.tbl = sims$thresholds,
h2 = 0.5, pid = "indiv_ID", fid = "fid", role = "role",
out = c("genetic"), tol = 0.01)
#> The number of workers is 1
#> # A tibble: 10 × 4
#>    fid    indiv_ID       est   var
#>    <chr>  <chr>        <dbl> <dbl>
#>  1 fid_1  fid_1    -0.00459  0.494
#>  2 fid_2  fid_2    -0.000385 0.499
#>  3 fid_3  fid_3    -0.00413  0.494
#>  4 fid_4  fid_4    -0.00214  0.497
#>  5 fid_5  fid_5    -0.000385 0.499
#>  6 fid_6  fid_6    -0.000242 0.500
#>  7 fid_7  fid_7    -0.00122  0.498
#>  8 fid_8  fid_8    -0.000612 0.499
#>  9 fid_9  fid_9    -0.00109  0.498
#> 10 fid_10 fid_10   -0.00511  0.493
```
