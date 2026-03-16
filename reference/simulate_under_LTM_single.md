# Simulate under the liability threshold model (single phenotype).

`simulate_under_LTM_single` simulates families and thresholds under the
liability threshold model for a given family structure and a single
phenotype. Please note that it is not possible to simulate different
family structures.

## Usage

``` r
simulate_under_LTM_single(
  fam_vec = c("m", "f", "s1", "mgm", "mgf", "pgm", "pgf"),
  n_fam = NULL,
  add_ind = TRUE,
  h2 = 0.5,
  n_sim = 1000,
  pop_prev = 0.1
)
```

## Arguments

- fam_vec:

  A vector of strings holding the different family members. All family
  members must be represented by strings from the following list: - `m`
  (Mother) - `f` (Father) - `c[0-9]*.[0-9]*` (Children) - `mgm`
  (Maternal grandmother) - `mgf` (Maternal grandfather) - `pgm`
  (Paternal grandmother) - `pgf` (Paternal grandfather) - `s[0-9]*`
  (Full siblings) - `mhs[0-9]*` (Half-siblings - maternal side) -
  `phs[0-9]*` (Half-siblings - paternal side) - `mau[0-9]*`
  (Aunts/Uncles - maternal side) - `pau[0-9]*` (Aunts/Uncles - paternal
  side). Defaults to `c("m","f","s1","mgm","mgf","pgm","pgf")`.

- n_fam:

  A named vector holding the desired number of family members. See
  [`setNames`](https://rdrr.io/r/stats/setNames.html). All names must be
  picked from the list mentioned above. Defaults to `NULL`.

- add_ind:

  A logical scalar indicating whether the genetic component of the full
  liability as well as the full liability for the underlying target
  individual should be included in the covariance matrix. Defaults to
  `TRUE`.

- h2:

  A number representing the liability-scale heritability for a single
  phenotype. Must be non-negative. Note that under the liability
  threshold model, the heritability must also be at most 1. Defaults to
  0.5.

- n_sim:

  A positive number representing the number of simulations. Defaults to
  1000.

- pop_prev:

  A positive number representing the population prevalence, i.e. the
  overall prevalence in the population. Must be smaller than 1. Defaults
  to 0.1.

## Value

If either `fam_vec` or `n_fam` is used as the argument, if it is of the
required format, if the liability-scale heritability `h2` is a number
satisfying \\0 \leq h^2\\, `n_sim` is a strictly positive number, and
`pop_prev` is a positive number that is at most one, then the output
will be a list holding two tibbles. The first tibble, `sim_obs`, holds
the simulated liabilities, the disease status and the current
age/age-of-onset for all family members in each of the `n_sim` families.
The second tibble, `thresholds`, holds the family identifier, the
personal identifier, the role (specified in fam_vec or n_fam) as well as
the lower and upper thresholds for all individuals in all families. Note
that this tibble has the format required in
[`estimate_liability`](https://emilmip.github.io/LTFGRS/reference/estimate_liability.md).
In addition, note that if neither `fam_vec` nor `n_fam` are specified,
the function returns the disease status, the current age/age-of-onset,
the lower and upper thresholds, as well as the personal identifier for a
single individual, namely the individual under consideration (called
`o`). If both `fam_vec` and `n_fam` are defined, the user is asked to '
decide on which of the two vectors to use.

## See also

[`construct_covmat`](https://emilmip.github.io/LTFGRS/reference/construct_covmat.md),
[`simulate_under_LTM_multi`](https://emilmip.github.io/LTFGRS/reference/simulate_under_LTM_multi.md),
[`simulate_under_LTM`](https://emilmip.github.io/LTFGRS/reference/simulate_under_LTM.md)

## Examples

``` r
simulate_under_LTM_single()
#> $sim_obs
#> # A tibble: 1,000 × 26
#>    fid         g      o       m       f     s1    mgm     mgf    pgm    pgf
#>    <chr>   <dbl>  <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>  <dbl>  <dbl>
#>  1 fid_1   0.321  0.536 -0.0243  1.59    0.827 -1.07   0.0377  0.375  0.225
#>  2 fid_2   0.694  0.797  0.517  -0.639   0.245  0.214  1.28   -0.236 -0.144
#>  3 fid_3  -1.18  -1.83  -0.0783  0.171   0.349 -1.46   0.828   1.05  -1.62 
#>  4 fid_4   0.229 -1.16   0.0353  0.237   1.79   2.19  -0.721   0.209  0.612
#>  5 fid_5   0.154 -1.14  -0.930  -0.384  -0.698 -1.36   0.392   0.751 -0.906
#>  6 fid_6  -0.665 -0.452 -1.63    0.416   1.52  -2.15  -0.986   0.134  1.72 
#>  7 fid_7   0.563  0.438 -0.979   1.59   -0.521  1.81   0.807   0.172  1.55 
#>  8 fid_8  -0.887 -1.42   1.14   -2.39   -0.307  0.960  0.390  -0.506 -0.407
#>  9 fid_9  -0.263 -0.929 -0.236  -0.0244  0.891 -0.336  0.320   0.341 -1.79 
#> 10 fid_10 -0.673  0.673  0.198  -0.448  -1.64  -0.710  2.13   -0.230 -1.52 
#> # ℹ 990 more rows
#> # ℹ 16 more variables: o_status <lgl>, m_status <lgl>, f_status <lgl>,
#> #   s1_status <lgl>, mgm_status <lgl>, mgf_status <lgl>, pgm_status <lgl>,
#> #   pgf_status <lgl>, o_aoo <dbl>, m_aoo <dbl>, f_aoo <dbl>, s1_aoo <dbl>,
#> #   mgm_aoo <dbl>, mgf_aoo <dbl>, pgm_aoo <dbl>, pgf_aoo <dbl>
#> 
#> $thresholds
#> # A tibble: 8,000 × 5
#>    fid    indiv_ID role  lower upper
#>    <chr>  <chr>    <chr> <dbl> <dbl>
#>  1 fid_1  fid_1_1  o      -Inf  2.95
#>  2 fid_2  fid_2_1  o      -Inf  3.35
#>  3 fid_3  fid_3_1  o      -Inf  2.51
#>  4 fid_4  fid_4_1  o      -Inf  3.31
#>  5 fid_5  fid_5_1  o      -Inf  3.14
#>  6 fid_6  fid_6_1  o      -Inf  3.06
#>  7 fid_7  fid_7_1  o      -Inf  2.47
#>  8 fid_8  fid_8_1  o      -Inf  3.38
#>  9 fid_9  fid_9_1  o      -Inf  2.59
#> 10 fid_10 fid_10_1 o      -Inf  3.17
#> # ℹ 7,990 more rows
#> 
simulate_under_LTM_single(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2),
c("m","mgm","mgf","mhs")))
#> $sim_obs
#> # A tibble: 1,000 × 20
#>    fid         g       o       m    mgm    mgf    mhs1    mhs2 o_status m_status
#>    <chr>   <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>   <dbl> <lgl>    <lgl>   
#>  1 fid_1 -1.38   -1.62   -1.35    1.79  -1.37  -1.71   -0.595  FALSE    FALSE   
#>  2 fid_2  1.09    0.548   0.814  -0.357  0.750  0.970   0.618  FALSE    FALSE   
#>  3 fid_3 -0.797  -1.22   -2.53   -2.15  -1.09  -0.584  -1.68   FALSE    FALSE   
#>  4 fid_4  0.0998  1.08    0.527  -0.148 -0.758 -0.175   0.511  FALSE    FALSE   
#>  5 fid_5 -1.05   -1.93   -0.0152 -0.693  0.841 -0.236   0.242  FALSE    FALSE   
#>  6 fid_6  0.139   0.678  -1.00   -0.528 -1.73  -0.332  -0.118  FALSE    FALSE   
#>  7 fid_7 -0.216  -1.56    0.671   0.103 -0.419 -0.0806 -0.0516 FALSE    FALSE   
#>  8 fid_8  0.387  -0.0113  1.42    0.355  0.218 -0.0232 -0.885  FALSE    TRUE    
#>  9 fid_9  0.557   1.84    0.823  -0.268  1.46   0.0325 -0.699  TRUE     FALSE   
#> 10 fid_…  0.128   0.515   0.756   0.673 -0.163 -1.71    0.0366 FALSE    FALSE   
#> # ℹ 990 more rows
#> # ℹ 10 more variables: mgm_status <lgl>, mgf_status <lgl>, mhs1_status <lgl>,
#> #   mhs2_status <lgl>, o_aoo <dbl>, m_aoo <dbl>, mgm_aoo <dbl>, mgf_aoo <dbl>,
#> #   mhs1_aoo <dbl>, mhs2_aoo <dbl>
#> 
#> $thresholds
#> # A tibble: 6,000 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o     -Inf     3.10
#>  2 fid_2  fid_2_1  o     -Inf     3.31
#>  3 fid_3  fid_3_1  o     -Inf     2.47
#>  4 fid_4  fid_4_1  o     -Inf     2.87
#>  5 fid_5  fid_5_1  o     -Inf     3.03
#>  6 fid_6  fid_6_1  o     -Inf     2.91
#>  7 fid_7  fid_7_1  o     -Inf     3.31
#>  8 fid_8  fid_8_1  o     -Inf     3.03
#>  9 fid_9  fid_9_1  o        1.85  1.85
#> 10 fid_10 fid_10_1 o     -Inf     2.87
#> # ℹ 5,990 more rows
#> 
simulate_under_LTM_single(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE,
h2 = 0.5, n_sim = 500, pop_prev = .05)
#> $sim_obs
#> # A tibble: 500 × 10
#>    fid          m       f      s1 m_status f_status s1_status m_aoo f_aoo s1_aoo
#>    <chr>    <dbl>   <dbl>   <dbl> <lgl>    <lgl>    <lgl>     <dbl> <dbl>  <dbl>
#>  1 fid_1   0.202  -0.753  -2.67   FALSE    FALSE    FALSE        44    40     17
#>  2 fid_2  -0.286   0.553  -0.0829 FALSE    FALSE    FALSE        65    62     38
#>  3 fid_3   0.0843  0.534  -0.326  FALSE    FALSE    FALSE        59    55     36
#>  4 fid_4  -1.43   -3.12   -2.17   FALSE    FALSE    FALSE        32    35     14
#>  5 fid_5   1.84    0.627   0.797  TRUE     FALSE    FALSE        65    48     19
#>  6 fid_6   0.261   2.17    0.359  FALSE    TRUE     FALSE        53    53     27
#>  7 fid_7  -0.297  -1.33    0.158  FALSE    FALSE    FALSE        56    45     26
#>  8 fid_8  -2.33    0.324   0.768  FALSE    FALSE    FALSE        46    40     20
#>  9 fid_9  -0.668   0.155  -1.04   FALSE    FALSE    FALSE        48    50     27
#> 10 fid_10 -0.607  -0.0234  0.953  FALSE    FALSE    FALSE        43    41     19
#> # ℹ 490 more rows
#> 
#> $thresholds
#> # A tibble: 1,500 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  m     -Inf     2.51
#>  2 fid_2  fid_2_1  m     -Inf     1.84
#>  3 fid_3  fid_3_1  m     -Inf     1.99
#>  4 fid_4  fid_4_1  m     -Inf     2.97
#>  5 fid_5  fid_5_1  m        1.84  1.84
#>  6 fid_6  fid_6_1  m     -Inf     2.18
#>  7 fid_7  fid_7_1  m     -Inf     2.08
#>  8 fid_8  fid_8_1  m     -Inf     2.44
#>  9 fid_9  fid_9_1  m     -Inf     2.36
#> 10 fid_10 fid_10_1 m     -Inf     2.55
#> # ℹ 1,490 more rows
#> 
simulate_under_LTM_single(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5,
n_sim = 200, pop_prev = 0.05)
#> Warning: Neither fam_vec nor n_fam is specified...
#> $sim_obs
#> # A tibble: 200 × 5
#>    fid         g      o o_status o_aoo
#>    <chr>   <dbl>  <dbl> <lgl>    <dbl>
#>  1 fid_1   0.377 -0.223 FALSE       10
#>  2 fid_2  -0.189 -0.622 FALSE       27
#>  3 fid_3  -0.901 -0.657 FALSE       36
#>  4 fid_4  -0.298 -1.36  FALSE       17
#>  5 fid_5  -0.718 -1.32  FALSE       33
#>  6 fid_6   0.886  0.775 FALSE       28
#>  7 fid_7   0.242  1.93  TRUE        61
#>  8 fid_8   0.547  1.32  FALSE       23
#>  9 fid_9   0.290  1.04  FALSE       13
#> 10 fid_10 -0.691 -0.766 FALSE       34
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 200 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o     -Inf     3.73
#>  2 fid_2  fid_2_1  o     -Inf     3.16
#>  3 fid_3  fid_3_1  o     -Inf     2.82
#>  4 fid_4  fid_4_1  o     -Inf     3.50
#>  5 fid_5  fid_5_1  o     -Inf     2.94
#>  6 fid_6  fid_6_1  o     -Inf     3.12
#>  7 fid_7  fid_7_1  o        1.93  1.93
#>  8 fid_8  fid_8_1  o     -Inf     3.30
#>  9 fid_9  fid_9_1  o     -Inf     3.63
#> 10 fid_10 fid_10_1 o     -Inf     2.90
#> # ℹ 190 more rows
#> 
```
