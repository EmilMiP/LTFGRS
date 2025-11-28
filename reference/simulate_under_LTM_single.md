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
#>    fid          g       o       m       f      s1    mgm     mgf     pgm    pgf
#>    <chr>    <dbl>   <dbl>   <dbl>   <dbl>   <dbl>  <dbl>   <dbl>   <dbl>  <dbl>
#>  1 fid_1   0.196   0.623  -0.832  -0.854  -0.300  -1.40   0.0735  1.86   -2.10 
#>  2 fid_2  -0.123   0.128   1.52    0.587   0.350  -1.19   1.63   -0.872   2.04 
#>  3 fid_3   0.544   0.258   1.23   -0.669  -1.00    1.39  -2.01   -0.164   0.804
#>  4 fid_4   0.0890 -0.0819 -0.134   0.0406 -0.931  -0.325  0.0917  1.07   -0.362
#>  5 fid_5  -0.0975  0.0639 -1.48   -1.50    0.891   0.146 -0.104  -1.71   -0.654
#>  6 fid_6   1.22    0.481  -1.23    1.40    0.0379  0.350 -1.18    1.58   -0.499
#>  7 fid_7   1.45    2.19   -0.0677  2.23    1.10    0.298  0.225   2.05    2.31 
#>  8 fid_8  -0.327  -0.482   0.285  -0.807  -1.02   -0.508 -0.620   0.391   0.181
#>  9 fid_9   0.563  -0.426  -0.284  -0.0952  0.679   0.349 -0.644  -0.662   0.745
#> 10 fid_10  0.241   1.75    0.307  -0.0360 -0.467  -0.131  1.11    0.0888 -0.809
#> # ℹ 990 more rows
#> # ℹ 16 more variables: o_status <lgl>, m_status <lgl>, f_status <lgl>,
#> #   s1_status <lgl>, mgm_status <lgl>, mgf_status <lgl>, pgm_status <lgl>,
#> #   pgf_status <lgl>, o_aoo <dbl>, m_aoo <dbl>, f_aoo <dbl>, s1_aoo <dbl>,
#> #   mgm_aoo <dbl>, mgf_aoo <dbl>, pgm_aoo <dbl>, pgf_aoo <dbl>
#> 
#> $thresholds
#> # A tibble: 8,000 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o     -Inf     2.91
#>  2 fid_2  fid_2_1  o     -Inf     2.91
#>  3 fid_3  fid_3_1  o     -Inf     3.28
#>  4 fid_4  fid_4_1  o     -Inf     2.99
#>  5 fid_5  fid_5_1  o     -Inf     3.35
#>  6 fid_6  fid_6_1  o     -Inf     3.24
#>  7 fid_7  fid_7_1  o        2.18  2.18
#>  8 fid_8  fid_8_1  o     -Inf     2.87
#>  9 fid_9  fid_9_1  o     -Inf     3.55
#> 10 fid_10 fid_10_1 o        1.74  1.74
#> # ℹ 7,990 more rows
#> 
simulate_under_LTM_single(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2),
c("m","mgm","mgf","mhs")))
#> $sim_obs
#> # A tibble: 1,000 × 20
#>    fid           g       o       m     mgm    mgf   mhs1      mhs2 o_status
#>    <chr>     <dbl>   <dbl>   <dbl>   <dbl>  <dbl>  <dbl>     <dbl> <lgl>   
#>  1 fid_1  -1.13    -0.712  -0.226   0.421  -1.90   0.727  0.000923 FALSE   
#>  2 fid_2  -0.210   -0.143   0.0987  0.603  -0.264  1.51   0.375    FALSE   
#>  3 fid_3   0.00321 -0.0803 -0.813  -0.339   1.81   0.586 -0.311    FALSE   
#>  4 fid_4   0.811    0.354  -0.495  -0.0963  0.601  0.944  0.730    FALSE   
#>  5 fid_5  -0.700   -1.85   -0.136   0.567  -0.467 -1.06   0.155    FALSE   
#>  6 fid_6  -0.139    0.319   0.191   1.05   -2.22   0.145 -0.829    FALSE   
#>  7 fid_7   0.991    1.31    0.130  -0.566  -0.461 -0.152 -1.80     TRUE    
#>  8 fid_8   0.412    0.794   2.09    0.967   1.02   0.533 -0.0306   FALSE   
#>  9 fid_9  -0.562    0.370  -0.514  -0.245   0.963 -0.206 -1.03     FALSE   
#> 10 fid_10 -0.361   -0.490   0.342  -0.998  -1.28   1.23  -0.828    FALSE   
#> # ℹ 990 more rows
#> # ℹ 11 more variables: m_status <lgl>, mgm_status <lgl>, mgf_status <lgl>,
#> #   mhs1_status <lgl>, mhs2_status <lgl>, o_aoo <dbl>, m_aoo <dbl>,
#> #   mgm_aoo <dbl>, mgf_aoo <dbl>, mhs1_aoo <dbl>, mhs2_aoo <dbl>
#> 
#> $thresholds
#> # A tibble: 6,000 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o     -Inf     3.52
#>  2 fid_2  fid_2_1  o     -Inf     2.76
#>  3 fid_3  fid_3_1  o     -Inf     2.68
#>  4 fid_4  fid_4_1  o     -Inf     3.24
#>  5 fid_5  fid_5_1  o     -Inf     3.52
#>  6 fid_6  fid_6_1  o     -Inf     2.87
#>  7 fid_7  fid_7_1  o        1.31  1.31
#>  8 fid_8  fid_8_1  o     -Inf     2.76
#>  9 fid_9  fid_9_1  o     -Inf     2.68
#> 10 fid_10 fid_10_1 o     -Inf     2.99
#> # ℹ 5,990 more rows
#> 
simulate_under_LTM_single(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE,
h2 = 0.5, n_sim = 500, pop_prev = .05)
#> $sim_obs
#> # A tibble: 500 × 10
#>    fid          m      f      s1 m_status f_status s1_status m_aoo f_aoo s1_aoo
#>    <chr>    <dbl>  <dbl>   <dbl> <lgl>    <lgl>    <lgl>     <dbl> <dbl>  <dbl>
#>  1 fid_1  -0.0137 -0.882  1.07   FALSE    FALSE    FALSE        45    51     25
#>  2 fid_2  -0.684   0.197 -0.433  FALSE    FALSE    FALSE        31    38     12
#>  3 fid_3  -1.46    0.218  0.933  FALSE    FALSE    FALSE        52    43     25
#>  4 fid_4   0.0930 -0.525  1.19   FALSE    FALSE    FALSE        65    56     35
#>  5 fid_5  -0.368  -1.92  -0.803  FALSE    FALSE    FALSE        51    39     21
#>  6 fid_6  -0.514   2.12   1.43   FALSE    TRUE     FALSE        58    55     32
#>  7 fid_7   0.423   1.44  -1.03   FALSE    FALSE    FALSE        41    38     11
#>  8 fid_8   0.674   0.307  1.31   FALSE    FALSE    FALSE        44    39     15
#>  9 fid_9   1.83   -0.469 -0.0785 TRUE     FALSE    FALSE        66    49     22
#> 10 fid_10  0.390  -0.942  0.406  FALSE    FALSE    FALSE        50    50     21
#> # ℹ 490 more rows
#> 
#> $thresholds
#> # A tibble: 1,500 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  m     -Inf     2.48
#>  2 fid_2  fid_2_1  m     -Inf     3.01
#>  3 fid_3  fid_3_1  m     -Inf     2.21
#>  4 fid_4  fid_4_1  m     -Inf     1.84
#>  5 fid_5  fid_5_1  m     -Inf     2.25
#>  6 fid_6  fid_6_1  m     -Inf     2.02
#>  7 fid_7  fid_7_1  m     -Inf     2.63
#>  8 fid_8  fid_8_1  m     -Inf     2.51
#>  9 fid_9  fid_9_1  m        1.83  1.83
#> 10 fid_10 fid_10_1 m     -Inf     2.29
#> # ℹ 1,490 more rows
#> 
simulate_under_LTM_single(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5,
n_sim = 200, pop_prev = 0.05)
#> Warning: Neither fam_vec nor n_fam is specified...
#> $sim_obs
#> # A tibble: 200 × 5
#>    fid          g       o o_status o_aoo
#>    <chr>    <dbl>   <dbl> <lgl>    <dbl>
#>  1 fid_1   0.131  -1.74   FALSE       32
#>  2 fid_2  -0.0343 -1.39   FALSE       22
#>  3 fid_3  -0.932   0.295  FALSE       35
#>  4 fid_4   0.320  -0.0842 FALSE       15
#>  5 fid_5  -0.907  -0.525  FALSE       11
#>  6 fid_6   0.366   1.79   TRUE        68
#>  7 fid_7  -1.43   -2.60   FALSE       17
#>  8 fid_8   0.341   0.466  FALSE       38
#>  9 fid_9   0.459   0.564  FALSE       30
#> 10 fid_10  0.0491 -1.00   FALSE       12
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 200 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o     -Inf     2.97
#>  2 fid_2  fid_2_1  o     -Inf     3.33
#>  3 fid_3  fid_3_1  o     -Inf     2.86
#>  4 fid_4  fid_4_1  o     -Inf     3.57
#>  5 fid_5  fid_5_1  o     -Inf     3.70
#>  6 fid_6  fid_6_1  o        1.79  1.79
#>  7 fid_7  fid_7_1  o     -Inf     3.50
#>  8 fid_8  fid_8_1  o     -Inf     2.75
#>  9 fid_9  fid_9_1  o     -Inf     3.05
#> 10 fid_10 fid_10_1 o     -Inf     3.67
#> # ℹ 190 more rows
#> 
```
