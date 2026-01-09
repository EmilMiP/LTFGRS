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
#>    fid          g       o      m       f      s1     mgm    mgf     pgm     pgf
#>    <chr>    <dbl>   <dbl>  <dbl>   <dbl>   <dbl>   <dbl>  <dbl>   <dbl>   <dbl>
#>  1 fid_1  -1.05   -1.40   -0.337 -1.63   -0.328   1.81   -2.15  -1.13   -0.158 
#>  2 fid_2   1.32    1.17    0.713 -0.747   1.67   -0.664   2.24   1.05    0.0421
#>  3 fid_3   0.408  -0.472  -1.25   1.07   -1.98   -0.287   0.681  0.427  -0.0474
#>  4 fid_4  -0.0738 -0.0329 -0.878 -0.169   0.0600  0.991  -0.442  0.310   0.265 
#>  5 fid_5  -1.36   -1.82    0.619  0.0920 -0.198  -1.67   -0.617  2.04   -0.306 
#>  6 fid_6  -0.437   0.669  -0.314  0.504  -0.960   1.20   -0.876  0.804   1.62  
#>  7 fid_7   0.498   1.52    1.06   0.230   0.488   1.83    2.09  -0.0585 -0.273 
#>  8 fid_8  -0.0524 -0.698  -0.899 -0.502  -0.710   0.410   0.200  1.34   -0.987 
#>  9 fid_9  -0.200  -0.279   0.516  0.319  -0.538  -0.707   0.701 -0.606   2.04  
#> 10 fid_10  0.0708  0.0331 -0.582  0.0307  0.936   0.0437 -0.855 -0.142   0.625 
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
#>  1 fid_1  fid_1_1  o     -Inf     3.35
#>  2 fid_2  fid_2_1  o     -Inf     3.24
#>  3 fid_3  fid_3_1  o     -Inf     3.03
#>  4 fid_4  fid_4_1  o     -Inf     2.87
#>  5 fid_5  fid_5_1  o     -Inf     3.55
#>  6 fid_6  fid_6_1  o     -Inf     2.83
#>  7 fid_7  fid_7_1  o        1.51  1.51
#>  8 fid_8  fid_8_1  o     -Inf     3.14
#>  9 fid_9  fid_9_1  o     -Inf     3.14
#> 10 fid_10 fid_10_1 o     -Inf     3.10
#> # ℹ 7,990 more rows
#> 
simulate_under_LTM_single(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2),
c("m","mgm","mgf","mhs")))
#> $sim_obs
#> # A tibble: 1,000 × 20
#>    fid          g       o       m    mgm    mgf     mhs1    mhs2 o_status
#>    <chr>    <dbl>   <dbl>   <dbl>  <dbl>  <dbl>    <dbl>   <dbl> <lgl>   
#>  1 fid_1  -0.113  -0.0132  0.261  -1.84   0.939  0.00177 -0.450  FALSE   
#>  2 fid_2   0.0784  0.0427  0.691  -0.257  1.60   0.346    0.227  FALSE   
#>  3 fid_3  -0.220  -0.912   0.0968  1.88   0.608 -0.144    1.35   FALSE   
#>  4 fid_4  -0.229  -0.802  -0.0115  0.502  0.837  0.628   -0.0958 FALSE   
#>  5 fid_5  -1.13   -0.468   0.277  -0.345 -1.06   0.204   -0.444  FALSE   
#>  6 fid_6   0.480   0.577   0.898  -2.04   0.481 -0.637    1.23   FALSE   
#>  7 fid_7   0.519   0.356  -0.743  -0.629 -0.115 -1.93    -0.410  FALSE   
#>  8 fid_8   1.04    2.08    0.932   0.915  0.386 -0.261   -1.01   TRUE    
#>  9 fid_9   0.355  -0.200  -0.0574  1.05  -0.131 -0.972   -0.394  FALSE   
#> 10 fid_10 -0.164   0.413  -1.04   -1.51   1.33  -0.925    0.543  FALSE   
#> # ℹ 990 more rows
#> # ℹ 11 more variables: m_status <lgl>, mgm_status <lgl>, mgf_status <lgl>,
#> #   mhs1_status <lgl>, mhs2_status <lgl>, o_aoo <dbl>, m_aoo <dbl>,
#> #   mgm_aoo <dbl>, mgf_aoo <dbl>, mhs1_aoo <dbl>, mhs2_aoo <dbl>
#> 
#> $thresholds
#> # A tibble: 6,000 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o     -Inf     2.68
#>  2 fid_2  fid_2_1  o     -Inf     3.24
#>  3 fid_3  fid_3_1  o     -Inf     3.52
#>  4 fid_4  fid_4_1  o     -Inf     2.87
#>  5 fid_5  fid_5_1  o     -Inf     3.06
#>  6 fid_6  fid_6_1  o     -Inf     2.76
#>  7 fid_7  fid_7_1  o     -Inf     2.68
#>  8 fid_8  fid_8_1  o        2.09  2.09
#>  9 fid_9  fid_9_1  o     -Inf     2.79
#> 10 fid_10 fid_10_1 o     -Inf     3.10
#> # ℹ 5,990 more rows
#> 
simulate_under_LTM_single(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE,
h2 = 0.5, n_sim = 500, pop_prev = .05)
#> $sim_obs
#> # A tibble: 500 × 10
#>    fid         m       f      s1 m_status f_status s1_status m_aoo f_aoo s1_aoo
#>    <chr>   <dbl>   <dbl>   <dbl> <lgl>    <lgl>    <lgl>     <dbl> <dbl>  <dbl>
#>  1 fid_1   0.168 -0.962  -0.0703 FALSE    FALSE    FALSE        65    56     35
#>  2 fid_2  -1.03   1.83    1.24   FALSE    TRUE     FALSE        51    66     21
#>  3 fid_3   0.486 -1.21    0.910  FALSE    FALSE    FALSE        58    60     32
#>  4 fid_4  -1.04  -1.91   -2.38   FALSE    FALSE    FALSE        41    38     11
#>  5 fid_5  -0.912 -0.380  -0.618  FALSE    FALSE    FALSE        44    39     15
#>  6 fid_6  -0.756 -0.0461 -0.691  FALSE    FALSE    FALSE        51    49     22
#>  7 fid_7   0.341 -0.940   0.0704 FALSE    FALSE    FALSE        50    50     21
#>  8 fid_8  -0.649  0.194  -0.980  FALSE    FALSE    FALSE        64    65     40
#>  9 fid_9   1.50  -0.194  -0.110  FALSE    FALSE    FALSE        51    49     24
#> 10 fid_10 -0.435  0.974  -0.347  FALSE    FALSE    FALSE        30    28     10
#> # ℹ 490 more rows
#> 
#> $thresholds
#> # A tibble: 1,500 × 5
#>    fid    indiv_ID role  lower upper
#>    <chr>  <chr>    <chr> <dbl> <dbl>
#>  1 fid_1  fid_1_1  m      -Inf  1.84
#>  2 fid_2  fid_2_1  m      -Inf  2.25
#>  3 fid_3  fid_3_1  m      -Inf  2.02
#>  4 fid_4  fid_4_1  m      -Inf  2.63
#>  5 fid_5  fid_5_1  m      -Inf  2.51
#>  6 fid_6  fid_6_1  m      -Inf  2.25
#>  7 fid_7  fid_7_1  m      -Inf  2.29
#>  8 fid_8  fid_8_1  m      -Inf  1.86
#>  9 fid_9  fid_9_1  m      -Inf  2.25
#> 10 fid_10 fid_10_1 m      -Inf  3.05
#> # ℹ 1,490 more rows
#> 
simulate_under_LTM_single(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5,
n_sim = 200, pop_prev = 0.05)
#> Warning: Neither fam_vec nor n_fam is specified...
#> $sim_obs
#> # A tibble: 200 × 5
#>    fid         g       o o_status o_aoo
#>    <chr>   <dbl>   <dbl> <lgl>    <dbl>
#>  1 fid_1   1.36   1.57   FALSE       15
#>  2 fid_2  -1.04  -1.14   FALSE       11
#>  3 fid_3  -0.367 -0.150  FALSE       24
#>  4 fid_4   1.20   1.13   FALSE       17
#>  5 fid_5  -0.814  0.189  FALSE       38
#>  6 fid_6   1.06   1.99   TRUE        59
#>  7 fid_7   0.806  0.0122 FALSE       12
#>  8 fid_8  -0.298 -0.0825 FALSE       15
#>  9 fid_9   0.656 -0.114  FALSE       39
#> 10 fid_10 -0.206  0.965  FALSE       30
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 200 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o     -Inf     3.57
#>  2 fid_2  fid_2_1  o     -Inf     3.70
#>  3 fid_3  fid_3_1  o     -Inf     3.26
#>  4 fid_4  fid_4_1  o     -Inf     3.50
#>  5 fid_5  fid_5_1  o     -Inf     2.75
#>  6 fid_6  fid_6_1  o        1.99  1.99
#>  7 fid_7  fid_7_1  o     -Inf     3.67
#>  8 fid_8  fid_8_1  o     -Inf     3.57
#>  9 fid_9  fid_9_1  o     -Inf     2.71
#> 10 fid_10 fid_10_1 o     -Inf     3.05
#> # ℹ 190 more rows
#> 
```
