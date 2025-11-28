# Simulate under the liability threshold model.

`simulate_under_LTM` simulates families and thresholds under the
liability threshold model for a given family structure and a variable
number of phenotypes.Please note that it is not possible to simulate
different family structures.

## Usage

``` r
simulate_under_LTM(
  fam_vec = c("m", "f", "s1", "mgm", "mgf", "pgm", "pgf"),
  n_fam = NULL,
  add_ind = TRUE,
  h2 = 0.5,
  genetic_corrmat = NULL,
  full_corrmat = NULL,
  phen_names = NULL,
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

  Either a number or a numeric vector holding the liability-scale
  heritability(ies) for one or more phenotypes. All entries in `h2` must
  be non-negative. Note that under the liability threshold model, the
  heritabilities must also be at most 1. Defaults to 0.5.

- genetic_corrmat:

  Either `NULL` or a numeric matrix holding the genetic correlations
  between the desired phenotypes. Must be specified, if
  `length(h2)`\\\>0\\, and will be ignored if `h2` is a number. All
  diagonal entries in `genetic_corrmat` must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric. Defaults to `NULL`.

- full_corrmat:

  Either `NULL` or a numeric matrix holding the full correlations
  between the desired phenotypes. Must be specified, if
  `length(h2)`\\\>0\\, and will be ignored if `h2` is a number. All
  diagonal entries in `full_corrmat` must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric. Defaults to `NULL`.

- phen_names:

  Either `NULL` or character vector holding the phenotype names. These
  names will be used to create the row and column names for the
  covariance matrix. Must be specified, if `length(h2)` \\\> 0\\, and
  will be ignored if `h2` is a number. If it is not specified, the names
  will default to phenotype1, phenotype2, etc. Defaults to `NULL`.

- n_sim:

  A positive number representing the number of simulations. Defaults to
  1000.

- pop_prev:

  Either a number or a numeric vector holding the population
  prevalence(s), i.e. the overall prevalence(s) in the population. All
  entries in `pop_prev` must be positive and smaller than 1. Defaults to
  0.1.

## Value

If either `fam_vec` or `n_fam` is used as the argument, if it is of the
required format, if the liability-scale heritability `h2` is a number
satisfying \\0 \leq h^2\\, `n_sim` is a strictly positive number, and
`pop_prev` is a positive number that is at most one, then the output
will be a list containing two tibbles. The first tibble, `sim_obs`,
holds the simulated liabilities, the disease status and the current
age/age-of-onset for all family members in each of the `n_sim` families.
The second tibble, `thresholds`, holds the family identifier, the
personal identifier, the role (specified in fam_vec or n_fam) as well as
the lower and upper thresholds for all individuals in all families. Note
that this tibble has the format required in
[`estimate_liability`](https://emilmip.github.io/LTFGRS/reference/estimate_liability.md).
If either `fam_vec` or `n_fam` is used as the argument and if it is of
the required format, if `genetic_corrmat` and `full_corrmat` are two
numeric and symmetric matrices satisfying that all diagonal entries are
one and that all off-diagonal entries are between -1 and 1, if the
liability-scale heritabilities in `h2_vec` are numbers satisfying \\0
\leq h^2_i\\ for all \\i \in \\1,...,n_pheno\\\\, `n_sim` is a strictly
positive number, and `pop_prev` is a positive numeric vector such that
all entries are at most one, then the output will be a list containing
the following lists. The first outer list, which is named after the
first phenotype in `phen_names`, holds the tibble `sim_obs`, which holds
the simulated liabilities, the disease status and the current
age/age-of-onset for all family members in each of the `n_sim` families
for the first phenotype. As the first outer list, the second outer list,
which is named after the second phenotype in `phen_names`, holds the
tibble `sim_obs`, which holds the simulated liabilities, the disease
status and the current age/age-of-onset for all family members in each
of the `n_sim` families for the second phenotype. There is a list
containing `sim_obs` for each phenotype in `phen_names`. The last list
entry, `thresholds`, holds the family identifier, the personal
identifier, the role (specified in fam_vec or n_fam) as well as the
lower and upper thresholds for all individuals in all families and all
phenotypes. Note that this tibble has the format required in
[`estimate_liability`](https://emilmip.github.io/LTFGRS/reference/estimate_liability.md).
Finally, note that if neither `fam_vec` nor `n_fam` are specified, the
function returns the disease status, the current age/age-of-onset, the
lower and upper thresholds, as well as the personal identifier for a
single individual, namely the individual under consideration (called
`o`). If both `fam_vec` and `n_fam` are defined, the user is asked to '
decide on which of the two vectors to use.

## Details

This function can be used to simulate the case-control status, the
current age and age-of-onset as well as the lower and upper thresholds
for a variable number of phenotypes for all family members in each of
the `n_sim` families. If `h2` is a number, `simulate_under_LTM`
simulates the case- control status, the current age and age-of-onset as
well as thresholds for a single phenotype. However, if `h2` is a numeric
vector, if `genetic_corrmat` and `full_corrmat` are two symmetric
correlation matrices, and if `phen_names` and `pop_prev` are to numeric
vectors holding the phenotype names and the population prevalences,
respectively, then `simulate_under_LTM` simulates the case-control
status, the current age and age-of-onset as well as thresholds for two
or more (correlated) phenotypes. The family members can be specified
using one of two possible formats.

## See also

[`construct_covmat`](https://emilmip.github.io/LTFGRS/reference/construct_covmat.md)
[`simulate_under_LTM_single`](https://emilmip.github.io/LTFGRS/reference/simulate_under_LTM_single.md)
[`simulate_under_LTM_multi`](https://emilmip.github.io/LTFGRS/reference/simulate_under_LTM_multi.md)

## Examples

``` r
simulate_under_LTM()
#> $sim_obs
#> # A tibble: 1,000 × 26
#>    fid          g       o      m        f      s1     mgm    mgf     pgm    pgf
#>    <chr>    <dbl>   <dbl>  <dbl>    <dbl>   <dbl>   <dbl>  <dbl>   <dbl>  <dbl>
#>  1 fid_1  -0.552  -0.583   1.49  -0.172    2.14    0.249   2.70   0.0504  0.415
#>  2 fid_2  -0.0492  0.0242  1.43   0.00244 -0.716  -0.730  -1.44  -0.647   0.259
#>  3 fid_3   0.565   0.133   0.184  0.125   -0.0207 -0.367   1.86  -0.810   0.559
#>  4 fid_4  -1.24   -2.07    0.763 -0.619    0.893   0.0873  2.41   0.194   0.892
#>  5 fid_5  -0.475  -1.07   -0.986  0.0768   1.08   -1.23    0.431  1.74   -0.886
#>  6 fid_6  -0.257  -0.843  -0.758 -1.44     0.751  -0.898   0.271  0.540  -1.42 
#>  7 fid_7   0.802   0.951   1.19  -0.895    0.905  -0.295   0.594 -0.744  -0.462
#>  8 fid_8   0.895   2.11    2.45   0.718    0.822  -0.0944  1.000  0.310   0.473
#>  9 fid_9   0.328   0.458   0.992  0.293    0.0837  0.282  -1.36   2.07   -1.23 
#> 10 fid_10 -0.0797 -0.311  -1.58  -0.762    0.561  -0.225   0.851  0.871  -1.51 
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
#>  1 fid_1  fid_1_1  o     -Inf     2.55
#>  2 fid_2  fid_2_1  o     -Inf     2.63
#>  3 fid_3  fid_3_1  o     -Inf     3.06
#>  4 fid_4  fid_4_1  o     -Inf     3.03
#>  5 fid_5  fid_5_1  o     -Inf     3.17
#>  6 fid_6  fid_6_1  o     -Inf     2.72
#>  7 fid_7  fid_7_1  o     -Inf     3.55
#>  8 fid_8  fid_8_1  o        2.09  2.09
#>  9 fid_9  fid_9_1  o     -Inf     2.76
#> 10 fid_10 fid_10_1 o     -Inf     3.28
#> # ℹ 7,990 more rows
#> 

genetic_corrmat <- matrix(0.4, 3, 3)
diag(genetic_corrmat) <- 1
full_corrmat <- matrix(0.6, 3, 3)
diag(full_corrmat) <- 1

simulate_under_LTM(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2,2),
c("m","mgm","mgf","s","mhs")))
#> $sim_obs
#> # A tibble: 1,000 × 26
#>    fid          g        o       m    mgm    mgf     s1      s2     mhs1    mhs2
#>    <chr>    <dbl>    <dbl>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl>    <dbl>   <dbl>
#>  1 fid_1   0.666   1.13e+0  0.843  -0.445 -1.17  -1.11  -0.0360 -0.00930 -1.19  
#>  2 fid_2   0.444  -3.54e-2 -0.0623 -0.423  0.220  1.70   0.542   0.184   -0.566 
#>  3 fid_3   0.194   9.44e-1 -0.367   0.627 -1.77   0.199  0.700   0.800    0.463 
#>  4 fid_4   0.325   2.80e-1 -0.122  -0.202  0.549  0.414 -0.454   0.751    0.604 
#>  5 fid_5  -1.34   -6.48e-1  0.376  -0.780 -0.641 -0.674 -2.14   -0.348    0.775 
#>  6 fid_6   0.0892  1.30e-1  0.985   1.49  -0.321  0.955 -0.898  -1.52    -0.0170
#>  7 fid_7   0.117   1.20e-1  0.396   0.468  1.08   1.33  -0.450  -0.317    1.20  
#>  8 fid_8   0.795   3.71e-4  0.851  -1.23  -0.425 -0.706 -1.51   -1.89     0.448 
#>  9 fid_9   1.08    1.42e+0 -1.18   -1.03   0.326  0.834  1.47    1.01     0.952 
#> 10 fid_10 -0.230  -4.07e-1  1.20   -0.883 -0.330 -0.538 -0.931   0.360    0.751 
#> # ℹ 990 more rows
#> # ℹ 16 more variables: o_status <lgl>, m_status <lgl>, mgm_status <lgl>,
#> #   mgf_status <lgl>, s1_status <lgl>, s2_status <lgl>, mhs1_status <lgl>,
#> #   mhs2_status <lgl>, o_aoo <dbl>, m_aoo <dbl>, mgm_aoo <dbl>, mgf_aoo <dbl>,
#> #   s1_aoo <dbl>, s2_aoo <dbl>, mhs1_aoo <dbl>, mhs2_aoo <dbl>
#> 
#> $thresholds
#> # A tibble: 8,000 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o     -Inf     2.47
#>  2 fid_2  fid_2_1  o     -Inf     3.38
#>  3 fid_3  fid_3_1  o     -Inf     2.68
#>  4 fid_4  fid_4_1  o     -Inf     2.68
#>  5 fid_5  fid_5_1  o     -Inf     2.51
#>  6 fid_6  fid_6_1  o     -Inf     3.06
#>  7 fid_7  fid_7_1  o     -Inf     2.83
#>  8 fid_8  fid_8_1  o     -Inf     3.55
#>  9 fid_9  fid_9_1  o        1.42  1.42
#> 10 fid_10 fid_10_1 o     -Inf     3.38
#> # ℹ 7,990 more rows
#> 

simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE,
genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, n_sim = 200)
#> $sim_obs
#> # A tibble: 200 × 10
#>    fid          m      f      s1 m_status f_status s1_status m_aoo f_aoo s1_aoo
#>    <chr>    <dbl>  <dbl>   <dbl> <lgl>    <lgl>    <lgl>     <dbl> <dbl>  <dbl>
#>  1 fid_1   1.03    0.926  0.944  FALSE    FALSE    FALSE        59    48     29
#>  2 fid_2   0.871  -0.906  0.0658 FALSE    FALSE    FALSE        60    52     33
#>  3 fid_3   1.02    1.07   2.11   FALSE    FALSE    TRUE         42    47     47
#>  4 fid_4   0.306  -0.529 -0.747  FALSE    FALSE    FALSE        60    58     33
#>  5 fid_5   1.61    2.28   0.555  TRUE     TRUE     FALSE        61    44     18
#>  6 fid_6   0.197  -0.563 -0.455  FALSE    FALSE    FALSE        51    57     27
#>  7 fid_7   0.471  -0.846 -1.09   FALSE    FALSE    FALSE        61    61     38
#>  8 fid_8   0.0529  0.374  0.520  FALSE    FALSE    FALSE        60    70     40
#>  9 fid_9   0.0799  1.24   0.484  FALSE    FALSE    FALSE        50    58     31
#> 10 fid_10 -1.35   -0.456 -1.50   FALSE    FALSE    FALSE        51    42     22
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 600 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  m     -Inf     1.68
#>  2 fid_2  fid_2_1  m     -Inf     1.64
#>  3 fid_3  fid_3_1  m     -Inf     2.34
#>  4 fid_4  fid_4_1  m     -Inf     1.64
#>  5 fid_5  fid_5_1  m        1.62  1.62
#>  6 fid_6  fid_6_1  m     -Inf     1.97
#>  7 fid_7  fid_7_1  m     -Inf     1.62
#>  8 fid_8  fid_8_1  m     -Inf     1.64
#>  9 fid_9  fid_9_1  m     -Inf     2.01
#> 10 fid_10 fid_10_1 m     -Inf     1.97
#> # ℹ 590 more rows
#> 

simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5,
n_sim = 200, pop_prev = 0.05)
#> Warning: Neither fam_vec nor n_fam is specified...
#> $sim_obs
#> # A tibble: 200 × 5
#>    fid         g      o o_status o_aoo
#>    <chr>   <dbl>  <dbl> <lgl>    <dbl>
#>  1 fid_1   0.838  1.87  TRUE        64
#>  2 fid_2  -0.213 -0.137 FALSE       37
#>  3 fid_3  -1.29  -1.45  FALSE       30
#>  4 fid_4   0.780  1.29  FALSE       28
#>  5 fid_5   0.652  0.731 FALSE       35
#>  6 fid_6  -0.211  0.244 FALSE       21
#>  7 fid_7  -0.155  0.310 FALSE       37
#>  8 fid_8   1.87   2.02  TRUE        58
#>  9 fid_9   1.35   2.14  TRUE        54
#> 10 fid_10 -0.592 -1.20  FALSE       33
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 200 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o        1.86  1.86
#>  2 fid_2  fid_2_1  o     -Inf     2.79
#>  3 fid_3  fid_3_1  o     -Inf     3.05
#>  4 fid_4  fid_4_1  o     -Inf     3.12
#>  5 fid_5  fid_5_1  o     -Inf     2.86
#>  6 fid_6  fid_6_1  o     -Inf     3.37
#>  7 fid_7  fid_7_1  o     -Inf     2.79
#>  8 fid_8  fid_8_1  o        2.02  2.02
#>  9 fid_9  fid_9_1  o        2.14  2.14
#> 10 fid_10 fid_10_1 o     -Inf     2.94
#> # ℹ 190 more rows
#> 
```
