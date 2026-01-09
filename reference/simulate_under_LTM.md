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
#>    fid           g       o       m      f     s1     mgm    mgf     pgm     pgf
#>    <chr>     <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>  <dbl>   <dbl>   <dbl>
#>  1 fid_1   0.831   -0.161   0.218   0.734 -1.04   0.245   1.75   2.11   -0.361 
#>  2 fid_2  -0.279   -0.0892  0.0166  0.599  0.114 -0.499  -0.125  0.255   1.19  
#>  3 fid_3  -0.0410  -0.171  -0.839   0.520  0.210 -0.807  -0.480 -1.38    0.248 
#>  4 fid_4   0.990    1.23    0.326   0.493  0.753  0.404   0.282 -0.0224  0.691 
#>  5 fid_5   0.391    1.29    0.452   0.797  0.104 -1.10    0.227  0.203  -0.0584
#>  6 fid_6  -0.00377 -1.21   -0.313  -1.25  -0.109  2.10    1.96  -0.120  -2.22  
#>  7 fid_7  -0.959   -1.10   -0.912  -0.773  0.210 -1.15    0.246 -0.332  -1.50  
#>  8 fid_8  -0.500   -0.248  -0.0187 -0.297 -0.413 -1.38   -0.175 -0.712  -0.348 
#>  9 fid_9  -0.0922   0.185  -1.59   -0.580 -2.26  -1.81   -0.483  0.122  -0.937 
#> 10 fid_10  1.25     0.977   0.855   0.384 -1.15   0.0627  1.51   0.624   1.50  
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
#>  1 fid_1  fid_1_1  o     -Inf     3.28
#>  2 fid_2  fid_2_1  o     -Inf     2.76
#>  3 fid_3  fid_3_1  o     -Inf     3.28
#>  4 fid_4  fid_4_1  o     -Inf     2.59
#>  5 fid_5  fid_5_1  o        1.29  1.29
#>  6 fid_6  fid_6_1  o     -Inf     3.24
#>  7 fid_7  fid_7_1  o     -Inf     2.91
#>  8 fid_8  fid_8_1  o     -Inf     2.68
#>  9 fid_9  fid_9_1  o     -Inf     3.35
#> 10 fid_10 fid_10_1 o     -Inf     2.63
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
#>    fid         g      o      m    mgm     mgf      s1      s2    mhs1    mhs2
#>    <chr>   <dbl>  <dbl>  <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#>  1 fid_1  -0.938 -0.573 -0.257 -1.39   0.527  -0.617  -0.474  -0.590  -0.0331
#>  2 fid_2   1.13   0.857  0.226 -0.568 -0.151   1.14   -0.240   0.410  -1.67  
#>  3 fid_3   0.404  0.792  0.969  0.553  0.577   0.322  -0.0994 -0.0188  0.631 
#>  4 fid_4   0.115 -0.411  0.465  0.615 -1.88    0.135   0.858  -0.740  -0.620 
#>  5 fid_5  -0.572 -1.79  -0.164  0.951 -0.0698  0.0229  0.911   1.37   -0.198 
#>  6 fid_6   0.156 -0.849 -1.48  -0.104 -0.242  -0.0978  0.0133  0.347   0.874 
#>  7 fid_7   0.722 -0.157 -0.195  1.06   1.89   -0.119   1.35   -1.17   -0.468 
#>  8 fid_8  -1.20  -1.99  -2.14   0.221  0.681   0.458  -2.04   -1.35   -0.179 
#>  9 fid_9   0.977  1.61   1.16   1.08   0.0250  0.0835  1.72   -0.787  -0.308 
#> 10 fid_10 -0.421 -0.921  0.243  0.681  0.306   0.227   0.957  -0.600  -0.653 
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
#>  1 fid_1  fid_1_1  o     -Inf     3.24
#>  2 fid_2  fid_2_1  o     -Inf     2.83
#>  3 fid_3  fid_3_1  o     -Inf     2.83
#>  4 fid_4  fid_4_1  o     -Inf     3.48
#>  5 fid_5  fid_5_1  o     -Inf     2.87
#>  6 fid_6  fid_6_1  o     -Inf     2.55
#>  7 fid_7  fid_7_1  o     -Inf     2.59
#>  8 fid_8  fid_8_1  o     -Inf     3.48
#>  9 fid_9  fid_9_1  o        1.62  1.62
#> 10 fid_10 fid_10_1 o     -Inf     3.10
#> # ℹ 7,990 more rows
#> 

simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE,
genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, n_sim = 200)
#> $sim_obs
#> # A tibble: 200 × 10
#>    fid          m      f      s1 m_status f_status s1_status m_aoo f_aoo s1_aoo
#>    <chr>    <dbl>  <dbl>   <dbl> <lgl>    <lgl>    <lgl>     <dbl> <dbl>  <dbl>
#>  1 fid_1  -1.18    0.645 -0.198  FALSE    FALSE    FALSE        55    61     31
#>  2 fid_2  -0.465   1.80  -0.615  FALSE    TRUE     FALSE        45    55     22
#>  3 fid_3   0.173  -1.97  -0.770  FALSE    FALSE    FALSE        32    42     12
#>  4 fid_4   0.902  -0.122  0.100  FALSE    FALSE    FALSE        39    47     20
#>  5 fid_5   0.227   1.86   0.152  FALSE    TRUE     FALSE        67    54     38
#>  6 fid_6  -1.02    0.108  0.0188 FALSE    FALSE    FALSE        53    57     34
#>  7 fid_7  -0.0401 -0.118 -0.953  FALSE    FALSE    FALSE        56    53     35
#>  8 fid_8  -0.548  -0.343 -0.456  FALSE    FALSE    FALSE        32    36     12
#>  9 fid_9  -0.853   0.345  0.0529 FALSE    FALSE    FALSE        60    66     38
#> 10 fid_10 -1.33   -1.93  -1.13   FALSE    FALSE    FALSE        35    39     17
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 600 × 5
#>    fid    indiv_ID role  lower upper
#>    <chr>  <chr>    <chr> <dbl> <dbl>
#>  1 fid_1  fid_1_1  m      -Inf  1.81
#>  2 fid_2  fid_2_1  m      -Inf  2.22
#>  3 fid_3  fid_3_1  m      -Inf  2.76
#>  4 fid_4  fid_4_1  m      -Inf  2.47
#>  5 fid_5  fid_5_1  m      -Inf  1.47
#>  6 fid_6  fid_6_1  m      -Inf  1.89
#>  7 fid_7  fid_7_1  m      -Inf  1.78
#>  8 fid_8  fid_8_1  m      -Inf  2.76
#>  9 fid_9  fid_9_1  m      -Inf  1.64
#> 10 fid_10 fid_10_1 m      -Inf  2.63
#> # ℹ 590 more rows
#> 

simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5,
n_sim = 200, pop_prev = 0.05)
#> Warning: Neither fam_vec nor n_fam is specified...
#> $sim_obs
#> # A tibble: 200 × 5
#>    fid         g      o o_status o_aoo
#>    <chr>   <dbl>  <dbl> <lgl>    <dbl>
#>  1 fid_1  -0.510 -1.47  FALSE       37
#>  2 fid_2  -0.434  0.307 FALSE       40
#>  3 fid_3   0.967  1.10  FALSE       22
#>  4 fid_4   0.149 -0.364 FALSE       33
#>  5 fid_5   0.124 -0.326 FALSE       24
#>  6 fid_6   1.03   2.32  TRUE        49
#>  7 fid_7   1.25   1.58  FALSE       11
#>  8 fid_8   1.06   0.244 FALSE       20
#>  9 fid_9  -0.715 -0.328 FALSE       28
#> 10 fid_10 -1.00  -1.37  FALSE       39
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 200 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o     -Inf     2.79
#>  2 fid_2  fid_2_1  o     -Inf     2.67
#>  3 fid_3  fid_3_1  o     -Inf     3.33
#>  4 fid_4  fid_4_1  o     -Inf     2.94
#>  5 fid_5  fid_5_1  o     -Inf     3.26
#>  6 fid_6  fid_6_1  o        2.32  2.32
#>  7 fid_7  fid_7_1  o     -Inf     3.70
#>  8 fid_8  fid_8_1  o     -Inf     3.40
#>  9 fid_9  fid_9_1  o     -Inf     3.12
#> 10 fid_10 fid_10_1 o     -Inf     2.71
#> # ℹ 190 more rows
#> 
```
