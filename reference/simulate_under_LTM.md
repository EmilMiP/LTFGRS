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
#>    fid         g       o      m        f      s1      mgm    mgf     pgm     pgf
#>    <chr>   <dbl>   <dbl>  <dbl>    <dbl>   <dbl>    <dbl>  <dbl>   <dbl>   <dbl>
#>  1 fid_1 -0.860   0.157  -1.55   1.24    -0.558  -0.00886 -2.53   0.590  -0.780 
#>  2 fid_2 -0.474   0.646   0.865  0.464    0.812  -0.313    1.10  -0.406   0.885 
#>  3 fid_3  1.58    2.97    0.225  1.34     0.917  -0.506    1.90   0.416   0.343 
#>  4 fid_4  1.70    1.52    1.53   0.600    0.334   1.41    -0.328  1.35    0.0340
#>  5 fid_5 -0.552  -0.583   1.49  -0.172    2.14    0.249    2.70   0.0504  0.415 
#>  6 fid_6 -0.0492  0.0242  1.43   0.00244 -0.716  -0.730   -1.44  -0.647   0.259 
#>  7 fid_7  0.565   0.133   0.184  0.125   -0.0207 -0.367    1.86  -0.810   0.559 
#>  8 fid_8 -1.24   -2.07    0.763 -0.619    0.893   0.0873   2.41   0.194   0.892 
#>  9 fid_9 -0.475  -1.07   -0.986  0.0768   1.08   -1.23     0.431  1.74   -0.886 
#> 10 fid_… -0.257  -0.843  -0.758 -1.44     0.751  -0.898    0.271  0.540  -1.42  
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
#>  1 fid_1  fid_1_1  o     -Inf     3.38
#>  2 fid_2  fid_2_1  o     -Inf     3.45
#>  3 fid_3  fid_3_1  o        2.99  2.99
#>  4 fid_4  fid_4_1  o        1.51  1.51
#>  5 fid_5  fid_5_1  o     -Inf     2.68
#>  6 fid_6  fid_6_1  o     -Inf     3.55
#>  7 fid_7  fid_7_1  o     -Inf     3.35
#>  8 fid_8  fid_8_1  o     -Inf     2.68
#>  9 fid_9  fid_9_1  o     -Inf     2.95
#> 10 fid_10 fid_10_1 o     -Inf     3.24
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
#>    fid          g       o      m    mgm     mgf     s1     s2    mhs1     mhs2
#>    <chr>    <dbl>   <dbl>  <dbl>  <dbl>   <dbl>  <dbl>  <dbl>   <dbl>    <dbl>
#>  1 fid_1   0.597  -0.0736  0.593  1.54  -0.0975  2.15  -0.622 -0.312  -0.688  
#>  2 fid_2   0.832   0.795   0.170 -1.32  -0.0514  0.145  0.222  0.960  -0.0251 
#>  3 fid_3  -0.415  -1.25   -1.26   0.592 -1.93    0.646 -1.55  -1.73   -0.00565
#>  4 fid_4  -0.437  -0.176   0.127  0.416 -2.05   -0.108  0.154 -1.71   -1.06   
#>  5 fid_5   0.0853  0.739   1.95   1.89   0.203   0.729  0.516  1.32    1.41   
#>  6 fid_6   0.451  -0.0504 -0.412  1.03  -0.383   1.90   0.538  0.322   1.82   
#>  7 fid_7   0.777  -0.190   0.466  1.26   0.809   0.885  0.855  0.514   0.722  
#>  8 fid_8   0.746   0.709  -0.512 -0.568  1.73    0.929  1.02   1.60    1.26   
#>  9 fid_9  -0.933  -0.470   0.644 -0.742  1.27   -0.323  0.564  0.0259 -0.945  
#> 10 fid_10 -0.321  -0.428   0.319  0.892 -0.701  -0.309 -0.411  0.120   0.688  
#> # ℹ 990 more rows
#> # ℹ 16 more variables: o_status <lgl>, m_status <lgl>, mgm_status <lgl>,
#> #   mgf_status <lgl>, s1_status <lgl>, s2_status <lgl>, mhs1_status <lgl>,
#> #   mhs2_status <lgl>, o_aoo <dbl>, m_aoo <dbl>, mgm_aoo <dbl>, mgf_aoo <dbl>,
#> #   s1_aoo <dbl>, s2_aoo <dbl>, mhs1_aoo <dbl>, mhs2_aoo <dbl>
#> 
#> $thresholds
#> # A tibble: 8,000 × 5
#>    fid    indiv_ID role  lower upper
#>    <chr>  <chr>    <chr> <dbl> <dbl>
#>  1 fid_1  fid_1_1  o      -Inf  2.95
#>  2 fid_2  fid_2_1  o      -Inf  2.43
#>  3 fid_3  fid_3_1  o      -Inf  2.83
#>  4 fid_4  fid_4_1  o      -Inf  2.43
#>  5 fid_5  fid_5_1  o      -Inf  2.59
#>  6 fid_6  fid_6_1  o      -Inf  2.72
#>  7 fid_7  fid_7_1  o      -Inf  3.52
#>  8 fid_8  fid_8_1  o      -Inf  3.28
#>  9 fid_9  fid_9_1  o      -Inf  2.47
#> 10 fid_10 fid_10_1 o      -Inf  3.42
#> # ℹ 7,990 more rows
#> 

simulate_under_LTM(fam_vec = c("m","f","s1"), n_fam = NULL, add_ind = FALSE,
genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, n_sim = 200)
#> $sim_obs
#> # A tibble: 200 × 10
#>    fid          m       f     s1 m_status f_status s1_status m_aoo f_aoo s1_aoo
#>    <chr>    <dbl>   <dbl>  <dbl> <lgl>    <lgl>    <lgl>     <dbl> <dbl>  <dbl>
#>  1 fid_1   0.971   1.48   -1.11  FALSE    TRUE     FALSE        48    66     21
#>  2 fid_2  -1.35   -0.244  -0.397 FALSE    FALSE    FALSE        65    58     35
#>  3 fid_3   0.183   1.14    0.694 FALSE    FALSE    FALSE        31    41     12
#>  4 fid_4  -1.08    0.318  -0.455 FALSE    FALSE    FALSE        66    64     39
#>  5 fid_5  -0.434   1.98    0.923 FALSE    TRUE     FALSE        51    51     27
#>  6 fid_6  -0.0308  1.24   -1.67  FALSE    FALSE    FALSE        41    36     14
#>  7 fid_7   0.0450  1.94    0.556 FALSE    TRUE     FALSE        35    52     17
#>  8 fid_8  -0.446  -0.0383 -0.843 FALSE    FALSE    FALSE        58    63     35
#>  9 fid_9  -0.485   0.485  -0.221 FALSE    FALSE    FALSE        41    39     14
#> 10 fid_10 -0.230  -0.0316 -1.08  FALSE    FALSE    FALSE        39    33     15
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 600 × 5
#>    fid    indiv_ID role  lower upper
#>    <chr>  <chr>    <chr> <dbl> <dbl>
#>  1 fid_1  fid_1_1  m      -Inf  2.09
#>  2 fid_2  fid_2_1  m      -Inf  1.51
#>  3 fid_3  fid_3_1  m      -Inf  2.79
#>  4 fid_4  fid_4_1  m      -Inf  1.49
#>  5 fid_5  fid_5_1  m      -Inf  1.97
#>  6 fid_6  fid_6_1  m      -Inf  2.39
#>  7 fid_7  fid_7_1  m      -Inf  2.63
#>  8 fid_8  fid_8_1  m      -Inf  1.71
#>  9 fid_9  fid_9_1  m      -Inf  2.39
#> 10 fid_10 fid_10_1 m      -Inf  2.47
#> # ℹ 590 more rows
#> 

simulate_under_LTM(fam_vec = c(), n_fam = NULL, add_ind = TRUE, h2 = 0.5,
n_sim = 200, pop_prev = 0.05)
#> Warning: Neither fam_vec nor n_fam is specified...
#> $sim_obs
#> # A tibble: 200 × 5
#>    fid         g       o o_status o_aoo
#>    <chr>   <dbl>   <dbl> <lgl>    <dbl>
#>  1 fid_1   1.31   0.0869 FALSE       30
#>  2 fid_2   0.280 -0.648  FALSE       13
#>  3 fid_3   0.406  0.487  FALSE       15
#>  4 fid_4   0.873  0.103  FALSE       21
#>  5 fid_5   0.578  0.407  FALSE       21
#>  6 fid_6   0.108  0.767  FALSE       40
#>  7 fid_7   0.603 -1.09   FALSE       19
#>  8 fid_8   0.133 -0.936  FALSE       23
#>  9 fid_9  -0.301 -0.173  FALSE       13
#> 10 fid_10  1.72   2.58   TRUE        42
#> # ℹ 190 more rows
#> 
#> $thresholds
#> # A tibble: 200 × 5
#>    fid    indiv_ID role    lower upper
#>    <chr>  <chr>    <chr>   <dbl> <dbl>
#>  1 fid_1  fid_1_1  o     -Inf     3.05
#>  2 fid_2  fid_2_1  o     -Inf     3.63
#>  3 fid_3  fid_3_1  o     -Inf     3.57
#>  4 fid_4  fid_4_1  o     -Inf     3.37
#>  5 fid_5  fid_5_1  o     -Inf     3.37
#>  6 fid_6  fid_6_1  o     -Inf     2.67
#>  7 fid_7  fid_7_1  o     -Inf     3.44
#>  8 fid_8  fid_8_1  o     -Inf     3.30
#>  9 fid_9  fid_9_1  o     -Inf     3.63
#> 10 fid_10 fid_10_1 o        2.59  2.59
#> # ℹ 190 more rows
#> 
```
