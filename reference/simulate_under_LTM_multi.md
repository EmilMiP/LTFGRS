# Simulate under the liability threshold model (multiple phenotypes).

`simulate_under_LTM_multi` simulates families and thresholds under the
liability threshold model for a given family structure and multiple
phenotypes. Please note that it is not possible to simulate different
family structures.

## Usage

``` r
simulate_under_LTM_multi(
  fam_vec = c("m", "f", "s1", "mgm", "mgf", "pgm", "pgf"),
  n_fam = NULL,
  add_ind = TRUE,
  genetic_corrmat = diag(3),
  full_corrmat = diag(3),
  h2_vec = rep(0.5, 3),
  phen_names = NULL,
  n_sim = 1000,
  pop_prev = rep(0.1, 3)
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

- genetic_corrmat:

  A numeric matrix holding the genetic correlations between the desired
  phenotypes. All diagonal entries must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric. Defaults to `diag(3)`.

- full_corrmat:

  A numeric matrix holding the full correlations between the desired
  phenotypes. All diagonal entries must be equal to one, while all
  off-diagonal entries must be between -1 and 1. In addition, the matrix
  must be symmetric. Defaults to `diag(3)`.

- h2_vec:

  A numeric vector holding the liability-scale heritabilities for a
  number of phenotype. All entries must be non-negative. Note that under
  the liability threshold model, the heritabilities must also be at
  most 1. Defaults to `rep(0.5,3)`.

- phen_names:

  A character vector holding the phenotype names. These names will be
  used to create the row and column names for the covariance matrix. If
  it is not specified, the names will default to phenotype1, phenotype2,
  etc. Defaults to `NULL`.

- n_sim:

  A positive number representing the number of simulations. Defaults to
  1000.

- pop_prev:

  A numeric vector holding the population prevalences, i.e. the overall
  prevalences in the population. All entries in `pop_prev` must be
  positive and smaller than 1. Defaults to `rep(.1,3)`.

## Value

If either `fam_vec` or `n_fam` is used as the argument and if it is of
the required format, if `genetic_corrmat` and `full_corrmat` are two
numeric and symmetric matrices satisfying that all diagonal entries are
one and that all off-diagonal entries are between -1 and 1, if the
liability-scale heritabilities in `h2_vec` are numbers satisfying \\0
\leq h^2_i\\ for all \\i \in \\1,...,n_pheno\\\\, `n_sim` is a strictly
positive number, and `pop_prev` is a positive numeric vector such that
all entries are at most one, then the output will be a list containing
lists for each phenotype. The first outer list, which is named after the
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

## See also

[`construct_covmat`](https://emilmip.github.io/LTFGRS/reference/construct_covmat.md)

## Examples

``` r
simulate_under_LTM_multi()
#> $phenotype1
#> $phenotype1$sim_obs
#> # A tibble: 1,000 × 26
#>    fid    g_phenotype1 o_phenotype1 m_phenotype1 f_phenotype1 s1_phenotype1
#>    <chr>         <dbl>        <dbl>        <dbl>        <dbl>         <dbl>
#>  1 fid_1       -0.753        -0.132       -1.09        -0.404        -0.270
#>  2 fid_2       -0.101         1.10        -0.767       -0.493        -1.14 
#>  3 fid_3       -0.224        -0.891       -0.132       -1.38         -0.474
#>  4 fid_4       -0.0705        0.741        0.210       -0.378         0.641
#>  5 fid_5        0.738         0.133       -0.897        1.42          0.761
#>  6 fid_6        0.316         1.69         0.551        0.305         0.884
#>  7 fid_7       -0.245        -1.29        -1.64        -1.95          1.38 
#>  8 fid_8       -0.321        -1.17         0.462       -0.611         1.99 
#>  9 fid_9        0.687         0.892        0.365       -0.399         1.69 
#> 10 fid_10      -0.696        -0.283        0.270       -1.24         -0.910
#> # ℹ 990 more rows
#> # ℹ 20 more variables: mgm_phenotype1 <dbl>, mgf_phenotype1 <dbl>,
#> #   pgm_phenotype1 <dbl>, pgf_phenotype1 <dbl>, o_phenotype1_status <lgl>,
#> #   m_phenotype1_status <lgl>, f_phenotype1_status <lgl>,
#> #   s1_phenotype1_status <lgl>, mgm_phenotype1_status <lgl>,
#> #   mgf_phenotype1_status <lgl>, pgm_phenotype1_status <lgl>,
#> #   pgf_phenotype1_status <lgl>, o_phenotype1_aoo <dbl>, …
#> 
#> 
#> $phenotype2
#> $phenotype2$sim_obs
#> # A tibble: 1,000 × 26
#>    fid    g_phenotype2 o_phenotype2 m_phenotype2 f_phenotype2 s1_phenotype2
#>    <chr>         <dbl>        <dbl>        <dbl>        <dbl>         <dbl>
#>  1 fid_1        -0.229       -0.542      -0.266        2.28           0.714
#>  2 fid_2         0.327       -0.980      -0.829       -0.223          0.960
#>  3 fid_3         0.208        0.756       0.729       -0.376         -0.903
#>  4 fid_4         0.245        0.590       1.34        -0.538          0.334
#>  5 fid_5         0.766        1.07       -0.0648       0.988         -0.205
#>  6 fid_6         0.793        1.47       -0.964       -0.557         -0.409
#>  7 fid_7        -0.201       -0.486       0.192        1.28           1.16 
#>  8 fid_8         0.631        0.920       0.134       -0.0381         1.23 
#>  9 fid_9        -0.367       -0.522       1.10         1.73           0.481
#> 10 fid_10        1.21         0.744       0.900        0.299          1.33 
#> # ℹ 990 more rows
#> # ℹ 20 more variables: mgm_phenotype2 <dbl>, mgf_phenotype2 <dbl>,
#> #   pgm_phenotype2 <dbl>, pgf_phenotype2 <dbl>, o_phenotype2_status <lgl>,
#> #   m_phenotype2_status <lgl>, f_phenotype2_status <lgl>,
#> #   s1_phenotype2_status <lgl>, mgm_phenotype2_status <lgl>,
#> #   mgf_phenotype2_status <lgl>, pgm_phenotype2_status <lgl>,
#> #   pgf_phenotype2_status <lgl>, o_phenotype2_aoo <dbl>, …
#> 
#> 
#> $phenotype3
#> $phenotype3$sim_obs
#> # A tibble: 1,000 × 26
#>    fid    g_phenotype3 o_phenotype3 m_phenotype3 f_phenotype3 s1_phenotype3
#>    <chr>         <dbl>        <dbl>        <dbl>        <dbl>         <dbl>
#>  1 fid_1       -0.830        -1.66        1.40          1.33         0.971 
#>  2 fid_2       -0.212         0.102       0.0530       -0.330        0.0383
#>  3 fid_3       -0.0909       -0.571      -0.155         0.461        0.476 
#>  4 fid_4       -0.102         1.30       -1.62          0.817       -0.333 
#>  5 fid_5        0.371        -1.33        1.01         -0.848       -0.658 
#>  6 fid_6        0.413         1.28       -1.53          1.68        -0.715 
#>  7 fid_7        0.998         0.270       0.848         0.725        1.85  
#>  8 fid_8       -0.370        -0.859      -0.703         0.307        1.47  
#>  9 fid_9        0.136         0.334       1.50         -0.395        0.841 
#> 10 fid_10      -0.199        -0.630       0.192        -1.20        -1.20  
#> # ℹ 990 more rows
#> # ℹ 20 more variables: mgm_phenotype3 <dbl>, mgf_phenotype3 <dbl>,
#> #   pgm_phenotype3 <dbl>, pgf_phenotype3 <dbl>, o_phenotype3_status <lgl>,
#> #   m_phenotype3_status <lgl>, f_phenotype3_status <lgl>,
#> #   s1_phenotype3_status <lgl>, mgm_phenotype3_status <lgl>,
#> #   mgf_phenotype3_status <lgl>, pgm_phenotype3_status <lgl>,
#> #   pgf_phenotype3_status <lgl>, o_phenotype3_aoo <dbl>, …
#> 
#> 
#> $thresholds
#> # A tibble: 8,000 × 9
#>    fid    indiv_ID role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>  <chr>    <chr>            <dbl>            <dbl>            <dbl>
#>  1 fid_1  fid_1_1  o              -Inf                2.83          -Inf   
#>  2 fid_2  fid_2_1  o              -Inf                2.95          -Inf   
#>  3 fid_3  fid_3_1  o              -Inf                2.76          -Inf   
#>  4 fid_4  fid_4_1  o              -Inf                3.38          -Inf   
#>  5 fid_5  fid_5_1  o              -Inf                3.42          -Inf   
#>  6 fid_6  fid_6_1  o                 1.68             1.68             1.47
#>  7 fid_7  fid_7_1  o              -Inf                3.14          -Inf   
#>  8 fid_8  fid_8_1  o              -Inf                3.14          -Inf   
#>  9 fid_9  fid_9_1  o              -Inf                2.72          -Inf   
#> 10 fid_10 fid_10_1 o              -Inf                3.21          -Inf   
#> # ℹ 7,990 more rows
#> # ℹ 3 more variables: upper_phenotype2 <dbl>, lower_phenotype3 <dbl>,
#> #   upper_phenotype3 <dbl>
#> 

genetic_corrmat <- matrix(0.4, 3, 3)
diag(genetic_corrmat) <- 1
full_corrmat <- matrix(0.6, 3, 3)
diag(full_corrmat) <- 1

simulate_under_LTM_multi(fam_vec = NULL, n_fam = stats::setNames(c(1,1,1,2,2),
c("m","mgm","mgf","s","mhs")))
#> $phenotype1
#> $phenotype1$sim_obs
#> # A tibble: 1,000 × 26
#>    fid    g_phenotype1 o_phenotype1 m_phenotype1 mgm_phenotype1 mgf_phenotype1
#>    <chr>         <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
#>  1 fid_1         0.128     -1.17          0.822          -0.832         0.0560
#>  2 fid_2        -0.351      0.217         1.25            0.854         0.566 
#>  3 fid_3         0.127      0.261        -1.31           -0.257         0.673 
#>  4 fid_4         2.03       2.14          0.0462          0.890         1.77  
#>  5 fid_5        -0.243     -0.248         2.01           -0.415        -1.36  
#>  6 fid_6        -0.674      0.00275      -0.591          -0.674         1.10  
#>  7 fid_7        -0.633     -1.05         -2.22           -0.151        -0.999 
#>  8 fid_8         0.315     -0.0979       -0.278           1.45          0.0554
#>  9 fid_9         0.163     -0.744        -0.649          -0.468         0.626 
#> 10 fid_10        0.269      0.728        -0.194           0.581        -0.280 
#> # ℹ 990 more rows
#> # ℹ 20 more variables: s1_phenotype1 <dbl>, s2_phenotype1 <dbl>,
#> #   mhs1_phenotype1 <dbl>, mhs2_phenotype1 <dbl>, o_phenotype1_status <lgl>,
#> #   m_phenotype1_status <lgl>, mgm_phenotype1_status <lgl>,
#> #   mgf_phenotype1_status <lgl>, s1_phenotype1_status <lgl>,
#> #   s2_phenotype1_status <lgl>, mhs1_phenotype1_status <lgl>,
#> #   mhs2_phenotype1_status <lgl>, o_phenotype1_aoo <dbl>, …
#> 
#> 
#> $phenotype2
#> $phenotype2$sim_obs
#> # A tibble: 1,000 × 26
#>    fid    g_phenotype2 o_phenotype2 m_phenotype2 mgm_phenotype2 mgf_phenotype2
#>    <chr>         <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
#>  1 fid_1        0.105         0.630       0.756          -0.228          0.665
#>  2 fid_2        0.896         1.75        1.12            1.28           0.490
#>  3 fid_3       -0.489         0.525      -0.688          -0.838         -0.273
#>  4 fid_4        0.554        -0.116      -1.61           -1.98          -1.82 
#>  5 fid_5       -1.41         -1.50       -0.0545          1.14          -0.125
#>  6 fid_6       -0.525         0.716       0.354           0.886          0.977
#>  7 fid_7        0.0385        0.218      -1.19            0.340          0.875
#>  8 fid_8       -0.784         0.269      -0.858          -0.558         -0.798
#>  9 fid_9        0.126         0.453      -0.654          -0.870          0.173
#> 10 fid_10       0.289         0.484       0.425          -1.88           3.10 
#> # ℹ 990 more rows
#> # ℹ 20 more variables: s1_phenotype2 <dbl>, s2_phenotype2 <dbl>,
#> #   mhs1_phenotype2 <dbl>, mhs2_phenotype2 <dbl>, o_phenotype2_status <lgl>,
#> #   m_phenotype2_status <lgl>, mgm_phenotype2_status <lgl>,
#> #   mgf_phenotype2_status <lgl>, s1_phenotype2_status <lgl>,
#> #   s2_phenotype2_status <lgl>, mhs1_phenotype2_status <lgl>,
#> #   mhs2_phenotype2_status <lgl>, o_phenotype2_aoo <dbl>, …
#> 
#> 
#> $phenotype3
#> $phenotype3$sim_obs
#> # A tibble: 1,000 × 26
#>    fid    g_phenotype3 o_phenotype3 m_phenotype3 mgm_phenotype3 mgf_phenotype3
#>    <chr>         <dbl>        <dbl>        <dbl>          <dbl>          <dbl>
#>  1 fid_1        -0.353       0.292      -0.00167        -0.0630         0.497 
#>  2 fid_2        -1.33       -1.45        0.0648          0.960         -1.07  
#>  3 fid_3        -0.541       0.0715      0.574           0.225          0.979 
#>  4 fid_4         1.11        1.60        0.720          -3.27           0.519 
#>  5 fid_5         0.417      -1.21       -0.855          -0.712          0.170 
#>  6 fid_6         0.802       0.997       1.14            0.604          0.721 
#>  7 fid_7         0.554       0.772       1.31            0.0657         1.83  
#>  8 fid_8        -0.615      -0.760      -1.43           -0.626         -1.47  
#>  9 fid_9         0.432       0.0293      0.262          -0.331         -0.0442
#> 10 fid_10       -0.330      -0.336       0.509           0.474         -0.198 
#> # ℹ 990 more rows
#> # ℹ 20 more variables: s1_phenotype3 <dbl>, s2_phenotype3 <dbl>,
#> #   mhs1_phenotype3 <dbl>, mhs2_phenotype3 <dbl>, o_phenotype3_status <lgl>,
#> #   m_phenotype3_status <lgl>, mgm_phenotype3_status <lgl>,
#> #   mgf_phenotype3_status <lgl>, s1_phenotype3_status <lgl>,
#> #   s2_phenotype3_status <lgl>, mhs1_phenotype3_status <lgl>,
#> #   mhs2_phenotype3_status <lgl>, o_phenotype3_aoo <dbl>, …
#> 
#> 
#> $thresholds
#> # A tibble: 8,000 × 9
#>    fid    indiv_ID role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>  <chr>    <chr>            <dbl>            <dbl>            <dbl>
#>  1 fid_1  fid_1_1  o              -Inf                3.38          -Inf   
#>  2 fid_2  fid_2_1  o              -Inf                2.91             1.74
#>  3 fid_3  fid_3_1  o              -Inf                2.79          -Inf   
#>  4 fid_4  fid_4_1  o                 2.13             2.13          -Inf   
#>  5 fid_5  fid_5_1  o              -Inf                3.24          -Inf   
#>  6 fid_6  fid_6_1  o              -Inf                3.31          -Inf   
#>  7 fid_7  fid_7_1  o              -Inf                2.91          -Inf   
#>  8 fid_8  fid_8_1  o              -Inf                2.79          -Inf   
#>  9 fid_9  fid_9_1  o              -Inf                3.17          -Inf   
#> 10 fid_10 fid_10_1 o              -Inf                3.24          -Inf   
#> # ℹ 7,990 more rows
#> # ℹ 3 more variables: upper_phenotype2 <dbl>, lower_phenotype3 <dbl>,
#> #   upper_phenotype3 <dbl>
#> 

simulate_under_LTM_multi(fam_vec = c("m","f","s1"), add_ind = FALSE,
genetic_corrmat = genetic_corrmat, full_corrmat = full_corrmat, n_sim = 100)
#> $phenotype1
#> $phenotype1$sim_obs
#> # A tibble: 100 × 10
#>    fid    m_phenotype1 f_phenotype1 s1_phenotype1 m_phenotype1_status
#>    <chr>         <dbl>        <dbl>         <dbl> <lgl>              
#>  1 fid_1        -0.834       1.77         -0.887  FALSE              
#>  2 fid_2        -0.993       0.682        -0.849  FALSE              
#>  3 fid_3        -0.914       0.363         0.156  FALSE              
#>  4 fid_4        -0.869       1.04          0.620  FALSE              
#>  5 fid_5         1.82       -1.76         -0.339  TRUE               
#>  6 fid_6        -0.523      -1.22          0.661  FALSE              
#>  7 fid_7         1.69        0.0377        0.323  TRUE               
#>  8 fid_8         0.519      -0.181         0.488  FALSE              
#>  9 fid_9        -0.787       0.195        -0.0356 FALSE              
#> 10 fid_10       -1.23       -1.43         -1.93   FALSE              
#> # ℹ 90 more rows
#> # ℹ 5 more variables: f_phenotype1_status <lgl>, s1_phenotype1_status <lgl>,
#> #   m_phenotype1_aoo <dbl>, f_phenotype1_aoo <dbl>, s1_phenotype1_aoo <dbl>
#> 
#> 
#> $phenotype2
#> $phenotype2$sim_obs
#> # A tibble: 100 × 10
#>    fid    m_phenotype2 f_phenotype2 s1_phenotype2 m_phenotype2_status
#>    <chr>         <dbl>        <dbl>         <dbl> <lgl>              
#>  1 fid_1         0.841        1.27          0.330 FALSE              
#>  2 fid_2        -2.21        -0.740        -1.61  FALSE              
#>  3 fid_3        -0.585       -0.163         0.449 FALSE              
#>  4 fid_4         0.706       -1.84         -1.17  FALSE              
#>  5 fid_5         0.729       -1.46          1.13  FALSE              
#>  6 fid_6        -0.264       -1.87         -0.432 FALSE              
#>  7 fid_7         2.34        -0.131         0.490 TRUE               
#>  8 fid_8        -0.125       -0.156         0.672 FALSE              
#>  9 fid_9         0.316       -0.383         0.517 FALSE              
#> 10 fid_10       -0.955       -1.40         -1.04  FALSE              
#> # ℹ 90 more rows
#> # ℹ 5 more variables: f_phenotype2_status <lgl>, s1_phenotype2_status <lgl>,
#> #   m_phenotype2_aoo <dbl>, f_phenotype2_aoo <dbl>, s1_phenotype2_aoo <dbl>
#> 
#> 
#> $phenotype3
#> $phenotype3$sim_obs
#> # A tibble: 100 × 10
#>    fid    m_phenotype3 f_phenotype3 s1_phenotype3 m_phenotype3_status
#>    <chr>         <dbl>        <dbl>         <dbl> <lgl>              
#>  1 fid_1         0.121        1.92         -0.906 FALSE              
#>  2 fid_2        -0.299        0.266        -1.02  FALSE              
#>  3 fid_3        -2.62        -0.952        -0.592 FALSE              
#>  4 fid_4         0.172       -0.149        -1.03  FALSE              
#>  5 fid_5        -0.238       -1.21         -0.407 FALSE              
#>  6 fid_6        -0.342       -1.89          0.391 FALSE              
#>  7 fid_7        -0.280       -0.702        -1.00  FALSE              
#>  8 fid_8         0.115       -0.506         0.729 FALSE              
#>  9 fid_9         1.24         0.919         0.515 FALSE              
#> 10 fid_10        0.764       -0.642        -0.795 FALSE              
#> # ℹ 90 more rows
#> # ℹ 5 more variables: f_phenotype3_status <lgl>, s1_phenotype3_status <lgl>,
#> #   m_phenotype3_aoo <dbl>, f_phenotype3_aoo <dbl>, s1_phenotype3_aoo <dbl>
#> 
#> 
#> $thresholds
#> # A tibble: 300 × 9
#>    fid    indiv_ID role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>  <chr>    <chr>            <dbl>            <dbl>            <dbl>
#>  1 fid_1  fid_1_1  m              -Inf                2.18          -Inf   
#>  2 fid_2  fid_2_1  m              -Inf                2.09          -Inf   
#>  3 fid_3  fid_3_1  m              -Inf                1.71          -Inf   
#>  4 fid_4  fid_4_1  m              -Inf                1.64          -Inf   
#>  5 fid_5  fid_5_1  m                 1.81             1.81          -Inf   
#>  6 fid_6  fid_6_1  m              -Inf                1.74          -Inf   
#>  7 fid_7  fid_7_1  m                 1.71             1.71             2.34
#>  8 fid_8  fid_8_1  m              -Inf                2.01          -Inf   
#>  9 fid_9  fid_9_1  m              -Inf                1.68          -Inf   
#> 10 fid_10 fid_10_1 m              -Inf                1.62          -Inf   
#> # ℹ 290 more rows
#> # ℹ 3 more variables: upper_phenotype2 <dbl>, lower_phenotype3 <dbl>,
#> #   upper_phenotype3 <dbl>
#> 

simulate_under_LTM_multi(fam_vec = c(), n_fam = NULL, add_ind = TRUE, n_sim = 150)
#> $phenotype1
#> $phenotype1$sim_obs
#> # A tibble: 150 × 5
#>    fid    g_phenotype1 o_phenotype1 o_phenotype1_status o_phenotype1_aoo
#>    <chr>         <dbl>        <dbl> <lgl>                          <dbl>
#>  1 fid_1       -0.449       0.214   FALSE                             37
#>  2 fid_2        0.161       0.600   FALSE                             20
#>  3 fid_3        1.05        1.46    TRUE                              68
#>  4 fid_4        0.751      -0.0429  FALSE                             21
#>  5 fid_5        0.911       1.40    TRUE                              72
#>  6 fid_6       -0.463      -2.00    FALSE                             31
#>  7 fid_7        0.690       0.163   FALSE                             26
#>  8 fid_8       -0.329       0.159   FALSE                             34
#>  9 fid_9       -0.0705     -0.00849 FALSE                             39
#> 10 fid_10       0.299       0.237   FALSE                             11
#> # ℹ 140 more rows
#> 
#> 
#> $phenotype2
#> $phenotype2$sim_obs
#> # A tibble: 150 × 5
#>    fid    g_phenotype2 o_phenotype2 o_phenotype2_status o_phenotype2_aoo
#>    <chr>         <dbl>        <dbl> <lgl>                          <dbl>
#>  1 fid_1      -0.874         -0.723 FALSE                             37
#>  2 fid_2       0.773          1.28  FALSE                             20
#>  3 fid_3       1.38           0.282 FALSE                             32
#>  4 fid_4       0.873          0.373 FALSE                             21
#>  5 fid_5      -0.547         -0.691 FALSE                             34
#>  6 fid_6      -0.0146         0.255 FALSE                             31
#>  7 fid_7      -0.812         -0.282 FALSE                             26
#>  8 fid_8      -0.893         -0.307 FALSE                             34
#>  9 fid_9       0.449          0.222 FALSE                             39
#> 10 fid_10      0.00115       -1.06  FALSE                             11
#> # ℹ 140 more rows
#> 
#> 
#> $phenotype3
#> $phenotype3$sim_obs
#> # A tibble: 150 × 5
#>    fid    g_phenotype3 o_phenotype3 o_phenotype3_status o_phenotype3_aoo
#>    <chr>         <dbl>        <dbl> <lgl>                          <dbl>
#>  1 fid_1        -0.489      -0.175  FALSE                             37
#>  2 fid_2         0.530       0.613  FALSE                             20
#>  3 fid_3         0.606       0.0891 FALSE                             32
#>  4 fid_4         0.538       0.608  FALSE                             21
#>  5 fid_5        -0.187      -0.0331 FALSE                             34
#>  6 fid_6         0.450       0.218  FALSE                             31
#>  7 fid_7        -0.135      -0.428  FALSE                             26
#>  8 fid_8        -0.124       0.757  FALSE                             34
#>  9 fid_9        -1.04       -1.15   FALSE                             39
#> 10 fid_10       -0.421      -0.461  FALSE                             11
#> # ℹ 140 more rows
#> 
#> 
#> $thresholds
#> # A tibble: 150 × 9
#>    fid    indiv_ID role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>  <chr>    <chr>            <dbl>            <dbl>            <dbl>
#>  1 fid_1  fid_1_1  o              -Inf                2.55             -Inf
#>  2 fid_2  fid_2_1  o              -Inf                3.21             -Inf
#>  3 fid_3  fid_3_1  o                 1.45             1.45             -Inf
#>  4 fid_4  fid_4_1  o              -Inf                3.17             -Inf
#>  5 fid_5  fid_5_1  o                 1.39             1.39             -Inf
#>  6 fid_6  fid_6_1  o              -Inf                2.79             -Inf
#>  7 fid_7  fid_7_1  o              -Inf                2.99             -Inf
#>  8 fid_8  fid_8_1  o              -Inf                2.68             -Inf
#>  9 fid_9  fid_9_1  o              -Inf                2.47             -Inf
#> 10 fid_10 fid_10_1 o              -Inf                3.52             -Inf
#> # ℹ 140 more rows
#> # ℹ 3 more variables: upper_phenotype2 <dbl>, lower_phenotype3 <dbl>,
#> #   upper_phenotype3 <dbl>
#> 
```
