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
#>  1 fid_1       -0.421        -0.270       0.0523      -1.14          0.556 
#>  2 fid_2       -0.751        -1.29        0.274       -1.26         -1.12  
#>  3 fid_3       -0.748        -0.545      -0.693        0.540         1.10  
#>  4 fid_4       -0.216         0.455      -1.57        -0.0379        1.17  
#>  5 fid_5        1.21          1.16        0.920        1.20         -0.117 
#>  6 fid_6        0.189         0.686      -0.940        0.0242       -0.133 
#>  7 fid_7       -0.505         1.30        0.808        1.05          1.44  
#>  8 fid_8        0.0510        1.78       -0.209       -0.791        -0.714 
#>  9 fid_9       -0.0823        1.28       -0.408       -0.299         0.730 
#> 10 fid_10      -0.726        -0.933       0.395       -0.106         0.0750
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
#>  1 fid_1         1.68        1.31         0.142         1.29          0.650
#>  2 fid_2        -0.233       0.734       -1.79         -0.167        -0.585
#>  3 fid_3        -0.257      -0.952        0.694         0.967        -0.529
#>  4 fid_4        -0.540      -0.141       -0.984        -0.587        -2.05 
#>  5 fid_5         0.517      -0.0906      -0.0594        0.915        -0.324
#>  6 fid_6        -0.674      -0.704        0.169        -0.571         0.354
#>  7 fid_7         1.40        1.66        -0.171         1.82          1.25 
#>  8 fid_8         0.100       0.954        0.722        -0.141         1.64 
#>  9 fid_9         1.43        1.01         0.0514        0.682         2.90 
#> 10 fid_10        0.692       1.26         0.614         2.29          2.43 
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
#>  1 fid_1        1.46         1.50         0.178        1.63          0.878 
#>  2 fid_2       -0.161       -0.0200       0.246        0.150        -0.0925
#>  3 fid_3        0.577        0.729        0.0663       0.705         1.28  
#>  4 fid_4        0.0959      -0.237       -0.922       -0.749        -0.622 
#>  5 fid_5       -0.960       -1.04        -2.45        -0.694        -0.813 
#>  6 fid_6        0.702       -0.388        0.661        0.0769        0.854 
#>  7 fid_7        0.633        1.61         0.131       -0.305        -0.230 
#>  8 fid_8        0.628        1.68         0.505       -0.874         0.514 
#>  9 fid_9        0.0466       0.601        0.0929       1.26         -1.49  
#> 10 fid_10      -0.855       -1.32         0.0930       0.631        -1.61  
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
#>  1 fid_1  fid_1_1  o              -Inf                3.14             1.31
#>  2 fid_2  fid_2_1  o              -Inf                3.14          -Inf   
#>  3 fid_3  fid_3_1  o              -Inf                2.72          -Inf   
#>  4 fid_4  fid_4_1  o              -Inf                3.21          -Inf   
#>  5 fid_5  fid_5_1  o              -Inf                3.42          -Inf   
#>  6 fid_6  fid_6_1  o              -Inf                2.63          -Inf   
#>  7 fid_7  fid_7_1  o                 1.30             1.30             1.64
#>  8 fid_8  fid_8_1  o                 1.78             1.78          -Inf   
#>  9 fid_9  fid_9_1  o              -Inf                2.55          -Inf   
#> 10 fid_10 fid_10_1 o              -Inf                3.42          -Inf   
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
#>  1 fid_1         0.210        0.335        1.83          -0.243        -0.0676
#>  2 fid_2         1.00         1.93         0.594          0.846        -0.462 
#>  3 fid_3         0.729        0.941        0.398          3.35         -1.10  
#>  4 fid_4         1.47         2.03        -0.313          1.30          1.04  
#>  5 fid_5        -1.31        -1.08         1.19           0.549         0.779 
#>  6 fid_6         0.446       -0.147       -1.31          -0.739        -0.136 
#>  7 fid_7        -0.667       -0.874       -0.846         -0.976         0.277 
#>  8 fid_8        -0.168       -0.987        2.10           0.658         0.131 
#>  9 fid_9         1.07         1.87         1.52           0.247         0.0695
#> 10 fid_10        0.139        0.913       -0.247          0.211         0.0226
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
#>  1 fid_1        0.680         1.10         1.21          -1.49           0.365
#>  2 fid_2        0.450         1.01         2.18           2.57          -0.633
#>  3 fid_3       -0.454        -0.745       -0.635          0.750         -1.39 
#>  4 fid_4       -0.992        -0.492       -0.853         -0.999          1.55 
#>  5 fid_5        0.443         1.08         1.33           0.258          1.94 
#>  6 fid_6        0.987         0.798        0.598          0.342          1.04 
#>  7 fid_7        0.780         0.899       -0.792         -1.17           0.181
#>  8 fid_8       -0.291         0.492        0.237         -0.486         -1.34 
#>  9 fid_9       -0.0233       -0.368       -1.40          -1.42           1.27 
#> 10 fid_10       1.42         -0.322        0.538         -1.15           0.144
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
#>  1 fid_1         0.152     -0.00889       -0.907        -0.900           0.193
#>  2 fid_2        -0.593      0.00699       -1.24          0.655           0.585
#>  3 fid_3         1.13       0.832          1.07          1.71            0.114
#>  4 fid_4         0.349      0.582          1.07          0.492           1.56 
#>  5 fid_5         0.214      1.31          -3.40         -0.701           0.798
#>  6 fid_6        -0.394     -1.97           0.111         1.79           -1.29 
#>  7 fid_7         1.22       0.707          1.47          1.81            0.879
#>  8 fid_8        -1.31      -1.49          -1.13          0.0580         -2.05 
#>  9 fid_9        -0.289     -0.861         -0.227        -1.20           -0.548
#> 10 fid_10       -0.120      0.0878        -0.262         1.26            0.385
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
#>  1 fid_1  fid_1_1  o              -Inf                2.79             -Inf
#>  2 fid_2  fid_2_1  o                 1.93             1.93             -Inf
#>  3 fid_3  fid_3_1  o              -Inf                3.24             -Inf
#>  4 fid_4  fid_4_1  o                 2.01             2.01             -Inf
#>  5 fid_5  fid_5_1  o              -Inf                2.43             -Inf
#>  6 fid_6  fid_6_1  o              -Inf                3.21             -Inf
#>  7 fid_7  fid_7_1  o              -Inf                2.91             -Inf
#>  8 fid_8  fid_8_1  o              -Inf                2.95             -Inf
#>  9 fid_9  fid_9_1  o                 1.89             1.89             -Inf
#> 10 fid_10 fid_10_1 o              -Inf                2.87             -Inf
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
#>  1 fid_1         1.01         0.750       -0.458  FALSE              
#>  2 fid_2         0.993       -0.109        0.551  FALSE              
#>  3 fid_3         0.922        1.78        -0.0302 FALSE              
#>  4 fid_4        -1.30        -0.520       -0.838  FALSE              
#>  5 fid_5        -1.75         0.125       -0.304  FALSE              
#>  6 fid_6        -0.344        0.936        0.825  FALSE              
#>  7 fid_7        -2.12         0.434        0.314  FALSE              
#>  8 fid_8        -0.787        0.325        0.472  FALSE              
#>  9 fid_9         0.252        1.07         0.0813 FALSE              
#> 10 fid_10        0.142        0.114        0.440  FALSE              
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
#>  1 fid_1        1.37         0.399         1.09   TRUE               
#>  2 fid_2        0.0556      -1.20         -1.63   FALSE              
#>  3 fid_3       -0.808       -0.241        -0.147  FALSE              
#>  4 fid_4        0.370       -0.830        -0.987  FALSE              
#>  5 fid_5        1.38         0.943         0.0300 TRUE               
#>  6 fid_6       -0.643        0.838         1.08   FALSE              
#>  7 fid_7       -1.11        -0.581         0.211  FALSE              
#>  8 fid_8       -0.805        0.440         0.0602 FALSE              
#>  9 fid_9        0.370        0.970         0.217  FALSE              
#> 10 fid_10      -0.318       -0.0820       -0.195  FALSE              
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
#>  1 fid_1         2.38         0.437         0.723 TRUE               
#>  2 fid_2        -1.04         0.834         0.367 FALSE              
#>  3 fid_3         0.130        0.349         1.13  FALSE              
#>  4 fid_4         0.963       -1.06          0.570 FALSE              
#>  5 fid_5        -0.877        0.312        -1.32  FALSE              
#>  6 fid_6        -1.06         1.31          0.459 FALSE              
#>  7 fid_7         0.781        0.790         0.538 FALSE              
#>  8 fid_8        -0.160        0.639         0.171 FALSE              
#>  9 fid_9         0.414        0.723         0.910 FALSE              
#> 10 fid_10       -0.768        1.35          1.17  FALSE              
#> # ℹ 90 more rows
#> # ℹ 5 more variables: f_phenotype3_status <lgl>, s1_phenotype3_status <lgl>,
#> #   m_phenotype3_aoo <dbl>, f_phenotype3_aoo <dbl>, s1_phenotype3_aoo <dbl>
#> 
#> 
#> $thresholds
#> # A tibble: 300 × 9
#>    fid    indiv_ID role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>  <chr>    <chr>            <dbl>            <dbl>            <dbl>
#>  1 fid_1  fid_1_1  m                 -Inf             1.74             1.37
#>  2 fid_2  fid_2_1  m                 -Inf             2.30          -Inf   
#>  3 fid_3  fid_3_1  m                 -Inf             2.01          -Inf   
#>  4 fid_4  fid_4_1  m                 -Inf             1.68          -Inf   
#>  5 fid_5  fid_5_1  m                 -Inf             1.62             1.38
#>  6 fid_6  fid_6_1  m                 -Inf             2.13          -Inf   
#>  7 fid_7  fid_7_1  m                 -Inf             1.74          -Inf   
#>  8 fid_8  fid_8_1  m                 -Inf             1.89          -Inf   
#>  9 fid_9  fid_9_1  m                 -Inf             2.34          -Inf   
#> 10 fid_10 fid_10_1 m                 -Inf             2.22          -Inf   
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
#>  1 fid_1         0.130      1.04    FALSE                             34
#>  2 fid_2         0.890      0.194   FALSE                             31
#>  3 fid_3         1.48       1.49    TRUE                              66
#>  4 fid_4        -0.683     -0.461   FALSE                             34
#>  5 fid_5        -0.572     -1.37    FALSE                             39
#>  6 fid_6        -0.653     -0.00127 FALSE                             11
#>  7 fid_7        -0.617     -0.692   FALSE                             28
#>  8 fid_8        -0.388     -0.660   FALSE                             14
#>  9 fid_9        -0.247     -0.833   FALSE                             35
#> 10 fid_10       -0.299     -0.569   FALSE                             35
#> # ℹ 140 more rows
#> 
#> 
#> $phenotype2
#> $phenotype2$sim_obs
#> # A tibble: 150 × 5
#>    fid    g_phenotype2 o_phenotype2 o_phenotype2_status o_phenotype2_aoo
#>    <chr>         <dbl>        <dbl> <lgl>                          <dbl>
#>  1 fid_1        -1.25        -0.838 FALSE                             34
#>  2 fid_2        -0.511        0.553 FALSE                             31
#>  3 fid_3        -0.684       -0.601 FALSE                             26
#>  4 fid_4         1.25         2.13  TRUE                              47
#>  5 fid_5         0.125       -0.431 FALSE                             39
#>  6 fid_6         1.31         0.185 FALSE                             11
#>  7 fid_7        -0.520       -0.522 FALSE                             28
#>  8 fid_8        -0.394       -0.472 FALSE                             14
#>  9 fid_9        -0.971        0.467 FALSE                             35
#> 10 fid_10       -0.437       -1.17  FALSE                             35
#> # ℹ 140 more rows
#> 
#> 
#> $phenotype3
#> $phenotype3$sim_obs
#> # A tibble: 150 × 5
#>    fid    g_phenotype3 o_phenotype3 o_phenotype3_status o_phenotype3_aoo
#>    <chr>         <dbl>        <dbl> <lgl>                          <dbl>
#>  1 fid_1      -0.00973        0.136 FALSE                             34
#>  2 fid_2       0.556          0.641 FALSE                             31
#>  3 fid_3       0.923          0.811 FALSE                             26
#>  4 fid_4       1.37           2.23  TRUE                              45
#>  5 fid_5      -0.213         -0.352 FALSE                             39
#>  6 fid_6       1.13           0.875 FALSE                             11
#>  7 fid_7      -0.686         -0.718 FALSE                             28
#>  8 fid_8      -0.157          0.272 FALSE                             14
#>  9 fid_9       0.505         -0.700 FALSE                             35
#> 10 fid_10     -0.518          0.718 FALSE                             35
#> # ℹ 140 more rows
#> 
#> 
#> $thresholds
#> # A tibble: 150 × 9
#>    fid    indiv_ID role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>  <chr>    <chr>            <dbl>            <dbl>            <dbl>
#>  1 fid_1  fid_1_1  o              -Inf                2.68          -Inf   
#>  2 fid_2  fid_2_1  o              -Inf                2.79          -Inf   
#>  3 fid_3  fid_3_1  o                 1.49             1.49          -Inf   
#>  4 fid_4  fid_4_1  o              -Inf                2.68             2.13
#>  5 fid_5  fid_5_1  o              -Inf                2.47          -Inf   
#>  6 fid_6  fid_6_1  o              -Inf                3.52          -Inf   
#>  7 fid_7  fid_7_1  o              -Inf                2.91          -Inf   
#>  8 fid_8  fid_8_1  o              -Inf                3.42          -Inf   
#>  9 fid_9  fid_9_1  o              -Inf                2.63          -Inf   
#> 10 fid_10 fid_10_1 o              -Inf                2.63          -Inf   
#> # ℹ 140 more rows
#> # ℹ 3 more variables: upper_phenotype2 <dbl>, lower_phenotype3 <dbl>,
#> #   upper_phenotype3 <dbl>
#> 
```
