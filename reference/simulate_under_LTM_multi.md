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
#>  1 fid_1         1.32         0.552        0.832      -0.447         1.35  
#>  2 fid_2        -0.256        0.512        1.91       -1.40         -0.544 
#>  3 fid_3         0.385        0.376       -0.297      -1.45         -1.68  
#>  4 fid_4        -0.209       -0.436        0.881      -0.580         0.201 
#>  5 fid_5         0.338       -0.366       -1.64        0.0718       -0.774 
#>  6 fid_6        -0.820       -1.94        -0.150      -0.828         1.34  
#>  7 fid_7         0.440       -0.247        0.638       1.38         -1.51  
#>  8 fid_8         0.148        0.571        0.832       0.224         1.44  
#>  9 fid_9         1.27         1.38        -0.775       1.81          0.0458
#> 10 fid_10        0.597        1.66         1.92       -0.0697       -0.464 
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
#>  1 fid_1        0.296       -0.174         0.738       1.72         -0.0785
#>  2 fid_2        0.0747       0.695         0.688      -1.07          0.290 
#>  3 fid_3        0.941        0.852         0.558      -0.0449        1.64  
#>  4 fid_4        0.146       -0.0504        0.933       0.239        -0.783 
#>  5 fid_5        0.419        1.10         -0.169      -0.488         0.892 
#>  6 fid_6       -0.331       -0.489        -0.317       1.14         -0.271 
#>  7 fid_7       -0.607       -0.599        -0.811      -0.544         1.54  
#>  8 fid_8        0.0910       0.641        -1.78        0.896        -1.44  
#>  9 fid_9       -0.500       -0.564        -0.664      -0.198        -1.39  
#> 10 fid_10      -0.486        0.113        -0.446       0.741         0.588 
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
#>  1 fid_1       -1.88        -1.69         -1.16        -1.76         -2.44 
#>  2 fid_2       -0.277        0.537         1.15        -0.760        -0.537
#>  3 fid_3       -0.857       -1.19         -0.865        0.981        -1.52 
#>  4 fid_4        0.633        1.36         -0.811        0.343         0.907
#>  5 fid_5        0.362        1.30          1.08         0.222         0.806
#>  6 fid_6        0.597       -0.130         0.188        0.780         0.954
#>  7 fid_7        0.106       -0.0109        0.239        1.14          1.44 
#>  8 fid_8        1.03         1.61          0.277       -0.334        -0.329
#>  9 fid_9       -0.442       -0.764         1.30         0.332         0.703
#> 10 fid_10      -0.0873       0.607         0.662       -1.18         -0.518
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
#>  1 fid_1  fid_1_1  o              -Inf                2.68             -Inf
#>  2 fid_2  fid_2_1  o              -Inf                2.95             -Inf
#>  3 fid_3  fid_3_1  o              -Inf                3.17             -Inf
#>  4 fid_4  fid_4_1  o              -Inf                3.38             -Inf
#>  5 fid_5  fid_5_1  o              -Inf                2.79             -Inf
#>  6 fid_6  fid_6_1  o              -Inf                2.63             -Inf
#>  7 fid_7  fid_7_1  o              -Inf                3.17             -Inf
#>  8 fid_8  fid_8_1  o              -Inf                2.55             -Inf
#>  9 fid_9  fid_9_1  o                 1.38             1.38             -Inf
#> 10 fid_10 fid_10_1 o                 1.64             1.64             -Inf
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
#>  1 fid_1         0.582       -0.122       -0.218          0.474         -1.15 
#>  2 fid_2         1.08         1.11         2.50           1.06           2.22 
#>  3 fid_3         0.318       -0.102        0.671          0.922          1.27 
#>  4 fid_4         0.764        1.21         0.518          0.948          2.11 
#>  5 fid_5        -0.666       -0.965       -0.558         -0.720         -0.618
#>  6 fid_6        -1.46        -2.27        -1.79          -0.225         -0.894
#>  7 fid_7         0.561        1.31         0.470          1.11           1.31 
#>  8 fid_8         0.651        1.13         1.33           0.535          0.310
#>  9 fid_9        -0.632        0.133        0.854          0.744         -0.963
#> 10 fid_10       -0.735       -0.785       -0.819          0.610          0.401
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
#>  1 fid_1       -0.322        0.675       -0.316          1.47            0.468
#>  2 fid_2       -0.590       -0.585        0.495         -1.18            0.542
#>  3 fid_3       -0.0328      -0.0405       0.270         -0.0274         -1.03 
#>  4 fid_4        0.530        1.09        -0.723          0.487          -1.19 
#>  5 fid_5        0.763        0.669        1.53           0.490           0.856
#>  6 fid_6       -0.489       -3.04         0.0696         0.357           0.697
#>  7 fid_7       -0.660       -0.803       -0.0668         1.74           -3.72 
#>  8 fid_8        0.759        0.752        0.515         -2.06            0.617
#>  9 fid_9        0.923        0.475        1.92           0.317           1.30 
#> 10 fid_10      -1.10        -0.984       -1.70          -1.23           -0.828
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
#>  1 fid_1         0.811        0.875       1.54           -0.968          0.745
#>  2 fid_2         0.283       -0.621       0.0988          0.136          1.89 
#>  3 fid_3         1.08         1.17        0.942           1.90           0.296
#>  4 fid_4        -0.744       -0.403       0.800           0.779          0.194
#>  5 fid_5         0.193        0.834       1.95            1.96          -0.405
#>  6 fid_6         0.976       -0.181      -1.57           -0.747          1.37 
#>  7 fid_7        -0.433       -0.607       0.810          -0.305         -1.38 
#>  8 fid_8        -1.32        -0.564      -0.978          -0.802         -0.811
#>  9 fid_9         0.122        1.24        0.227          -1.21           2.45 
#> 10 fid_10       -0.442       -0.494       0.890           1.62           1.42 
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
#>  1 fid_1  fid_1_1  o              -Inf                3.42             -Inf
#>  2 fid_2  fid_2_1  o              -Inf                2.95             -Inf
#>  3 fid_3  fid_3_1  o              -Inf                2.99             -Inf
#>  4 fid_4  fid_4_1  o              -Inf                2.83             -Inf
#>  5 fid_5  fid_5_1  o              -Inf                3.24             -Inf
#>  6 fid_6  fid_6_1  o              -Inf                2.55             -Inf
#>  7 fid_7  fid_7_1  o                 1.31             1.31             -Inf
#>  8 fid_8  fid_8_1  o              -Inf                2.95             -Inf
#>  9 fid_9  fid_9_1  o              -Inf                3.42             -Inf
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
#>  1 fid_1        -0.544        1.66          0.825 FALSE              
#>  2 fid_2        -1.21         1.89          0.351 FALSE              
#>  3 fid_3         1.75         0.466         3.35  TRUE               
#>  4 fid_4        -0.592        1.06          1.73  FALSE              
#>  5 fid_5        -0.329       -0.182         1.72  FALSE              
#>  6 fid_6        -1.29        -1.07          0.427 FALSE              
#>  7 fid_7        -0.525       -1.34          0.144 FALSE              
#>  8 fid_8        -0.187       -0.441         0.687 FALSE              
#>  9 fid_9        -1.40         1.55         -1.52  FALSE              
#> 10 fid_10       -0.804       -0.803        -1.52  FALSE              
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
#>  1 fid_1        -0.752       0.129         0.444  FALSE              
#>  2 fid_2         2.26        0.218         0.923  TRUE               
#>  3 fid_3         0.385       0.0167        2.04   FALSE              
#>  4 fid_4         0.849      -0.380        -0.0785 FALSE              
#>  5 fid_5        -0.871       0.705         1.17   FALSE              
#>  6 fid_6        -0.706      -2.17         -0.990  FALSE              
#>  7 fid_7         0.102      -1.28         -0.330  FALSE              
#>  8 fid_8         0.919       0.285        -1.43   FALSE              
#>  9 fid_9        -0.949       0.590        -1.70   FALSE              
#> 10 fid_10        0.610      -0.603        -1.71   FALSE              
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
#>  1 fid_1        -0.994       -0.518        0.218  FALSE              
#>  2 fid_2        -0.351        1.66         0.541  FALSE              
#>  3 fid_3         0.382        0.569        2.63   FALSE              
#>  4 fid_4         1.18         0.188        0.134  FALSE              
#>  5 fid_5         1.01         0.532        2.20   FALSE              
#>  6 fid_6        -1.08        -0.645        0.187  FALSE              
#>  7 fid_7         0.565       -2.77        -1.19   FALSE              
#>  8 fid_8        -0.959        0.303       -0.0441 FALSE              
#>  9 fid_9         1.52         0.190       -0.971  TRUE               
#> 10 fid_10       -0.955       -0.788       -2.12   FALSE              
#> # ℹ 90 more rows
#> # ℹ 5 more variables: f_phenotype3_status <lgl>, s1_phenotype3_status <lgl>,
#> #   m_phenotype3_aoo <dbl>, f_phenotype3_aoo <dbl>, s1_phenotype3_aoo <dbl>
#> 
#> 
#> $thresholds
#> # A tibble: 300 × 9
#>    fid    indiv_ID role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>  <chr>    <chr>            <dbl>            <dbl>            <dbl>
#>  1 fid_1  fid_1_1  m              -Inf                1.71          -Inf   
#>  2 fid_2  fid_2_1  m              -Inf                2.26             2.26
#>  3 fid_3  fid_3_1  m                 1.74             1.74          -Inf   
#>  4 fid_4  fid_4_1  m              -Inf                2.47          -Inf   
#>  5 fid_5  fid_5_1  m              -Inf                2.22          -Inf   
#>  6 fid_6  fid_6_1  m              -Inf                1.68          -Inf   
#>  7 fid_7  fid_7_1  m              -Inf                2.34          -Inf   
#>  8 fid_8  fid_8_1  m              -Inf                1.78          -Inf   
#>  9 fid_9  fid_9_1  m              -Inf                1.71          -Inf   
#> 10 fid_10 fid_10_1 m              -Inf                1.89          -Inf   
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
#>  1 fid_1         0.126       0.495  FALSE                             20
#>  2 fid_2         0.805       0.747  FALSE                             14
#>  3 fid_3         0.848       1.71   TRUE                              58
#>  4 fid_4         0.596       0.0853 FALSE                             38
#>  5 fid_5         0.773       0.551  FALSE                             35
#>  6 fid_6        -0.322       1.26   FALSE                             20
#>  7 fid_7         0.212       0.236  FALSE                             24
#>  8 fid_8        -0.449       0.214  FALSE                             18
#>  9 fid_9         0.161       0.600  FALSE                             15
#> 10 fid_10        1.05        1.46   TRUE                              68
#> # ℹ 140 more rows
#> 
#> 
#> $phenotype2
#> $phenotype2$sim_obs
#> # A tibble: 150 × 5
#>    fid    g_phenotype2 o_phenotype2 o_phenotype2_status o_phenotype2_aoo
#>    <chr>         <dbl>        <dbl> <lgl>                          <dbl>
#>  1 fid_1        -0.133      -0.617  FALSE                             20
#>  2 fid_2         0.804       0.885  FALSE                             14
#>  3 fid_3         0.939       0.765  FALSE                             38
#>  4 fid_4        -2.54       -2.86   FALSE                             38
#>  5 fid_5        -0.684      -0.399  FALSE                             35
#>  6 fid_6         0.456      -0.255  FALSE                             20
#>  7 fid_7        -0.822       0.0780 FALSE                             24
#>  8 fid_8        -0.874      -0.723  FALSE                             18
#>  9 fid_9         0.773       1.28   FALSE                             15
#> 10 fid_10        1.38        0.282  FALSE                             18
#> # ℹ 140 more rows
#> 
#> 
#> $phenotype3
#> $phenotype3$sim_obs
#> # A tibble: 150 × 5
#>    fid    g_phenotype3 o_phenotype3 o_phenotype3_status o_phenotype3_aoo
#>    <chr>         <dbl>        <dbl> <lgl>                          <dbl>
#>  1 fid_1         0.898      1.24    FALSE                             20
#>  2 fid_2        -0.121     -0.00243 FALSE                             14
#>  3 fid_3         1.26       1.57    TRUE                              63
#>  4 fid_4        -0.181     -0.0490  FALSE                             38
#>  5 fid_5         0.358     -0.479   FALSE                             35
#>  6 fid_6        -0.676     -0.265   FALSE                             20
#>  7 fid_7        -1.83      -2.00    FALSE                             24
#>  8 fid_8        -0.489     -0.175   FALSE                             18
#>  9 fid_9         0.530      0.613   FALSE                             15
#> 10 fid_10        0.606      0.0891  FALSE                             18
#> # ℹ 140 more rows
#> 
#> 
#> $thresholds
#> # A tibble: 150 × 9
#>    fid    indiv_ID role  lower_phenotype1 upper_phenotype1 lower_phenotype2
#>    <chr>  <chr>    <chr>            <dbl>            <dbl>            <dbl>
#>  1 fid_1  fid_1_1  o              -Inf                3.21             -Inf
#>  2 fid_2  fid_2_1  o              -Inf                3.42             -Inf
#>  3 fid_3  fid_3_1  o                 1.71             1.71             -Inf
#>  4 fid_4  fid_4_1  o              -Inf                2.51             -Inf
#>  5 fid_5  fid_5_1  o              -Inf                2.63             -Inf
#>  6 fid_6  fid_6_1  o              -Inf                3.21             -Inf
#>  7 fid_7  fid_7_1  o              -Inf                3.06             -Inf
#>  8 fid_8  fid_8_1  o              -Inf                3.28             -Inf
#>  9 fid_9  fid_9_1  o              -Inf                3.38             -Inf
#> 10 fid_10 fid_10_1 o                 1.45             1.45             -Inf
#> # ℹ 140 more rows
#> # ℹ 3 more variables: upper_phenotype2 <dbl>, lower_phenotype3 <dbl>,
#> #   upper_phenotype3 <dbl>
#> 
```
