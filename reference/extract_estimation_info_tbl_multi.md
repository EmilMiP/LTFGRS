# Title Internal Function used to extact input needed from multi-trait tibble input for liability estimation

Title Internal Function used to extact input needed from multi-trait
tibble input for liability estimation

## Usage

``` r
extract_estimation_info_tbl_multi(
  .tbl,
  cur_fid,
  h2,
  fid,
  pid,
  role,
  useMixture,
  phen_names,
  genetic_corrmat,
  full_corrmat,
  add_ind = TRUE
)
```

## Arguments

- .tbl:

  .tbl input from estimate_liability

- cur_fid:

  current family ID being worked on

- h2:

  vector of heritability value from estimate_liability

- fid:

  name of family ID column

- pid:

  name of personal ID column

- role:

  name of role column

- useMixture:

  whether mixture model input is returned

- phen_names:

  vector of phenotype names as given to estimate_liability

- genetic_corrmat:

  genetic correlation matrix as given to estimate_liability

- full_corrmat:

  full correlation matrix as given to estimate_liability

- add_ind:

  Whether the genetic liability be added. Default is TRUE.

## Value

list with two elements: tbl (tibble with all relevant information) and
cov (covariance matrix) estimated through construct_covmat()
