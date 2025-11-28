# Title Internal Function used to extact input needed for liability estimation

Title Internal Function used to extact input needed for liability
estimation

## Usage

``` r
extract_estimation_info_tbl(
  .tbl,
  cur_fid,
  h2,
  fid,
  pid,
  role,
  useMixture,
  add_ind = TRUE
)
```

## Arguments

- .tbl:

  .tbl input from estimate_liability

- cur_fid:

  current family ID being worked on

- h2:

  heritability value from estimate_liability

- fid:

  name of family ID column

- pid:

  name of personal ID column

- role:

  name of role column

- useMixture:

  whether mixture model input is returned

- add_ind:

  Whether the genetic liability be added. Default is TRUE.

## Value

list with two elements: tbl (tibble with all relevant information) and
cov (covariance matrix) estimated through construct_covmat()
