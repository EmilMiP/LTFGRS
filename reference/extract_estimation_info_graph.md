# Title Internal Function used to extact input needed from graph input for liability estimation

Title Internal Function used to extact input needed from graph input for
liability estimation

## Usage

``` r
extract_estimation_info_graph(
  cur_fam_graph,
  cur_fid,
  h2,
  pid,
  useMixture,
  add_ind = TRUE
)
```

## Arguments

- cur_fam_graph:

  neightbourhood graph of degree n around proband

- cur_fid:

  proband ID

- h2:

  heritability value from estimate_liability

- pid:

  Name of column of personal ID

- useMixture:

  whether mixture input is returned

- add_ind:

  Whether the genetic liability be added. Default is TRUE.

## Value

list with two elements: tbl (tibble with all relevant information) and
cov (covariance matrix) estimated through
graph_based_covariance_construction()
