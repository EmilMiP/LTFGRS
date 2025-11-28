# Title Internal Function used to extact input needed from multi-trait graph input for liability estimation

Title Internal Function used to extact input needed from multi-trait
graph input for liability estimation

## Usage

``` r
extract_estimation_info_graph_multi(
  cur_fam_graph,
  fid,
  pid,
  cur_fid,
  h2_vec,
  genetic_corrmat,
  phen_names,
  useMixture,
  add_ind = TRUE
)
```

## Arguments

- cur_fam_graph:

  neightbourhood graph of degree n around proband

- fid:

  Name of column of family ID

- pid:

  Name of column of personal ID

- cur_fid:

  proband ID

- h2_vec:

  vector of heritability values from estimate_liability

- genetic_corrmat:

  genetic correlation matrix as given to estimate_liability

- phen_names:

  vector of phenotype names as given to estimate_liability

- useMixture:

  whether mixture model is used

- add_ind:

  Whether the genetic liability be added. Default is TRUE.

## Value

list with three elements: tbl (tibble with all relevant information),
cov (covariance matrix) estimated through
graph_based_covariance_construction_multi(), and newOrder (order of
individuals in covariance matrix)
