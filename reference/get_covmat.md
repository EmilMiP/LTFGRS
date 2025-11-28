# Construct kinship matrix from graph

construct the kinship matrix from a graph representation of a family,
centered on an index person (proband).

## Usage

``` r
get_covmat(fam_graph, h2, index_id = NA, add_ind = TRUE, fix_diag = TRUE)
```

## Arguments

- fam_graph:

  graph.

- h2:

  heritability.

- index_id:

  proband id. Only used in conjuction with add_ind = TRUE.

- add_ind:

  add genetic liability to the kinship matrix. Defaults to true.

- fix_diag:

  Whether to set diagonal to 1 for all entries except for the genetic
  liability.

## Value

A kinship matrix.

## Examples

``` r
fam <- data.frame(
i = c(1, 2, 3, 4),
f = c(3, 0, 4, 0),
m = c(2, 0, 0, 0)
)

thresholds <- data.frame(
  i = c(1, 2, 3, 4),
  lower = c(-Inf, -Inf, 0.8, 0.7),
  upper = c(0.8, 0.8, 0.8, 0.7)
)

graph <- prepare_graph(fam, icol = "i", fcol = "f", mcol = "m", node_attributes = thresholds)

get_covmat(graph, h2 = 0.5, index_id = "1")
#>         1    2    3     4   1_g
#> 1   1.000 0.25 0.25 0.125 0.500
#> 2   0.250 1.00 0.00 0.000 0.250
#> 3   0.250 0.00 1.00 0.250 0.250
#> 4   0.125 0.00 0.25 1.000 0.125
#> 1_g 0.500 0.25 0.25 0.125 0.500
```
