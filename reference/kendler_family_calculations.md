# Title Helper function for Kendler's FGRS

Title Helper function for Kendler's FGRS

## Usage

``` r
kendler_family_calculations(
  tbl,
  cov,
  pid,
  cur_dad_id,
  cur_mom_id,
  env_cor_sib = 1,
  env_cor_f = 1,
  env_cor_m = 1
)
```

## Arguments

- tbl:

  tibble with columns cip, lower, upper, and pid (the personal
  identifier column).

- cov:

  Kinship matrix with proband as first row and column

- pid:

  column name of personal identifier

- cur_dad_id:

  ID of father (not column name, but the actual ID)

- cur_mom_id:

  ID of mother (not column name, but the actual ID)

- env_cor_sib:

  Cohabitation effect, i.e. Factor by which the siblings are weighted.
  Defaults to 1.

- env_cor_f:

  Cohabitation effect, i.e. Factor by which the father is weighted.
  Defaults to 1.

- env_cor_m:

  Cohabitation effect, i.e. Factor by which the mother is weighted.
  Defaults to 1.

## Value

A tibble with family specific values required for Kendler's FGRS
calculation.

## Examples

``` r
# See Vignettes.
```
