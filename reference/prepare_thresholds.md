# Calculate (personalised) thresholds based on CIPs.

This function prepares input for `estimate_liability` by calculating
thresholds based on stratified cumulative incidence proportions (CIPs)
with options for interpolation for ages between CIP values. Given a
tibble with families and family members and (stratified) CIPs,
personalised thresholds will be calculated for each individual present
in `.tbl`. An individual may be in multiple families, but only once in
the same family.

## Usage

``` r
prepare_thresholds(
  .tbl,
  CIP,
  age_col,
  CIP_merge_columns = c("sex", "birth_year", "age"),
  CIP_cip_col = "cip",
  Kpop = "useMax",
  K_i_col = "K_i",
  K_pop_col = "K_pop",
  status_col = "status",
  thr_col = "thr",
  lower_col = "lower",
  upper_col = "upper",
  lower_equal_upper = FALSE,
  personal_thr = FALSE,
  interpolation = NULL,
  bst.params = list(max_depth = 10, base_score = 0, nthread = 4, min_child_weight = 10),
  min_CIP_value = 1e-05,
  xgboost_itr = 30
)
```

## Arguments

- .tbl:

  Tibble with family and personal id columns, as well as
  CIP_merge_columns and status.

- CIP:

  Tibble with population representative cumulative incidence
  proportions. CIP must contain columns from `CIP_merge_columns` and
  `cIP_cip_col`.

- age_col:

  Name of column with age at the end of follow-up or age at diagnosis
  for cases.

- CIP_merge_columns:

  The columns the CIPs are subset by, e.g. CIPs by birth_year, sex.

- CIP_cip_col:

  Name of column with CIP values.

- Kpop:

  Takes either "useMax" to use the maximum value in the CIP strata as
  population prevalence, or a tibble with population prevalence values
  based on other information. If a tibble is provided, it must contain
  columns from `.tbl` and a column named "K_pop" with population
  prevalence values. Defaults to "UseMax".

- K_i_col:

  Name of column to create for individual cumulative incidence
  proportion. Defaults to "K_i".

- K_pop_col:

  Name of column to create for population prevalence within an
  individuals CIP strata. Defaults to "K_pop".

- status_col:

  Column that contains the status of each family member. Coded as 0 or
  FALSE (control) and 1 or TRUE (case).

- thr_col:

  Name of column to create for threshold. Defaults to "thr".

- lower_col:

  Name of column to create for lower threshold. Defaults to "lower".

- upper_col:

  Name of column to create for upper threshold. Defaults to "upper".

- lower_equal_upper:

  Should the upper and lower threshold be the same for cases? Can be
  used if CIPs are detailed, e.g. stratified by birth year and sex.

- personal_thr:

  Should thresholds be based on stratified CIPs or population
  prevalence?

- interpolation:

  Type of interpolation, defaults to NULL.

- bst.params:

  List of parameters to pass on to xgboost. See xgboost documentation
  for details.

- min_CIP_value:

  Minimum cip value to allow. Too low values may lead to numerical
  instabilities.

- xgboost_itr:

  Number of iterations to run xgboost for.

## Value

Tibble with (personlised) thresholds for each family member (lower &
upper), the calculated cumulative incidence proportion for each
individual (K_i), and population prevalence within an individuals CIP
strata (K_pop; max value in stratum). The threshold and other
potentially relevant information can be added to the family graphs with
`familywise_attach_attributes`.

## Examples

``` r
tbl = data.frame(
fid = c(1, 1, 1, 1),
pid = c(1, 2, 3, 4),
role = c("o", "m", "f", "pgf"),
sex = c(1, 0, 1, 1),
status = c(0, 0, 1, 1),
age = c(22, 42, 48, 78),
birth_year = 2023 - c(22, 42, 48, 78),
aoo = c(NA, NA, 43, 45))

cip = data.frame(
age = c(22, 42, 43, 45, 48, 78),
birth_year = c(2001, 1981, 1975, 1945, 1975, 1945),
sex = c(1, 0, 1, 1, 1, 1),
cip = c(0.1, 0.2, 0.3, 0.3, 0.3, 0.4))

prepare_thresholds(.tbl = tbl, CIP = cip, age_col = "age", interpolation = NA)
#>   fid pid role sex status age birth_year aoo cip       thr     lower     upper
#> 1   1   1    o   1      0  22       2001  NA 0.1 1.2815516      -Inf 1.2815516
#> 2   1   2    m   0      0  42       1981  NA 0.2 0.8416212      -Inf 0.8416212
#> 3   1   3    f   1      1  48       1975  43 0.3 0.5244005 0.5244005       Inf
#> 4   1   4  pgf   1      1  78       1945  45 0.4 0.2533471 0.2533471       Inf
```
