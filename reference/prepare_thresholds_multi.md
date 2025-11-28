# Prepare thresholds for multiple phenotypes

This function is a wrapper around `prepare_thresholds` to prepare
thresholds for multiple phenotypes at once.

## Usage

``` r
prepare_thresholds_multi(
  .tbl,
  CIP_list,
  phen_names,
  CIP_merge_columns = c("sex", "birth_year", "age"),
  CIP_cip_col = "cip",
  age_col = "age",
  age_eof_base = "age_eof",
  status_col_base = "status",
  lower_base = "lower",
  upper_base = "upper",
  thr_col = "thr",
  K_i_col = "K_i",
  K_pop_col = "K_pop",
  Kpop_list = "useMax",
  personal_thr = TRUE,
  lower_equal_upper = FALSE,
  interpolation = "xgboost",
  bst.params = list(max_depth = 10, base_score = 0, nthread = 4, min_child_weight = 10),
  min_CIP_value = 1e-05,
  xgboost_itr = 30
)
```

## Arguments

- .tbl:

  Tibble with family and personal id columns, as well as
  CIP_merge_columns and status for each phenotype.

- CIP_list:

  List of tibbles with population representative cumulative incidence
  proportions. Each tibble must contain columns from `CIP_merge_columns`
  and `cIP_cip_col`.

- phen_names:

  Vector of phenotype names. Used to identify status columns and to name
  output columns.

- CIP_merge_columns:

  The columns the CIPs are subset by, e.g. CIPs by birth_year, sex. and
  age_col.

- CIP_cip_col:

  Name of column with CIP values.

- age_col:

  Name of column with age at the end of follow-up or age at diagnosis
  for cases.

- age_eof_base:

  Base name of age at end of follow-up column. The actual column name is
  constructed by appending the phenotype name. Defaults to "age_eof".

- status_col_base:

  Base name of status column. The actual column name is constructed by
  appending the phenotype name. Defaults to "status".

- lower_base:

  Base name of lower threshold column. The actual column name is
  constructed by appending the phenotype name. Defaults to "lower".

- upper_base:

  Base name of upper threshold column. The actual column name is
  constructed by appending the phenotype name. Defaults to "upper".

- thr_col:

  Base name of threshold column. The actual column name is constructed
  by appending the phenotype name. Defaults to "thr".

- K_i_col:

  Base name of individual cumulative incidence proportion column. The
  actual column name is constructed by appending the phenotype name.
  Defaults to "K_i".

- K_pop_col:

  Base name of population prevalence column. The actual column name is
  constructed by appending the phenotype name. Defaults to "K_pop".

- Kpop_list:

  List of population prevalence tibbles or "useMax" for each phenotype.
  If a tibble is provided, it must contain columns from `.tbl` and a
  column named "K_pop" with population prevalence values. Defaults to
  "UseMax".

- personal_thr:

  Should thresholds be based on stratified CIPs or population
  prevalence?

- lower_equal_upper:

  Should the upper and lower threshold be the same for cases? Can be
  used if CIPs are detailed, e.g. stratified by birth year and sex.

- interpolation:

  Type of interpolation, defaults to "xgboost".

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
strata (K_pop; max value in stratum) for each phenotype. The threshold
and other potentially relevant information can be added to the family
graphs with `familywise_attach_attributes`.

## Examples

``` r
# TODO: create simple example
```
