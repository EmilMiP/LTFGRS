# Familywise Censoring for Multiple Outcomes

This function applies familywise censoring across multiple outcomes in a
dataset. If a target outcome is specified, all outcomes are censored
based on the end of follow-up time of the target outcome, then
familywise censoring is performed for each outcome individually.

## Usage

``` r
familywise_censoring_multi(
  family_graphs,
  tbl,
  start,
  end_base,
  phen_names,
  target_outcome = NULL,
  status_col_base = "status",
  aod_col_base = "aod",
  age_eof_col_base = "age_eof",
  fam_graph_col = "fam_graph",
  fid = "fid",
  pid = "pid",
  simplify = TRUE,
  merge_by = pid
)
```

## Arguments

- family_graphs:

  A tibble containing family graphs for each family.

- tbl:

  A tibble containing individual-level data with multiple outcomes.

- start:

  The column name representing the start of follow-up time.

- end_base:

  The base name for the end of follow-up time columns (e.g., "end" for
  columns like "end_outcome1", "end_outcome2").

- phen_names:

  A vector of phenotype names corresponding to the different outcomes.

- target_outcome:

  An optional string specifying the target outcome for censoring. If
  provided, all outcomes will be censored based on this outcome's end of
  follow-up time.

- status_col_base:

  The base name for the status columns (default is "status").

- aod_col_base:

  The base name for the age of diagnosis columns (default is "aod").

- age_eof_col_base:

  The base name for the age at end of follow-up columns (default is
  "age_eof").

- fam_graph_col:

  The column name in \`family_graphs\` that contains the family graph
  (default is "fam_graph").

- fid:

  The column name representing family ID (default is "fid").

- pid:

  The column name representing individual ID (default is "pid").

- simplify:

  A logical indicating whether to simplify the output by removing
  columns not specific to an outcome (default is TRUE).

- merge_by:

  The column name(s) to merge results by (default is \`pid\`).

## Value

A tibble containing the familywise censored results for all outcomes.

## Examples

``` r
# TODO: Add examples
```
