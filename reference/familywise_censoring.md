# Censor Family Onsets for Multiple Families

This fucntion is a wrapper around `censor_family_onsets`. This functions
accepts a tibble with family graphs from `get_family_graphs`. It censors
the onset times for each individual in the family graph based on the
proband's end of follow-up. Returns a formatted output.

## Usage

``` r
familywise_censoring(
  family_graphs,
  tbl,
  start,
  end,
  event,
  status_col = "status",
  aod_col = "aod",
  age_eof_col = "age",
  fam_graph_col = "fam_graph",
  fid = "fid",
  pid = "pid",
  merge_by = pid
)
```

## Arguments

- family_graphs:

  Tibble with fid and family graphs columns.

- tbl:

  Tibble with information on each considered individual.

- start:

  Column name of start of follow up, typically date of birth.

- end:

  Column name of the personalised end of follow up.

- event:

  Column name of the event.

- status_col:

  Column name of the status (to be created). Defaults to "status".

- aod_col:

  Column name of the age of diagnosis (to be created). Defaults to
  "aod".

- age_eof_col:

  Column name of the age at the end of follow up (to be created).
  Defaults to "age_eof".

- fam_graph_col:

  Column name of family graphs in the 'family_graphs' object. Defaults
  to "fam_graph".

- fid:

  Family id, typically the name of the proband that a family graph is
  centred on. Defaults to "fid".

- pid:

  Personal identifier for each individual. Allows for multiple instances
  of the same individual across families. Defaults to "pid".

- merge_by:

  Column names to merge by. If different names are used for family
  graphs and tbl, a named vector can be specified: setNames(c("id"),
  c("pid")). Note id is the column name in tbl and pid is the column
  name in family_graphs. The column names used should reference the
  personal identifier.

## Value

A tibble with family ids and updated status, age of diagnosis, and age
at end of follow-up for each individual in the family based on the
proband's end of follow-up.

## Examples

``` r
# See Vignettes.
```
