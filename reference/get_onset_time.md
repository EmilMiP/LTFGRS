# Calculate age of diagnosis, age at end of follow up, and status

Calculate age of diagnosis, age at end of follow up, and status

## Usage

``` r
get_onset_time(
  tbl,
  start,
  end,
  event,
  status_col = "status",
  aod_col = "aod",
  age_eof_col = "age"
)
```

## Arguments

- tbl:

  tibble with start, end, and event as columns

- start:

  start of follow up, typically birth date, must be a date column

- end:

  end of follow up, must be a date column

- event:

  event of interest, typically date of diagnosis, must be a date column

- status_col:

  column name of status column to be created. Defaults to "status".

- aod_col:

  column name of age of diagnosis column to be created. Defaults to
  "aod".

- age_eof_col:

  column name of age at end of follow-up column to be created. Defaults
  to "age_eof".

## Value

tibble with added status, age of diagnosis, and age at end of follow-up

## Examples

``` r
# See vignettes.
```
