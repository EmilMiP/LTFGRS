# (Simplified) Kendler's FGRS

Function to calculate the simplified version of Kendler's FGRS based on
family data.

## Usage

``` r
kendler_simplified(
  .tbl = NULL,
  family_graphs = NULL,
  family_graphs_col = "fam_graph",
  pid = "pid",
  fid = "fid",
  role = NULL,
  dadcol,
  momcol,
  env_cor_sib = 0,
  env_cor_f = 0,
  env_cor_m = 0
)
```

## Arguments

- .tbl:

  A matrix, list or data frame that can be converted into a tibble. Must
  have at least five columns that hold the family identifier, the
  personal identifier, the role and the lower and upper thresholds. Note
  that the role must be one of the following abbreviations

  - `g` (Genetic component of full liability)

  - `o` (Full liability)

  - `m` (Mother)

  - `f` (Father)

  - `c[0-9]*.[0-9]*` (Children)

  - `mgm` (Maternal grandmother)

  - `mgf` (Maternal grandfather)

  - `pgm` (Paternal grandmother)

  - `pgf` (Paternal grandfather)

  - `s[0-9]*` (Full siblings)

  - `mhs[0-9]*` (Half-siblings - maternal side)

  - `phs[0-9]*` (Half-siblings - paternal side)

  - `mau[0-9]*` (Aunts/Uncles - maternal side)

  - `pau[0-9]*` (Aunts/Uncles - paternal side).

  Defaults to `NULL`. If `.tbl` is provided, `family_graphs` must be
  `NULL`.

- family_graphs:

  A tibble with columns pid and family_graph_col, dadcol, and momcol.
  See prepare_graph for construction of the graphs. The family graphs
  Defaults to NULL.

- family_graphs_col:

  Name of column with family graphs in family_graphs. Defaults to
  "fam_graph".

- pid:

  A string holding the name of the column in `.tbl` (or `family` and
  `threshs`) that hold the personal identifier(s). Defaults to "PID".

- fid:

  A string holding the name of the column in `.tbl` or `family` that
  holds the family identifier. Defaults to "fid".

- role:

  A string holding the name of the column in `.tbl` that holds the role.
  Each role must be chosen from the following list of abbreviations

  - `g` (Genetic component of full liability)

  - `o` (Full liability)

  - `m` (Mother)

  - `f` (Father)

  - `c[0-9]*.[0-9]*` (Children)

  - `mgm` (Maternal grandmother)

  - `mgf` (Maternal grandfather)

  - `pgm` (Paternal grandmother)

  - `pgf` (Paternal grandfather)

  - `s[0-9]*` (Full siblings)

  - `mhs[0-9]*` (Half-siblings - maternal side)

  - `phs[0-9]*` (Half-siblings - paternal side)

  - `mau[0-9]*` (Aunts/Uncles - maternal side)

  - `pau[0-9]*` (Aunts/Uncles - paternal side).

  Defaults to "role".

- dadcol:

  column name of father in family_graphs or .tbl.

- momcol:

  column name of mother in family_graphs or .tbl.

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

A tibble with summary values used to calculate the simplified kendler
FGRS and the FGRS itself.

## Details

The coding of the cohabitation effects differ slightly from the one
suggested by Kendler et al. Here, it is coded as env_cor\_\\ =
env_eff\_\\ / (gen_eff\_\\ + env_eff\_\\), while the original
implementation suggestes coding it as gen_eff\_\\ / (gen_eff\_\\ +
env_eff\_\\), i.e. the two are related as 1 - env_cor\_\\.

## Examples

``` r
# See Vignettes.
```
