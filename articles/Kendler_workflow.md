# Kendler workflow

``` r
library(LTFGRS)
library(dplyr)
library(lubridate)
library(rmarkdown)
```

***All data is simulated and is purely for demonstration purposes - This
includes the CIP.***

This vignette is intended to give users of the package an overview of
how to estimate a simplified version of the kendler’s FGRS. See
<https://doi.org/10.1001/jamapsychiatry.2021.0336> for more details on
the full method. The implementation provided here is heavily inspired by
the PAFGRS package (<https://github.com/BioPsyk/PAFGRS>).

The purpose of this vignette is to illustrate how to estimate the
simplified Kendler FGRS in practice. The simulated data is intended to
represent real-world data formats, while not being real data. As such,
generating the mock data may be more complex than what is otherwise
strictly necessary.

The vignette will cover the following steps:

- Simulate mock trio, phenotype, and CIP data
- Perform steps needed for FGRS calculations:
  - Assign thresholds to individuals based on CIP and phenotype data
  - Create a population graph from the trio data with attached node
    attributes
  - Automatic identification of n-degree relatives
- Estimate the (simplified) Kendler FGRS
  - Using both family graph and long format input

## Simulate mock trio, phenotype, and CIP data

We will set some population parameters for the simulation. The
parameters are as follows:

``` r
# population parameters and seed
set.seed(555)
h2 = .5 # heritability
K = .3 # population prevalence
```

### Cumulative incidence proportions (CIP)

One of the key required input variables of LTFGRS is the population
representative stratified cumulative incidence proportions (CIP) data.
LTFGRS is able to utilise the population representative stratified CIPs
to personalise thresholds for the liability-based predictors. The CIPs
are typically obtained from large population registers or other sources
that allow for population representative estimates. Here, we simulate a
format similar to how stratified CIPs may be stored. We assume the CIPs
have been stratified by sex and birth year. The population
representative stratified CIPs has the interpretation of being the
proportion of individuals born in a given year and sex that has been
diagnosed with the outcome of interest by age $x$.

``` r
# assuming we have been provided a CIP object of the following style:
CIP = expand.grid(list(age = 1:100,
                       birth_year = 1900:2024,
                       sex = 0:1)) %>%
  group_by(sex, birth_year) %>%
  mutate(cip = (1:n() - 1)/n() * K) %>%
  ungroup() %>% 
  print(n = 10)
#> # A tibble: 25,000 × 4
#>      age birth_year   sex   cip
#>    <int>      <int> <int> <dbl>
#>  1     1       1900     0 0    
#>  2     2       1900     0 0.003
#>  3     3       1900     0 0.006
#>  4     4       1900     0 0.009
#>  5     5       1900     0 0.012
#>  6     6       1900     0 0.015
#>  7     7       1900     0 0.018
#>  8     8       1900     0 0.021
#>  9     9       1900     0 0.024
#> 10    10       1900     0 0.027
#> # ℹ 24,990 more rows
```

### Trio information

The trio information presented here is a manually constructed to
resemble a typical way the trio data may be stored. The names are chosen
such that they resemble the relationship to the proband. This means
there are simple names such as “dad”, “mom”, or “sib”. There are also
more complex names such as “pgf” for paternal grand father, “muncle” for
maternal uncle, “hsmcousin” for half-sibling maternal cousin, etc. The
suffixes “H” and “W” mean husband and wife, respectively.

``` r
# hand curated trio information, taken from LTFHPlus vignette:
# https://emilmip.github.io/LTFHPlus/articles/FromTrioToFamilies.html
family = tribble(
  ~id, ~momcol, ~dadcol,
  "pid", "mom", "dad",
  "sib", "mom", "dad",
  "mhs", "mom", "dad2",
  "phs", "mom2", "dad",
  "mom", "mgm", "mgf",
  "dad", "pgm", "pgf",
  "dad2", "pgm2", "pgf2",
  "paunt", "pgm", "pgf",
  "pacousin", "paunt", "pauntH",
  "hspaunt", "pgm", "newpgf",
  "hspacousin", "hspaunt", "hspauntH",
  "puncle", "pgm", "pgf",
  "pucousin", "puncleW", "puncle",
  "maunt", "mgm", "mgf",
  "macousin", "maunt", "mauntH",
  "hsmuncle", "newmgm", "mgf",
  "hsmucousin", "hsmuncleW", "hsmuncle"
)
```

In addition, we ensure that every individual has their own unique row in
the trio data:

``` r
# in real trio data, we would have a row for every individual in the population, and NA if parents are unknown:
to_add_momcol = tibble(id = unique(family$momcol[!(family$momcol %in% family$id)]),
                       momcol = NA,
                       dadcol = NA)
to_add_dadcol = tibble(id = unique(family$dadcol[!(family$dadcol %in% family$id)]),
                       momcol = NA,
                       dadcol = NA)
trio = bind_rows(family, to_add_momcol, to_add_dadcol)
paged_table(trio)
```

### Phenotype data

We will simulate a phenotypes based on the liability thershold model. We
simulate liability based on the family structure defined above to assign
a case-control outcome to each individual. Then other covariates such as
sex and age are randomly assigned. To get the case-control status, we
first generate a (population) graph, calculate a kinship matrix based on
the heritability and kinship coefficient, and finally, draw liabilities
from a multivariate normal with the calculated kinship matrix as
covariance matrix. This is to ensure the liabilities are correlated
according to the family structure.

``` r
# creating a graph for the family
graph = prepare_graph(.tbl = trio, icol = "id", mcol = "momcol", fcol = "dadcol")
# calculating the kinship matrix based on the graph
cov_mat = get_covmat(fam_graph = graph, h2 = h2, index_id = "pid")
# creating a phenotype for the family
liabs = MASS::mvrnorm(n = 1, mu = rep(0, nrow(cov_mat)), Sigma = cov_mat)
```

Next, we will create the mock phenotype data:

``` r
# these values are simulated only for illustrative purposes and not to make sense(!)
pheno = tibble(
  id = names(liabs),
  status = liabs > qnorm(K, lower.tail = F),
  # no consideration for generation etc in fdato or birth_year:
  fdato = dmy(paste0(sample(1:28, length(liabs), replace = TRUE), "/", sample(1:12, length(liabs), replace = T), "/", sample(1940:2000, length(liabs), replace = TRUE))),
  birth_year = year(fdato),
  # age of onset only after fdato
  adhd = purrr::map2_chr(.x = status, .y = birth_year,
      ~ if(.x) paste0(sample(1:28, 1), "/", sample(1:12, 1), "/", sample((.y + 1):2010, 1)) else NA),
  # end of follow up assigned here
  indiv_eof = dmy("31/12/2010")) %>% # blanket time stop, meant to simulate end of registers
  mutate(
    # Assigning sex to each individual
    sex = case_when(
    id %in% family$momcol ~ 1,
    id %in% family$dadcol ~ 0,
    TRUE ~ sample(0:1, n(), replace = TRUE)),
    # converting to date format
    adhd = dmy(adhd),
    # eof either blanket time stop or event date
    indiv_eof = pmin(indiv_eof, adhd, na.rm = TRUE),
    # calculating age at the end of follow up
    age = as.numeric(difftime(indiv_eof, fdato, units = "days")) / 365.25) %>% 
  filter(id != "pid_g") # remove the genetic liability of the proband
paged_table(pheno)
```

The mock phenotype data is intended to resemble a format that can
typically be derived from most register or bio bank phenotypes where the
age of diagnosis is available. Columns of interest are:

- `fdato`: the birth date of the individual
- `birth_year`: the birth year of the individual
- `adhd`: the outcome of interest,
- `indiv_eof`: the personalised end of follow up for a given individual.
  It may be different for each individual due to any number of censoring
  or competing events.

## Preparing for `kendler_simplified()`

In a real world scenario, we will ***not*** have access to all of the
information used above. We will assume that the objects `CIP`, `trio`,
and `pheno` are the only information available to the user. These
objects hold information that can often be extracted from population
registers or bio banks.

- `CIP`: The `CIP` object carry information about the prevalence of the
  outcome of interest in the population and therefore also on how each
  participant fits into the population distribution.
- `trio`: The `trio` object holds the trio information, i.e. information
  about the family structure and how each individual is related to each
  other. In a real world scenario this object may contains millions of
  unique individuals.
- `pheno`: The `pheno` object holds phenotypic information on each
  individual present in the trio information.

### Preparing thresholds and cumulative incidence proportions

Due to data privacy, it is possible to encounter CIPs values that are
only provided at set values, e.g. a CIP value for each whole year by
birth year and sex, such as what is shown in the `CIP` object. However,
the observed ages (or age of diagnosis) are typically not integer
values. This means we may need to approximate the CIP values between the
provided values. We offer an XGboost based approach to interpolate the
CIPs between the provided values.

``` r
thresholds = prepare_thresholds(
  .tbl = pheno,
  CIP = CIP,
  age_col = "age" ,
  status_col = "status",
  lower_equal_upper = FALSE,
  personal_thr = TRUE, 
  interpolation = "xgboost"
)
```

The function
[`prepare_thresholds()`](https://emilmip.github.io/LTFGRS/reference/prepare_thresholds.md)
has more functionality than is strictly needed for the
[`kendler_simplified()`](https://emilmip.github.io/LTFGRS/reference/kendler_simplified.md).
See documentation for details. It returns a tibble with the following
additional columns:

- `K_i`: The CIP value for the individual. `K_i` is predicted if
  interpolation is used.
- `K_pop`: The population prevalence. Presently, it is calculated as the
  maximum CIP value within the CIP stratum an individual belongs to,
  e.g. for a male born in $2000$, `K_pop` is the maximum CIP value
  observed among males born in the year $2000$. Alternatively, acquired
  user-specified values through the `Kpop` argument.
- `thr`: The liability threshold used to determine case-control status.
  `thr` is used to determine the upper and lower thresholds of an
  individual.
- `lower`: lower threshold of an individual.
- `upper`: upper threshold of an individual.

For the simplified Kendler FGRS, only `lower`, `upper`, and `K_i` are
needed.

### Population graph

With the `trio` object, we can construct a population graph. The
population graph holds all familial connects identified in the trio
information and will form the basis of how families are identified. In
real-world applications, the population graph may contain millions of
individuals. To illustrate all the required steps, we will recalculate
the population graph, but add `sex`, `lower`, `upper`, and `K_i` as a
node attribute.

``` r
graph = prepare_graph(.tbl = trio, 
                      icol = "id",
                      mcol = "momcol", 
                      fcol = "dadcol",
                      node_attributes = select(thresholds, id, sex, K_i, lower, upper))
```

### Automatic identification of n-degree relatives

When we want to calculate a family genetic risk score, we need to create
a pedigree based on the proband and relations should be relative to the
proband. We are interested in identifying all family members up to some
degree of relatedness, $n$, without having to manually find all of these
family members. Manually identifying family members up to degree $4$ is
both time consuming and error prone. We have implemented an automatic
detection of family members that utilise a graph based on all
individuals in the trio information (ideally population registers) and
neighbourhood graphs. In short, we create a pedigree (directed graph)
with every individual in the trio data and copy sections around a
proband with all individuals that are $n$ steps away from the proband in
the graph (This is a neighbourhood graphs of degree n, here called a
family graph). Below, we only identify all second degree family members:

``` r
# Identify family members of degree n
family_graphs = get_family_graphs(pop_graph = graph,
                                  ndegree = 2,
                                  proband_vec = pheno$id,
                                  fid = "fid",
                                  fam_graph_col = "fam_graph")
```

For Kendler’s FGRS, we need to know who the parents are for each proband
to account for cohabitation. We can attach this information to the
`family_graphs` object from the `trio` object:

``` r
family_graphs = left_join(family_graphs, trio, by = c("fid" = "id"))
family_graphs %>% print(n = 4)
#> # A tibble: 31 × 4
#>   fid   fam_graph momcol dadcol
#>   <chr> <list>    <chr>  <chr> 
#> 1 dad   <igraph>  pgm    pgf   
#> 2 mom   <igraph>  mgm    mgf   
#> 3 dad2  <igraph>  pgm2   pgf2  
#> 4 mom2  <igraph>  NA     NA    
#> # ℹ 27 more rows
```

This means the `family_graphs` object now contains 4 columns:

- `fid`: The family id. Here, the proband id the graph is centred on.
- `fam_graph`: the family graph, in iGraph format
- `momcol`: the mother id of the proband
- `dadcol`: the father id of the proband

The family graph contains information on each family member, namely the
sex, K_i, lower, and upper thresholds, stored as node attributes. A
typical print of an igraph can be seen below:

``` r
family_graphs$fam_graph[[1]]
#> IGRAPH 4fbc9c0 DN-- 13 23 -- 
#> + attr: name (v/c), sex (v/n), K_i (v/n), lower (v/n), upper (v/n)
#> + edges from 4fbc9c0 (vertex names):
#>  [1] mom   ->sib      mom   ->pid      pgm   ->hspaunt  pgm   ->paunt   
#>  [5] pgm   ->puncle   pgm   ->dad      mom2  ->phs      paunt ->puncle  
#>  [9] paunt ->pacousin paunt ->dad      sib   ->pid      puncle->pucousin
#> [13] puncle->paunt    puncle->dad      pid   ->sib      pgf   ->paunt   
#> [17] pgf   ->puncle   pgf   ->dad      dad   ->paunt    dad   ->sib     
#> [21] dad   ->puncle   dad   ->pid      dad   ->phs
```

## Estimating the simplified Kendler FGRS

We have now prepared the family graph input needed to estimate the
simplified Kendler FGRS. Below, we calculate the FGRS with a
cohabitation effect of $0.5$ for both parents and siblings.

``` r
kendler_fgrs = kendler_simplified(
  family_graphs = family_graphs, 
  family_graphs_col = "fam_graph",
  pid = "pid",
  fid = "fid",
  dadcol = "dadcol", 
  momcol = "momcol", 
  env_cor_sib = .5,
  env_cor_f = .5, 
  env_cor_m = .5)
paged_table(kendler_fgrs)
```

It is also possible to use a long format input instead of the family
graphs. Getting this type of input correct can be tedious, since it
requires manually identifying family members. A simple example is shown
below that is generated from the family_graphs object:

``` r
# extract graph
pid_graph = (family_graphs %>% filter(fid == "pid") %>% pull(fam_graph))[[1]] 
# extract attributes from graph and format for tbl input:
tbl_input = igraph::vertex.attributes(pid_graph) %>% 
  as_tibble() %>% 
  rename(pid = name) %>% 
  mutate(fid = "pid",
         # roles are relative to the proband of the family. Only using pid here:
         role = stringr::str_replace(pid, "uncle|aunt", "au"),
         role = case_when(
           role == "mom" ~ "m",
           role == "dad" ~ "f",
           role == "sib" ~ "s",
           role == "pid"  ~ "o",
           TRUE ~ role
         )) %>% 
  relocate(fid, pid, role) %>% 
  left_join(trio, by = c("fid" = "id"))
# more than 1 pau role, making sure we know they are different:
which_pau = which(tbl_input$role == "pau")
tbl_input$role[which_pau] = paste0(tbl_input$role[which_pau], 1:2)
paged_table(tbl_input)
```

The `tbl_input` object now holds all the required information to
estimate the simplified Kendler FGRS for the individual `pid` input data
in the long format.

``` r
kendler_simplified(.tbl = tbl_input,
                   role = "role",
                   dadcol = "dadcol",
                   momcol = "momcol")
#> # A tibble: 1 × 6
#>   pid       S  varZ    mZ sum_r   FGRS
#>   <chr> <dbl> <dbl> <dbl> <dbl>  <dbl>
#> 1 pid   0.149  1.11 0.458  3.75 0.0954
```

Note: The population normalisation done in
[`kendler_simplified()`](https://emilmip.github.io/LTFGRS/reference/kendler_simplified.md)
is performed with assumed variances of $1$ to match the underlying
normal distribution instead of the observed variances, because only one
individual is used to estimate the FGRS. If two or more individuals are
used, the observed variance is used for normalisation instead.
