# LTFGRS Workflow: Prepare Genetic Liability for GWAS

**This vignette is heavily inspired by and re-uses code from the “LTFGRS
Workflow: Prepare Genetic Liability for Prediction” vignette.**. ***All
data is simulated and is purely for demonstration purposes - This
includes the CIP.***

This is going to be a vignette similar to “LTFGRS Workflow: Prepare
Genetic Liability for Prediction”, but this one focuses on preparing a
genetic liability score for use as a refined outcome for a genome-wide
association study (GWAS). If your intended use case is prediction, the
please follow the steps outlined in that vignette instead, since some
additional steps may be required.

In this vignette, we will

- simulate mock data
- calculate genetic liability in increasing degree of complexity,
  non-personalised thresholds, personalised based on stratified CIPs,
  mixture model. For each:
  - derive thresholds (and K_i and K_pop if needed)
  - attach information to pop_graph when constructing it
  - get family graphs of n’th degree

First, we load the required packages.

``` r
library(LTFGRS)
library(dplyr)
library(lubridate)
library(rmarkdown)
```

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
maternal uncle, “hsmcousin” for half-sibiling maternal cousin, etc. The
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

### Phenotype data

We will simulate a liability based on the family structure defined above
to assign a case-control outcome to each individual. Then other
covariates such as sex and age are randomly assigned. To get the
case-control status, we first generate a (population) graph, calculate a
kinship matrix based on the heritability and kinship coefficient, and
finally, draw liabilities from a multivariate normal with the calculated
kinship matrix as covariance matrix.

``` r
# creating a graph for the family
graph = prepare_graph(.tbl = family, icol = "id", mcol = "momcol", fcol = "dadcol")
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
  # no consideration for generation etc in sex, fdato or birth_year:
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
typically be derived from most register or biobank phenotypes where the
age of diagnosis is available. Columns of interest are:

- `fdato`: the birth date of the individual
- `birth_year`: the birth year of the individual
- `adhd`: the outcome of interest,
- `indiv_eof`: the personalised end of follow up for a given individual.
  it may be different for each individual due to any number of censoring
  or competing events.

## Preparing for `estimate_liability()`

In a real world scenario, we will not have access to all of the
information used above. We will assume that the objects `CIP`, `family`,
and `pheno` are the only information available to the user. These
objects hold information that can often be extracted from population
registers or bio banks.

- `CIP`: The `CIP` object carry information about the prevalence of the
  outcome of interest in the population and therefore also on how each
  participant fits into the population distribution.
- `family`: The `family` object holds the trio information,
  i.e. information about the family structure and how each individual is
  related to each other. In a real world scenario this object may
  contains millions of unique individuals.
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

The resulting `thresholds` object holds the lower and upper thresholds
for each individual, as well as the individual prevalence, `K_i`, and
the population prevalence, `K_pop`, based on the provided CIP data. The
thresholds are personalised for each individual based on their age, sex,
and birth year. Notably, if there is a high degree of confidence in the
accuracy of the population representative CIPs stratified by birth year
and sex (and if possible, other defining features), then the upper and
lower thresholds may be fixed at the same value. This is done by setting
`lower_equal_upper = TRUE`. The `thresholds` object can be inspected
below:

``` r
paged_table(thresholds)
```

### Population graph

With the `family` object, which holds the trio information, we can
construct a population graph. The population graph holds all familial
connects identified in the trio information and will form the basis of
how families are identified. In real-world applications, the population
graph may contain millions of individuals. Here, we construct the
population graph with all required information already attaches, namely
`lower`, `upper`, `K_i`, and `K_pop`.

``` r
graph = prepare_graph(.tbl = family, 
                      icol = "id",
                      mcol = "momcol", 
                      fcol = "dadcol",
                      node_attributes = select(thresholds, id, lower, upper, K_i, K_pop))
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
family graph).

``` r
# Identify family members of degree n
family_graphs = get_family_graphs(pop_graph = graph,
                                  ndegree = 1,
                                  proband_vec = pheno$id,
                                  fid = "fid",
                                  fam_graph_col = "fam_graph")
family_graphs %>% print(n = 4)
#> # A tibble: 31 × 2
#>   fid   fam_graph
#>   <chr> <list>   
#> 1 dad   <igraph> 
#> 2 mom   <igraph> 
#> 3 dad2  <igraph> 
#> 4 mom2  <igraph> 
#> # ℹ 27 more rows
```

The function
[`get_family_graphs()`](https://emilmip.github.io/LTFGRS/reference/get_family_graphs.md)
will return a formatted tibble. The output will have two columns
specified with the arguments, `fid` and `fam_graph_col`. `fid` is the
ids of the provided probands, who are also the individuals the
neighbourhood (family) graphs are centred on. Note: In the example
above, we have only one family graph for a given proband, however, an
individual may still appear in several family graphs as a relative.
E.g., a parent with two children may appearing in the family graph of
both of their children. `fam_graph_col` holds the family graphs and are
in the format of igraph. Operations on this level will not be required
for the average user. An igraph object is shown here for context:

``` r
family_graphs$fam_graph[[1]]
#> IGRAPH 98c42df DN-- 8 17 -- 
#> + attr: name (v/c), lower (v/n), upper (v/n), K_i (v/n), K_pop (v/n)
#> + edges from 98c42df (vertex names):
#>  [1] pgm   ->paunt  pgm   ->puncle pgm   ->dad    paunt ->puncle paunt ->dad   
#>  [6] sib   ->pid    puncle->paunt  puncle->dad    pid   ->sib    pgf   ->paunt 
#> [11] pgf   ->puncle pgf   ->dad    dad   ->paunt  dad   ->sib    dad   ->puncle
#> [16] dad   ->pid    dad   ->phs
```

## Estimating genetic liabilities with `estimate_liability()`

The function
[`estimate_liability()`](https://emilmip.github.io/LTFGRS/reference/estimate_liability.md)
is used to estimate the genetic liability. The function accepts two
types of input, here we will only focus on the graph-based input
generated above. The graph-based input offer the best flexibility and
scalability. The function has two arguments that are worth pointing out.

The first is `method`, which specifies the method used to estimate the
genetic liability. Currently, two methods are supported. The first is a
Gibbs sampler that samples from a truncated multivariate normal
distribution, `method = "Gibbs"`. The second is an iterative
Pearson-Aitken approach, `method = "PA"`. Generally speaking, the
Pearson-Aitken approach is faster.

The second argument is `useMixture`, which specifies whether to use the
mixture model or not. The mixture model is currently only supported with
`method = "PA"`. The mixture model considers the genetic liability of
controls as a mixture of the truncated normal for cases and controls,
rather than just the distribution of controls. This accounts for the
possibility that some controls are undiagnosed cases and accounts for it
in the genetic liability estimate.

### FGRS with personalised threshold and PA

``` r
ltfgrs_pa = estimate_liability(family_graphs = family_graphs,
                               h2 = h2, 
                               fid = "fid",
                               pid = "pid",
                               family_graphs_col = "fam_graph",
                               method = "PA", # <- METHOD
                               useMixture = F)
#> The number of workers is 1
paged_table(ltfgrs_pa)
```

When using `method = "PA"`, an iterative conditioning is performed,
which means the resulting estimate and uncertainty of the estimate is
the expected mean value and variance of the last iteration, which is the
proband’s genetic liability. This is highlighted by the use of `var` in
the output.

### FGRS with personalised threshold and Gibbs

``` r
ltfgrs_gibbs = estimate_liability(family_graphs = family_graphs, 
                                  h2 = h2,
                                  fid = "fid",
                                  pid = "pid",
                                  family_graphs_col = "fam_graph",
                                  method = "Gibbs", # <- METHOD
                                  useMixture = F)
#> The number of workers is 1
paged_table(ltfgrs_gibbs)
```

### FGRS with mixture model

Finally, it is also possible to use the PA estimation method to estimate
the genetic liability using the mixture model. This is done by setting
`useMixture = TRUE`. It requires that each individual has the upper and
lower thresholds, as well as `K_i` and `K_pop` attached as node
attributes in the family graphs.

``` r
ltfgrs_mixture = estimate_liability(family_graphs = family_graphs,
                                    h2 = h2, 
                                    fid = "fid",
                                    pid = "pid",
                                    family_graphs_col = "fam_graph",
                                    method = "PA", # <- METHOD
                                    useMixture = TRUE) # <- useMixture model 
#> The number of workers is 1
paged_table(ltfgrs_mixture)
```

### Parallelisation

The function
[`estimate_liability()`](https://emilmip.github.io/LTFGRS/reference/estimate_liability.md)
is able to use the `future` package to parallelise the estimation of the
genetic liability. This is done by setting a suitable `plan` with the
future backend. A plan suitable for most needs is
`plan(multisession, workers = NCORES)`, which means that the function
will run in parallel on the local PC utilising `NCORES`-cores. Other
parallelisation options exist, but they are all handled by the future
suit of packages.
