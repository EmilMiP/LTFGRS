# LTFGRS Workflow: Prepare Genetic Liability for Multi-trait GWAS

**This vignette is heavily inspired by and re-uses code from the “LTFGRS
Workflow: Prepare Genetic Liability for GWAS” vignette.**

***All data is simulated and is purely for demonstration purposes - This
includes the CIP.***

First, we load the required packages.

``` r
library(LTFGRS)
library(dplyr)
library(lubridate)
library(rmarkdown)
library(stringr)
```

## Simulate mock trio, phenotype, and CIP data

We will set some population parameters for the simulation. The
parameters are as follows:

``` r
set.seed(555)
Ntrait = 3 # number of traits to consider
h2_vec = rep(0.5, Ntrait) # heritability
corMat = matrix(0.3, nrow = Ntrait, ncol = Ntrait) # genetic correlation matrix
diag(corMat) = 1
K = .3 # population prevalence
all_outcomes = paste0("disorder_", c("1", "2", "3")) # disorder names
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

Since we are simulating multiple phenotypes, we need a CIP for each
phenotype. We can create a list of CIPs for each phenotype as follows:

``` r
# creating dummy CIP list for multiple traits
CIP_list = lapply(1:Ntrait, function(i) {
  CIP
})
```

When applying to real data, the order of the CIPs should match the order
of the phenotypes in the phenotype data. This is important for the
function
[`prepare_thresholds_multi()`](https://emilmip.github.io/LTFGRS/reference/prepare_thresholds_multi.md)
to correctly match the CIPs to the phenotypes.

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

### Phenotype data

We will simulate a liability based on the family structure defined above
to assign a case-control outcome to each individual. Then other
covariates such as sex and age are randomly assigned. To get the
case-control status, we first generate a (population) graph, calculate a
kinship matrix based on the heritability, genetic correlation, and
kinship coefficient, and finally, draw liabilities from a multivariate
normal with the calculated kinship matrix as covariance matrix. However,
since it is a multi-trait application, the covariance matrix derived
from the trio information is $N_{trait}$ times larger in all directions
compared to a single-trait application.

``` r

# creating a graph for the family
graph = prepare_graph(.tbl = family, icol = "id", mcol = "momcol", fcol = "dadcol")
# calculating the kinship matrix based on the graph
cov_mat_obj = graph_based_covariance_construction_multi(
  fid = "fam",
  pid = "pid",
  cur_proband_id = "pid",
  cur_family_graph = graph,
  h2_vec = h2_vec,
  genetic_corrmat = corMat,
  useMixture = FALSE,
  phen_names = c("1", "2", "3"),
  add_ind = TRUE)

# creating a phenotype for the family
liabs = MASS::mvrnorm(n = 1, mu = rep(0, nrow(cov_mat_obj$cov)), Sigma = cov_mat_obj$cov)
```

``` r
pheno = tibble(
  id = names(liabs),
  status = liabs > qnorm(K, lower.tail = F),
  # no consideration for generation etc in sex, fdato or birth_year:
  sex = rep(sample(0:1, size = length(liabs)/Ntrait, replace = TRUE), Ntrait),
  fdato = rep(dmy(paste0(sample(1:28, length(liabs)/Ntrait, replace = T), "/", sample(1:12, length(liabs)/Ntrait, replace = T), "/", sample(1940:2000, length(liabs)/Ntrait, replace = T))), Ntrait),
  birth_year = year(fdato),
  # age of onset only after fdato
  disorder = purrr::map2_chr(.x = status, .y = birth_year,
                         ~ if(.x) paste0(sample(1:28, 1), "/", sample(1:12, 1), "/", sample((.y + 1):2010, 1)) else NA),
  # end of follow up assigned here
  indiv_eof = dmy("31/12/2010")) %>% # blanket time stop, meant to simulate end of registers
  mutate(
    # converting to date format
    disorder = dmy(disorder),
    # eof either blanket time stop or event date
    indiv_eof = pmin(indiv_eof, disorder, na.rm = T),
     # calculating age at the end of follow up
    age_eof = as.numeric(difftime(indiv_eof, fdato, units = "days")) / 365.25) %>%
  filter(str_detect(id, "pid_g", negate = T))  # remove the genetic liability of the proband
paged_table(pheno)
```

An isolated look at one individual shows the additional information we
have on an individual compared to single-trait:

``` r
pheno %>% filter(str_detect(id, "pid"))
#> # A tibble: 3 × 8
#>   id    status   sex fdato      birth_year disorder   indiv_eof  age_eof
#>   <chr> <lgl>  <int> <date>          <dbl> <date>     <date>       <dbl>
#> 1 pid_1 TRUE       1 1960-10-03       1960 1980-09-22 1980-09-22    20.0
#> 2 pid_2 FALSE      1 1960-10-03       1960 NA         2010-12-31    50.2
#> 3 pid_3 TRUE       1 1960-10-03       1960 1972-06-24 1972-06-24    11.7
```

Next, we format the phenotype data into a format that can be derived
from most biobanks or electronic health records. The format is such that
each row is an individual and columns specify information on that
individual. Phenotype specific information follow the naming convention
of `<base_string>_<disorder_name>`. For example, `age_eof_disorder_1` is
the age at end of follow up for disorder 1.

``` r
# helper function - from: https://stackoverflow.com/questions/7963898/extracting-the-last-n-characters-from-a-string-in-r
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
pheno = pheno %>%
  mutate(trait = substrRight(id, 1),
         id = substr(id, 1, nchar(id) - 2)) %>%
  tidyr::pivot_wider(names_from = trait,
                     values_from = c(disorder, indiv_eof, status, age_eof),
                     names_glue = "{.value}_{trait}")
# sorting out names, such that they follow the required pattern with suffix ending in one of the disorder names:
colnames(pheno)[-(1:7)] = str_replace_all(colnames(pheno)[-(1:7)], c("_1" = "_disorder_1", "_2" = "_disorder_2", "_3" = "_disorder_3"))
```

The formatted phenotype data can be inspected below:

``` r
paged_table(pheno)
```

## Preparing for `estimate_liability()`

In a real world scenario, we will not have access to all of the
information used above. We will assume that the objects `CIP_list`,
`family`, and `pheno` are the only information available to the user.
These objects hold information that can often be extracted from
population registers or bio banks.

- `CIP_list`: The `CIP_list` object carry information about the
  prevalence all considered outcomes of interest in the population and
  therefore also on how each participant fits into the population
  distribution.
- `family`: The `family` object holds the trio information,
  i.e. information about the family structure and how each individual is
  related to each other. In a real world scenario this object may
  contains millions of unique individuals.
- `pheno`: The `pheno` object holds phenotypic information on each
  individual present in the trio information and on each outcome being
  considered.

Next, we calculate the personalised thresholds for each individual and
each phenotype. The function
[`prepare_thresholds_multi()`](https://emilmip.github.io/LTFGRS/reference/prepare_thresholds_multi.md)
will return a tibble with the calculated thresholds for each individual
and each phenotype. It is a wrapper function that applies the function
[`prepare_thresholds()`](https://emilmip.github.io/LTFGRS/reference/prepare_thresholds.md)
for each phenotype and combines the results into a single tibble. Due to
data privacy, it is possible to encounter CIPs values that are only
provided at set values, e.g. a CIP value for each whole year by birth
year and sex, such as what is shown in the `CIP` object. However, the
observed ages (or age of diagnosis) are typically not integer values.
This means we may need to approximate the CIP values between the
provided values. We offer an XGboost based approach to interpolate the
CIPs between the provided values.

``` r
# thresholds: personalised, fixed upper and lower threshold
multi_thrs_fixed = prepare_thresholds_multi(
  .tbl = pheno, CIP_list = CIP_list, phen_names = all_outcomes, personal_thr = TRUE, lower_equal_upper = TRUE
)
```

The resulting `multi_thrs_fixed` object holds the following information
on each provided disorder: the lower and upper thresholds for each
individual, the individual prevalence, `K_i`, and the population
prevalence, `K_pop`, based on the provided CIP data. The thresholds are
personalised for each individual based on their age, sex, and birth
year. Notably, if there is a high degree of confidence in the accuracy
of the population representative CIPs stratified by birth year and sex
(and if possible, other defining features), then the upper and lower
thresholds may be fixed at the same value. This is done by setting
`lower_equal_upper = TRUE`. The `multi_thrs_fixed` object can be
inspected below:

``` r
paged_table(multi_thrs_fixed)
```

### Population graph

With the `family` object, which holds the trio information, we can
construct a population graph. The population graph holds all familial
connects identified in the trio information and will form the basis of
how families are identified. In real-world applications, the population
graph may contain millions of individuals. Here, we construct the
population graph with the threshold and prevalence information attached
for each disorder at creation. This means each individual in the graph
will have `lower_disorder_1`, `lower_disorder_2`, `lower_disorder_3`,
and so on for the `upper`, `K_i`, and `K_pop` columns too.

We generate the (population) graph

``` r
graph = prepare_graph(.tbl = family, 
                      icol = "id",
                      mcol = "momcol", 
                      fcol = "dadcol",
                      node_attributes = select(multi_thrs_fixed, id,
                                               contains("lower"), 
                                               contains("upper"),
                                               contains("K_i"),
                                               contains("K_pop")))
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
                                  proband_vec = multi_thrs_fixed$id,
                                  fid = "fid",
                                  fam_graph_col = "fam_graph")
family_graphs %>% print(n = 4)
#> # A tibble: 31 × 2
#>   fid      fam_graph
#>   <chr>    <list>   
#> 1 mhs      <igraph> 
#> 2 mauntH   <igraph> 
#> 3 dad2     <igraph> 
#> 4 pacousin <igraph> 
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
both of their children, and as a proband with a graph centred on them.
`fam_graph_col` holds the family graphs and are in the format of igraph.
Operations on this level will not be required for the average user. An
igraph object is shown here for context:

``` r
family_graphs$fam_graph[[1]]
#> IGRAPH e02f16a DN-- 3 2 -- 
#> + attr: name (v/c), lower_disorder_1 (v/n), lower_disorder_2 (v/n),
#> | lower_disorder_3 (v/n), upper_disorder_1 (v/n), upper_disorder_2
#> | (v/n), upper_disorder_3 (v/n), K_i_disorder_1 (v/n), K_i_disorder_2
#> | (v/n), K_i_disorder_3 (v/n), K_pop_disorder_1 (v/n), K_pop_disorder_2
#> | (v/n), K_pop_disorder_3 (v/n)
#> + edges from e02f16a (vertex names):
#> [1] dad2->mhs mom ->mhs
```

``` r
PA_MT = estimate_liability(family_graphs = family_graphs,
                           h2 = h2_vec,
                           genetic_corrmat = corMat,
                           phen_names = all_outcomes,
                           full_corrmat = diag(Ntrait),
                           target_phenotype = "disorder_1")
#> The number of workers is 1

Gibbs_MT = estimate_liability(family_graphs = family_graphs,
                              h2 = h2_vec,
                              genetic_corrmat = corMat,
                              phen_names = all_outcomes,
                              full_corrmat = diag(Ntrait),
                              method = "Gibbs")
#> The number of workers is 1
```

***OBS:*** The PA estimation approach currently only supports returning
a genetic liability estimate for one of the provided outcomes at a time.
The disorder that is returned in specified with the argument
`target_phenotype`. The Gibbs sampling approach will return a genetic
liability estimate for all provided outcomes at the same time. This
means that if you want to get a genetic liability estimate for multiple
outcomes, you will need to run the PA estimation approach multiple
times, once for each outcome, while the Gibbs sampling approach can be
run once to get estimates for all outcomes at the same time.

The `PA_MT` estimate can be inspected below:

``` r
paged_table(PA_MT)
```

The `Gibbs_MT` estimate can be inspected below:

``` r
paged_table(Gibbs_MT)
```

### Mixture model genetic liability estimates

Following the same steps as above and with the same type of input data,
we we can also estimate multi-trait genetic liability estimates with the
mixture model approach:

``` r
# thresholds: personalised, fixed upper and lower threshold
multi_thrs_mixture = prepare_thresholds_multi(
  .tbl = pheno, 
  CIP_list = CIP_list, 
  phen_names = all_outcomes, 
  personal_thr = FALSE, 
  lower_equal_upper = FALSE
)

# creating a graph for the family with mixture model thresholds
graph_mixture = prepare_graph(.tbl = family, 
                              icol = "id",
                              mcol = "momcol",
                              fcol = "dadcol",
                              node_attributes = select(multi_thrs_mixture, id,
                                                       contains("lower"), 
                                                       contains("upper"),
                                                       contains("K_i"),
                                                       contains("K_pop")))
# Identify family members of degree n
family_graphs_mixture = get_family_graphs(pop_graph = graph_mixture,
                                          ndegree = 1,
                                          proband_vec = multi_thrs_mixture$id,
                                          fid = "fid",
                                          fam_graph_col = "fam_graph_mixture")


PA_MT_mix = estimate_liability(family_graphs = family_graphs_mixture,
                               family_graphs_col = "fam_graph_mixture",
                               h2 = h2_vec,
                               genetic_corrmat = corMat,
                               phen_names = all_outcomes,
                               full_corrmat = diag(Ntrait),
                               target_phenotype = "disorder_1", 
                               useMixture = TRUE)
#> The number of workers is 1
paged_table(PA_MT_mix)
```
