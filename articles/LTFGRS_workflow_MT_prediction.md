# LTFGRS Workflow: Prepare Genetic Liability for Multi-Trait Prediction

**This vignette is heavily inspired by and re-uses code from the “LTFGRS
Workflow: Prepare Genetic Liability for Prediction” vignette.** ***All
data is simulated and is purely for demonstration purposes - This
includes the CIP.***

This vignette is intended to give users of the package an overview of
how to use the different functions of the package and how they are
intended to work together to calculate a multi-trait genetic liability.
The liability-based predictors can be used for prediction and as the
outcome in a GWAS. If you want a refined outcome to be used in GWAS,
please see the vignette titled “LTFGRS Workflow: Prepare Genetic
Liability for multi-trait GWAS”. This vignette focuses on preparing for
prediction. This purpose requires the user to consider additional steps,
such as censoring of future events and masking the outcome of the target
individual. The vignette will cover the following steps:

- Simulate mock trio, phenotype, and CIP data
- Convert mock trio and phenotype data into a format suitable for
  estimate_liability()
  - Automatic identification of n-degree relatives
  - Censoring of future events on a family basis
  - Assign thresholds on censored families
  - Assign thresholds to family graphs
- Estimate the liability-based predictors using estimate_liability()
  - Both Gibbs and PA based estimators will be covered

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
```

    ## # A tibble: 25,000 × 4
    ##      age birth_year   sex   cip
    ##    <int>      <int> <int> <dbl>
    ##  1     1       1900     0 0    
    ##  2     2       1900     0 0.003
    ##  3     3       1900     0 0.006
    ##  4     4       1900     0 0.009
    ##  5     5       1900     0 0.012
    ##  6     6       1900     0 0.015
    ##  7     7       1900     0 0.018
    ##  8     8       1900     0 0.021
    ##  9     9       1900     0 0.024
    ## 10    10       1900     0 0.027
    ## # ℹ 24,990 more rows

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
    indiv_eof = pmin(indiv_eof, disorder, na.rm = T)) %>% 
  filter(str_detect(id, "pid_g", negate = T)) %>%   # remove the genetic liability of the proband
  # status is not required here. Only used it to generate the age of diagnosis.
  select(-contains("status"))
paged_table(pheno)
```

An isolated look at one individual shows the additional information we
have on an individual compared to single-trait:

``` r
pheno %>% filter(str_detect(id, "pid"))
```

    ## # A tibble: 3 × 6
    ##   id      sex fdato      birth_year disorder   indiv_eof 
    ##   <chr> <int> <date>          <dbl> <date>     <date>    
    ## 1 pid_1     1 1960-10-03       1960 1980-09-22 1980-09-22
    ## 2 pid_2     1 1960-10-03       1960 NA         2010-12-31
    ## 3 pid_3     1 1960-10-03       1960 1972-06-24 1972-06-24

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
                     values_from = c(disorder, indiv_eof),
                     names_glue = "{.value}_{trait}")
# sorting out names, such that they follow the required pattern with suffix ending in one of the disorder names:
colnames(pheno)[-(1:7)] = str_replace_all(colnames(pheno)[-(1:7)], 
                                          c("_1" = "_disorder_1", 
                                            "_2" = "_disorder_2", 
                                            "_3" = "_disorder_3"))
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

### Population graph

With the `family` object, which holds the trio information, we can
construct a population graph. The population graph holds all familial
connects identified in the trio information and will form the basis of
how families are identified. In real-world applications, the population
graph may contain millions of individuals. We generate the (population)
graph

``` r
graph = prepare_graph(.tbl = family, 
                      icol = "id",
                      mcol = "momcol", 
                      fcol = "dadcol")
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
```

    ## # A tibble: 31 × 2
    ##   fid   fam_graph
    ##   <chr> <list>   
    ## 1 pid   <igraph> 
    ## 2 dad   <igraph> 
    ## 3 mom   <igraph> 
    ## 4 dad2  <igraph> 
    ## # ℹ 27 more rows

### Censoring of future events on a family basis

The purpose of the genetic liability estimated here is prediction. In
epidemiology (and many other fields) there is an emphasis on ensuring
future events are not used to base predictions on. Hence, we need to
ensure that, within a family, no events that happen after the end of
follow up of the proband is used to estimate the genetic liability of
the proband. In real world analysis, the end of follow up can be due to
the proband being diagnosed or any censoring event, such as end of
register follow up, emigration, or death. The function
[`familywise_censoring_multi()`](https://emilmip.github.io/LTFGRS/reference/familywise_censoring_multi.md)
is a wrapper function around
[`familywise_censoring()`](https://emilmip.github.io/LTFGRS/reference/familywise_censoring.md).
It applied the the proband’s end of follow up to all family members and
all outcomes. As such, it offers a way to censor future events on a
family basis, by censoring all events that happen after the end of
follow up of the proband (`indiv_eof`) of the **target outcome**,
specified with the `target outcome` argument.

``` r
# calculate family specific thresholds

info = familywise_censoring_multi(
  family_graphs = family_graphs,
  tbl = pheno %>% rename(pid = id),
  start = "fdato",
  end_base = "indiv_eof",
  phen_names = all_outcomes,
  target_outcome = "disorder_1",
  status_col_base = "status",
  aod_col_base = "aod",
  age_eof_col_base = "age_eof",
  fam_graph_col = "fam_graph",
  fid = "fid",
  pid = "pid",
  simplify = TRUE # simplifies output by removing all columns not specific to an outcome.
)

paged_table(info)
```

The simplified return of function
[`familywise_censoring_multi()`](https://emilmip.github.io/LTFGRS/reference/familywise_censoring_multi.md)
will return a tibble with the following additional columns for **each
outcome**, following the naming convention of
`<base_string>_<disorder_name>`. The end of follow-up applied across a
family (and all disorders) is the one specified by `target_outcome`.
Hence, events that happen in the proband for a different outcome than
the `target_outcome`, but before the end of follow-up for the
`target_outcome`, will not be censored. The columns generated for a
given outcome are as follows:

- `status`: Assigned case-control status based on the family-wise
  censoring time, i.e. if the event happened after the end of follow up
  of the family, the status is set to `FALSE` (0).
- `aod`: The age of diagnosis, NA for controls.
- `age`: The age at the end of the family-wise follow up.

The above information is used to calculate the personalised thresholds
for each individual while also accounting for the fact that each family
may have differing levels of information available to them.

### Assign thresholds to censored families

Once the `status`, `aod`, and `age` are known for each outcome, we can
assign thresholds to each family and their family members. Notably, this
means that if an individual appears in multiple families, e.g. as a
proband and as a relative to a different proband, that individual may
have multiple (potentially different) thresholds assigned to them. The
use of `fid` and `pid` helps ensure that each individual can still be
uniquely identified.

Due to data privacy, it is possible to encounter CIPs values that are
only provided at set values, e.g. a CIP value for each whole year by
birth year and sex, such as what is shown in with the `CIP` object.
However, the observed ages (or age of diagnosis) are typically not
integer values. This means we may need to approximate the CIP values
between the provided values. We offer an XGboost based approach to
interpolate the CIPs between the provided values.

``` r
# thresholds: personalised, fixed upper and lower threshold
multi_thrs = prepare_thresholds_multi(
  .tbl = info,
  CIP_list = CIP_list, 
  phen_names = all_outcomes, 
  personal_thr = TRUE, 
  lower_equal_upper = TRUE
)
```

    ## Warning in prepare_thresholds(.tbl = .tbl %>% rename(`:=`(!!as.symbol(age_col),
    ## : prepare_thresholds: Some ages are negative. This may be due to the end of
    ## follow-up happening before the birth of an individual. Setting their thresholds
    ## to be uninformative. Please check the input data.
    ## Warning in prepare_thresholds(.tbl = .tbl %>% rename(`:=`(!!as.symbol(age_col),
    ## : prepare_thresholds: Some ages are negative. This may be due to the end of
    ## follow-up happening before the birth of an individual. Setting their thresholds
    ## to be uninformative. Please check the input data.
    ## Warning in prepare_thresholds(.tbl = .tbl %>% rename(`:=`(!!as.symbol(age_col),
    ## : prepare_thresholds: Some ages are negative. This may be due to the end of
    ## follow-up happening before the birth of an individual. Setting their thresholds
    ## to be uninformative. Please check the input data.

Note: Currently, it is possible to experience negative ages at the end
of follow up. This is due to the family-wise end of follow up ending
before an individual is born, e.g. proband is diagnosed in their youth,
then has a child later in life. These individuals will get a
non-informative threshold ($- \infty$ to max of `min_CIP_value` and the
predicted CIP \[`K_i`\]). In other words, these individuals will have no
impact on the estimate.

The function
[`prepare_thresholds_multi()`](https://emilmip.github.io/LTFGRS/reference/prepare_thresholds_multi.md)
is a wrapper around
[`prepare_thresholds()`](https://emilmip.github.io/LTFGRS/reference/prepare_thresholds.md)
and as such handles the application of
[`prepare_thresholds()`](https://emilmip.github.io/LTFGRS/reference/prepare_thresholds.md)
and formatting of its result into one tibble. The columns generated for
a given outcome follow the naming convention of
`<base_string>_<disorder_name>`.
[`prepare_thresholds()`](https://emilmip.github.io/LTFGRS/reference/prepare_thresholds.md)
has several options that are worth pointing out. The first is
`lower_equal_upper`, which is used to determine if the `upper` and
`lower` thresholds should be the same for cases or not. This may be
useful if the CIP values are considered very accurate, as it may lead to
more accurate genetic liability estimates. The second is `personal_thr`,
which specifies if thresholds should be based on `K_i` or `K_pop`.
Basing the thresholds on `K_i` yields personalised thresholds that are
based on the stratification of the CIPs. With the argument `Kpop`, it is
possible to determine how the `K_pop` values are calculated. The current
default option for `Kpop` is `"useMax`, which calculates `K_pop` as the
maximum within each strata provided in the CIPs. Alternatively, a tibble
can be provided with the `Kpop` argument, which shares columns with
`.tbl`, e.g. sex, such that user-specific `K_pop` values can be
specified.

The function
[`prepare_thresholds_multi()`](https://emilmip.github.io/LTFGRS/reference/prepare_thresholds_multi.md)
will return a tibble with the following additional columns for **each
outcome**, following the naming convention of
`<base_string>_<disorder_name>`:

- `K_i`: The CIP value for the individual. `K_i` is predicted if
  interpolation is used.
- `K_pop`: The population prevalence. Currently calculated as the
  maximum CIP value within the CIP stratum an individual belongs to,
  e.g. for a male born in $2000$, `K_pop` is the maximum CIP value
  observed among males born in the year $2000$. Alternatively, acquired
  user-specified values through the `Kpop` argument.
- `thr`: The liability threshold used to determine case-control status.
  `thr` is used to determine the upper and lower thresholds of an
  individual.
- `lower`: lower threshold of an individual.
- `upper`: upper threshold of an individual.

If the mixture model is not used to calculate the genetic liability,
only `lower` and `upper` are needed.

``` r
paged_table(multi_thrs)
```

### Assign thresholds to family graphs

The function
[`estimate_liability()`](https://emilmip.github.io/LTFGRS/reference/estimate_liability.md)
can use the family graphs to calculate all necessary values for the
genetic liability, if the thresholds are stored as attributes in the
family graphs. We can attach the thresholds from `multi_thrs` to
`family_graphs` with the function
[`familywise_attach_attributes()`](https://emilmip.github.io/LTFGRS/reference/familywise_attach_attributes.md).
The function merges on the `fid` column and attaches any columns in
`fam_attr` that are specified in `cols_to_attach`. Since the purpose of
the genetic liability we are estimating is prediction, we do not wish to
use the information from the proband. We mask the proband’s information
by setting column with `upper` and `lower` in their name to $\infty$ and
$- \infty$ (non-informative values) and all other values to be `NA`. The
argument `proband_cols_to_censor` can then be used to censor a subset or
all columns. If columns with `K_i` or `K_pop` in the name are in
`proband_cols_to_censor`, they will be set to `NA`, which are
uninformative in the PA algorithm and will not be problematic.

``` r
# attach family specific thresholds
ltfgrs_input = familywise_attach_attributes(
  family_graphs = family_graphs,
  fam_attr = multi_thrs,
  fam_graph_col = "fam_graph",
  attached_fam_graph_col = "masked_fam_graph",
  cols_to_attach = str_subset(colnames(multi_thrs), "lower|upper|K_i|K_pop"),
  proband_cols_to_censor = str_subset(colnames(multi_thrs), "lower|upper|K_i|K_pop"))
ltfgrs_input %>% print(n = 4)
```

    ## # A tibble: 31 × 2
    ##   fid   masked_fam_graph
    ##   <chr> <list>          
    ## 1 pid   <igraph>        
    ## 2 dad   <igraph>        
    ## 3 mom   <igraph>        
    ## 4 dad2  <igraph>        
    ## # ℹ 27 more rows

The format is similar to the one used in
[`get_family_graphs()`](https://emilmip.github.io/LTFGRS/reference/get_family_graphs.md).
It is worth noting that the argument `proband_cols_to_censor` is only
required if the purpose of the resulting genetic liability is
prediction.

``` r
ltfgrs_input$masked_fam_graph[[1]]
```

    ## IGRAPH f97b4a6 DN-- 4 6 -- 
    ## + attr: name (v/c), K_i_disorder_1 (v/n), K_pop_disorder_1 (v/n),
    ## | lower_disorder_1 (v/n), upper_disorder_1 (v/n), K_i_disorder_2 (v/n),
    ## | K_pop_disorder_2 (v/n), lower_disorder_2 (v/n), upper_disorder_2
    ## | (v/n), K_i_disorder_3 (v/n), K_pop_disorder_3 (v/n), lower_disorder_3
    ## | (v/n), upper_disorder_3 (v/n)
    ## + edges from f97b4a6 (vertex names):
    ## [1] dad->sib dad->pid mom->sib mom->pid sib->pid pid->sib

The second row of the igraph (starting with “+attr:”) shows the
attributes that are available to each node in the graph.

**OBS:** There is one important decision to make. When calculating a
genetic liability for prediction, we do not want to use future events to
predict past events and we do not want to use the proband’s information
to predict the proband’s genetic liability. Hence, we need to censor the
proband’s information and all information that happens after the end of
follow up of the proband. This is done with the `proband_cols_to_censor`
argument in
[`familywise_attach_attributes()`](https://emilmip.github.io/LTFGRS/reference/familywise_attach_attributes.md).
However, with the introduction of multiple traits, it opens the
opportunity to perform familywise_censoring with respect to a
`target_outcome`, where events in the secondary outcomes (not the one
chosen for `target_outcome`) has events in the proband. There are
currently no clear guidelines on whether these events should be censored
or not. If the events are not censored, it may lead to more accurate
genetic liability estimates, but it may also lead to information leakage
or information that is very closely tied to the target outcome. A
cautious approach would be to censor all events in the proband,
including all events in the secondary outcomes. It may lead to less
accurate genetic liability estimates. The decision on whether to censor
these events or not is left to the discretion of the user. In the
example above, all information on the proband is censored.

## Estimating genetic liabilities with `estimate_liability()`

The function
[`estimate_liability()`](https://emilmip.github.io/LTFGRS/reference/estimate_liability.md)
is used to estimate the genetic liability. The function accepts two
types of input, here we will only focus on the graph-based input
generated above. The graph-based input offer the best flexibility and
scalability. The function has two arguments that are worth pointing out.

The first is `method`, which specifies the estimation method used to
estimate the genetic liability. Currently, two methods are supported.
The first is a Gibbs sampler that samples from a truncated multivariate
normal distribution, `method = "Gibbs"`. The second is an iterative
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
ltfgrs_pa = estimate_liability(family_graphs = ltfgrs_input,
                                h2 = h2_vec,
                               genetic_corrmat = corMat,
                               full_corrmat = diag(Ntrait),
                               phen_names = all_outcomes,
                               fid = "fid",
                               pid = "pid",
                               family_graphs_col = "masked_fam_graph",
                               method = "PA", # <- METHOD
                               target_phenotype = "disorder_1",
                               useMixture = F)
```

    ## The number of workers is 1

``` r
paged_table(ltfgrs_pa)
```

When using `method = "PA"`, an iterative conditioning is performed,
which means the resulting estimate and uncertainty of the estimate is
the expected mean value and variance of the last iteration, which is the
proband’s genetic liability. This is highlighted by the use of `var` in
the output.

***OBS:*** The PA estimation approach currently only supports returning
a genetic liability estimate for one of the provided outcomes at a time.
The disorder that is returned in specified with the argument
`target_phenotype`. The Gibbs sampling approach will return a genetic
liability estimate for all provided outcomes at the same time. This
means that if you want to get a genetic liability estimate for multiple
outcomes, you will need to run the PA estimation approach multiple
times, once for each outcome, while the Gibbs sampling approach can be
run once to get estimates for all outcomes at the same time.

### FGRS with personalised threshold and Gibbs

``` r
ltfgrs_gibbs = estimate_liability(family_graphs = ltfgrs_input,
                                  h2 = h2_vec,
                                  genetic_corrmat = corMat,
                                  full_corrmat = diag(Ntrait),
                                  phen_names = all_outcomes,
                                  fid = "fid",
                                  pid = "pid",
                                  family_graphs_col = "masked_fam_graph",
                                  method = "Gibbs", # <- METHOD
                                  useMixture = F)
```

    ## The number of workers is 1

``` r
paged_table(ltfgrs_gibbs)
```

When using `method = Gibbs`, samples are drawn from a truncated
multivariate normal distribution until the convergence criteria is met.
The output is the mean of genetic liability and the uncertainty is the
standard error of the mean. The standard error is denoted with the `se`
in the column name. Notably, the uncertainties with the Gibbs and PA
methods are different and should not be directly compared.

### FGRS with mixture model

The mixture model is currently only supported with the PA estimation
method. The mixture model considers the genetic liability of controls as
a mixture of the truncated normal for cases and controls, rather than
just the distribution of controls. This accounts for the possibility
that some controls are undiagnosed cases and accounts for it in the
genetic liability estimate. Below is a full run of preparation for a
mixture model application. In the data preparation steps, set
`personal_thr = FALSE` and `lower_equal_upper = FALSE` when calculating
thresholds. Make sure the `K_i` and `K_pop` are attached to the family
graph too. Finally, set `useMixture = TRUE` in
[`estimate_liability()`](https://emilmip.github.io/LTFGRS/reference/estimate_liability.md).

``` r
# creating a graph for the family with mixture model thresholds
graph_mixture = prepare_graph(.tbl = family, 
                              icol = "id",
                              mcol = "momcol",
                              fcol = "dadcol")

# Identify family members of degree n
family_graphs_mixture = get_family_graphs(pop_graph = graph_mixture,
                                          ndegree = 1,
                                          proband_vec = pheno$id,
                                          fid = "fid",
                                          fam_graph_col = "fam_graph")


info_mixture = familywise_censoring_multi(
  family_graphs = family_graphs_mixture,
  tbl = pheno %>% rename(pid = id),
  start = "fdato",
  end_base = "indiv_eof",
  phen_names = all_outcomes,
  target_outcome = "disorder_1",
  status_col_base = "status",
  aod_col_base = "aod",
  age_eof_col_base = "age_eof",
  fam_graph_col = "fam_graph",
  fid = "fid",
  pid = "pid",
  simplify = TRUE # simplifies output by removing all columns not specific to an outcome.
)


# thresholds: personalised, fixed upper and lower threshold
multi_thrs_mixture = prepare_thresholds_multi(
  .tbl = info_mixture, 
  CIP_list = CIP_list, 
  phen_names = all_outcomes, 
  personal_thr = FALSE, #<--
  lower_equal_upper = FALSE #<---
)


# attaching familywise censored thresholds to family graphs
ltfgrs_input = familywise_attach_attributes(
  family_graphs = family_graphs_mixture,
  fam_attr = multi_thrs_mixture,
  fam_graph_col = "fam_graph",
  attached_fam_graph_col = "masked_fam_graph_mixture",
  cols_to_attach = str_subset(colnames(multi_thrs), "lower|upper|K_i|K_pop"), #<-- includes K_i and K_pop
  proband_cols_to_censor = str_subset(colnames(multi_thrs), "lower|upper|K_i|K_pop")) #<--


PA_MT_mix = estimate_liability(family_graphs = ltfgrs_input,
                               family_graphs_col = "masked_fam_graph_mixture",
                               h2 = h2_vec,
                               genetic_corrmat = corMat,
                               phen_names = all_outcomes,
                               full_corrmat = diag(Ntrait),
                               target_phenotype = "disorder_1", 
                               useMixture = TRUE)#<-- using mixture model
paged_table(PA_MT_mix)
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
