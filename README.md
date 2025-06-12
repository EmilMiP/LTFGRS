# LTFGRS

LTFGRS is an R package that provides a unified interface to perform several phenotype-based family genetic risk scores (FGRS) as well as functions to perform the required data preparation and management steps.

The implemented methods are [LT-FH++](https://doi.org/10.1016/j.ajhg.2022.01.009), [PA-FGRS](https://pubmed.ncbi.nlm.nih.gov/39471805/), and [Kendler's FGRS](https://pubmed.ncbi.nlm.nih.gov/33881469/). 
The primary focus of the package is LT-FH++ and PA-FGRS. They share many similarities in their models, since they are both based on the same underlying liability threshold model with modifications. This package allows users to pick and chose which elements of those methods they want to use.
Notable properties of the package are:

- Automatic identification of family members up to n'th degree.
- personalised thresholds for each included family that can account for:
  - age of diagnosis
  - censoring
  - cohort effects
  - sex
  - family-wise censoring
- mixture distribution to handle partially observed cases.

# Installation
You can install the development version of LTFGRS from GitHub with:

```{r, eval=FALSE}
remotes::install_github("EmilMiP/LTFGRS")
```



<!-- badges: start -->
[![R-CMD-check](https://github.com/EmilMiP/LTFGRS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EmilMiP/LTFGRS/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->


to-do:

- [x] Move parts of LTFHPlus that should be reused to this package
- [x] Add PA algorithm
- [x] Add kendler
- [x] Ensure PA-FGRS can be run with unified (graph) input
- [x] Add a description of the functions
- [ ] Add (linear) interpolation to input function
- [x] Remove references to LT-FH++ in the code and descriptions
- [x] Rename prepare_LTFHPlus_input function 
- [ ] Add tests for critical functions
- [ ] Formalise README.md 
- [ ] Formalise description of the package
- [ ] Add a vignettes on how to use aspects of the package to perform PA-FGRS, LT-FH++, kendler, etc.
- [ ] Add a vignette on how to use the package to perform simulations
- [ ] Formalise pkgdown site
- [ ] Ensure only required packages are in the DESCRIPTION file
- [x] Add functions that allow us to manage thresholds and censor event
- [x] Add functions that allow us to remove / add / modify exising node attributes
- [x] Add functions that calculate thresholds given CIP data with interpolation options
- [x] Replace 'list.vertex.attributes()' (deprecated) with 'vertex_attr_names()'(new function) in the code
- [x] replace 'get.vertex.attribute()' (deprecated) with 'vertex_attr()' (new function) in the code
- [ ] ensure all functions have at least one example
- [ ] ensure the documentation of estimate_liability() reflects the different names of the returned column names depending on estimation method
- [x] Should Lucas be a contributor on the package?
- [x] Rename prepare_LTFHPlus_input to something else. prepare_thresholds?
- [x] Rename "use_fixed_case_thr" in prepare_LTFHPlus_input to "use_fixed_case_thr" to something more informative + update documentation to reflect what fixed means
- [ ] Add a wrapper function to format input with no family-wise censoring
- [ ] more flexibility to k_pop in prepare_thresholds. Add K_pop = "useMax" or K_pop = tbl (if some specific values are desired) instead?
- [x] rename censor_family_onsets_per_family. familywise_censoring?
- [x] rename assign_family_specific_thresholds or fam_graph_attach_attributes to something else?
- [x] replace all arguments of "fam_id" with "fid" to standardise all inputs.
