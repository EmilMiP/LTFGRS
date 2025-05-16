To be filled out.


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
- [ ] Remove references to LT-FH++ in the code and descriptions
- [ ] Rename prepare_LTFHPlus_input function 
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
- [ ] Should Lucas be a contributor on the package?
- [ ] Rename "use_fixed_case_thr" in prepare_LTFHPlus_input to "use_fixed_case_thr" to something more informative + update documentation to reflect what fixed means
- [ ] Add a wrapper function to format input with no family-wise censoring
- [ ] more flexibility to k_pop in prepare_LTFHPlus_input. Add K_pop = "useMax" or K_pop = tbl (if some specific values are desired) instead?
- [ ] rename prepare_LTFHPlus_input to something else. prepare_thresholds?
- [ ] rename censor_family_onsets_per_family. familywise_censoring?
- [ ] rename assign_family_specific_thresholds or fam_graph_attach_attributes to something else?
- [ ] replace all arguments of "fam_id" with "fid" to standardise all inputs. 
