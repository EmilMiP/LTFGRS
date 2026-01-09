# LTFGRS

LTFGRS is an R package that provides a unified interface to perform
several phenotype-based family genetic risk scores (FGRS) as well as
functions to perform the required data preparation and management steps.

The implemented methods are
[LT-FH++](https://doi.org/10.1016/j.ajhg.2022.01.009),
[PA-FGRS](https://pubmed.ncbi.nlm.nih.gov/39471805/), and [Kendler’s
FGRS](https://pubmed.ncbi.nlm.nih.gov/33881469/). The primary focus of
the package is LT-FH++ and PA-FGRS. They share many similarities in
their models, since they are both based on the same underlying liability
threshold model with modifications. This package allows users to pick
and chose which elements of those methods they want to use. Notable
properties of the package are:

- Automatic identification of family members up to n’th degree.
- personalised thresholds for each included family that can account for:
  - age of diagnosis
  - censoring
  - cohort effects
  - sex
  - family-wise censoring
- mixture distribution to handle partially observed controls.

# Installation

You can install the development version of LTFGRS from GitHub with:

    remotes::install_github("EmilMiP/LTFGRS")

You can install the LTFGRS from CRAN with:

    install.packages("LTFGRS")

# Documentation and Tutorials

Documentation and tutorials for the package can be found in the [pkgdown
site](https://emilmip.github.io/LTFGRS/).

# Citation

If you use LTFGRS to estimate genetic liabilities, please cite at least
one of the following papers depending on which method you used:

[Pedersen et al, Accounting for age of onset and family history improved
power in genome-wide association studies. AJHG,
2022.](https://doi.org/10.1016/j.ajhg.2022.01.009)

[Krebs et al., Genetic liability estimated from large-scale family data
improves genetic prediction, risk score profiling, and gene mapping for
major depression. AJHG,
2024](https://doi.org/10.1016/j.ajhg.2024.09.009)

[Kendler et al., Family Genetic Risk Scores and the Genetic Architecture
of Major Affective and Psychotic Disorders in a Swedish National Sample.
JAMA Psychiatry,
2021.](https://doi.org/10.1001/jamapsychiatry.2021.0336)

If automatic identification of family members is used, please cite:

[Pedersen et al., Automatic detection of n-degree family members.
Frontiers in Genetics,
2025.](https://doi.org/10.3389/fgene.2025.1708315)
