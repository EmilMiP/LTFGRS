% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_to_trio.R
\name{graph_to_trio}
\alias{graph_to_trio}
\title{Convert from igraph to trio information}
\usage{
graph_to_trio(
  graph,
  id = "id",
  dadid = "dadid",
  momid = "momid",
  sex = "sex",
  fixParents = TRUE
)
}
\arguments{
\item{graph}{An igraph graph object.}

\item{id}{Column of proband id. Defaults to id.}

\item{dadid}{Column of father id. Defaults to dadid.}

\item{momid}{Column of mother id. Defaults to momid.}

\item{sex}{Column of sex in igraph attributes. Defaults to sex.}

\item{fixParents}{Logical. If TRUE, the kinship2's fixParents will be run on the trio information before returning. Defaults to TRUE.}
}
\value{
A tibble with trio information.
}
\description{
This function converts an igraph object to a trio information format.
}
\details{
The sex column is required in the igraph attributes. The sex information is used to determine who is the mother and father in the trio.
}
\examples{
if (FALSE) {

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


thrs =  tibble(
  id = family \%>\% select(1:3) \%>\% unlist() \%>\% unique(),
  lower = sample(c(-Inf, 2), size = length(id), replace = TRUE),
  upper = sample(c(2, Inf), size = length(id), replace = TRUE),
 sex = case_when(
   id \%in\% family$momcol ~ "F",
    id \%in\% family$dadcol ~ "M",
    TRUE ~ NA)) \%>\%
  mutate(sex = sapply(sex, function(x) ifelse(is.na(x),
  sample(c("M", "F"), 1), x)))
graph = prepare_graph(.tbl = family,
icol = "id", fcol = "dadcol", mcol = "momcol", node_attributes = thrs)
}

}
