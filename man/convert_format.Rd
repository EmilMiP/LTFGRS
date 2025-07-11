% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Prepare_Graph_Input.R
\name{convert_format}
\alias{convert_format}
\title{Attempts to convert the list entry input format to a long format}
\usage{
convert_format(family, threshs, personal_id_col = "pid", role_col = NULL)
}
\arguments{
\item{family}{a tibble with two entries, family id and personal id. personal id should end in "_role", if a role column is not present.}

\item{threshs}{thresholds, with a personal id (without role) as well as the lower and upper thresholds}

\item{personal_id_col}{column name that holds the personal id}

\item{role_col}{column name that holds the role}
}
\value{
returns a format similar to \code{prepare_thresholds}, which is used by \code{estimate_liability}
}
\description{
Attempts to convert the list entry input format to a long format
}
\examples{
family <- data.frame(
fid = c(1, 1, 1, 1),
pid = c(1, 2, 3, 4),
role = c("o", "m", "f", "pgf")
)

threshs <- data.frame(
  pid = c(1, 2, 3, 4),
  lower = c(-Inf, -Inf, 0.8, 0.7),
  upper = c(0.8, 0.8, 0.8, 0.7)
)

convert_format(family, threshs)

}
