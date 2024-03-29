% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{update_corecmotif}
\alias{update_corecmotif}
\title{Update a CoRecMotif}
\usage{
update_corecmotif(corecmotif, keep_match = FALSE)
}
\arguments{
\item{corecmotif}{\link[=CoRecMotif-class]{CoRecMotif}. The
\link[=CoRecMotif-class]{CoRecMotif} to update.}

\item{keep_match}{\code{logical(1)}. Should the \code{match_motif}, \code{match_pvalue}, and
\code{match_cluster} slots be kept? If \code{FALSE}, they will be reset to \code{NA}.
(Default: FALSE)}
}
\value{
An updated \link[=CoRecMotif-class]{CoRecMotif}.
}
\description{
Updates all the values that are calculated automatically when a
\link[=CoRecMotif-class]{CoRecMotif} is created.
}
\examples{
print("FILL THIS IN")
}
