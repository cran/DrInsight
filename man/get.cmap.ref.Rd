% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readin.cmap.R
\name{get.cmap.ref}
\alias{get.cmap.ref}
\title{Read In and Process CMap Reference Profiles}
\usage{
get.cmap.ref(cmap.data.path = NULL, probe.to.genes = NULL, drug.info = NULL)
}
\arguments{
\item{cmap.data.path}{The local path and the name of the downloaded data matrix file. The data matrix file should be in txt format.}

\item{probe.to.genes}{The ID converter between Affymetrix probe IDs (the IDs used in CMap data matrix) and official gene symbol.
The packge comes with an embeded probe.to.genes file that can be directly used. User can use their own converter file with two columns named "ID", and "Gene.Symbol'.}

\item{drug.info}{The drug instance information from CMap data. The deafault drug.info file comes with DrInsight package.}
}
\description{
This function allows user to load in the CMap drug rank matrix user downloaded from CMap website: https://portals.broadinstitute.org/cmap/ (data matrix in the "downloads" section).
}
\keyword{CMap}
\keyword{Profiles}
\keyword{Reference}
