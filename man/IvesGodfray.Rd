\name{IvesGodfray}
\alias{IvesGodfray}
\docType{data}
\title{ Host-parasitoid food web data }
\description{
  Data on the structure of a host-parasitoid food web from Ives & Godfray (2006). Includes information on phylogenetic covariances among 12 leaf-mining moth hosts and 27 species of parasitoid wasps.
}
\usage{data(phylocom)}
\format{
    A list with three elements:
    \itemize{
      \item{ host }{ Phylogenetic variance/covariance matrix for 12 leaf-mining moth hosts }
      \item{ parasitoid }{ Phylogenetic variance/covariance matrix for 27 species of parasitoid wasps }
      \item{ interactions }{ Matrix describing interactions between hosts and parasitoids }
      }
}

\source{
Ives A.R. & Godfray H.C. (2006) Phylogenetic analysis of trophic associations. The American Naturalist, 168, E1-E14
}

\seealso{ \code{\link{pblm}} }

\keyword{datasets}
