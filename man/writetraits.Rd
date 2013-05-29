\name{writetraits}
\alias{writetraits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write a Phylocom traits formatted file}
\description{
  Write a Phylocom traits formatted file
}
\usage{
writetraits(trt, file = "", bin = NULL, sigd = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trt}{ Data frame containing trait data }
  \item{file}{ Filename path }
  \item{bin}{ Vector index of trait columns to be treated as binary }
  \item{sigd}{ Significant digits for output }
}

\references{ Webb, C.O., Ackerly, D.D., and Kembel, S.W. 2008. Phylocom: software for the analysis of phylogenetic community structure and trait evolution. Version 4.0.1. \url{http://www.phylodiversity.net/phylocom/}. }
\author{ David Ackerly <dackerly@berkeley.edu> and Steven Kembel <steve.kembel@gmail.com> }
\keyword{file}
