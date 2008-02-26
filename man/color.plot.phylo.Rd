\name{color.plot.phylo}
\alias{color.plot.phylo}
\title{Color tip labels based on trait}
\description{
  Plots a phylogeny with tip labels colored to indicate continuous or discrete trait values}
\usage{
  color.plot.phylo(
  phylo, df, trait, taxa.names,
  num.breaks = ifelse(is.factor(df[,trait]),
      length(levels(df[,trait])), 12),
  col.names = rainbow(num.breaks),
  leg.title = NULL,
  main = trait,
  ...)}
\arguments{
  \item{phylo}{ ~~Describe \code{phylo} here~~ }
  \item{df}{ ~~Describe \code{df} here~~ }
  \item{trait}{ ~~Describe \code{trait} here~~ }
  \item{taxa.names}{ ~~Describe \code{taxa.names} here~~ }
  \item{num.breaks}{ ~~Describe \code{num.breaks} here~~ }
  \item{col.names}{ ~~Describe \code{col.names} here~~ }
  \item{leg.title}{ ~~Describe \code{leg.title} here~~ }
  \item{main}{ ~~Describe \code{main} here~~ }
  \item{...}{ ~~Describe \code{...} here~~ }
}
\details{
  ...
}
\value{
  The command is invoked for it's side effect, a plot of the \code{phylo} with tips colored based on \code{trait}
}

\author{ Peter Cowan <pdc@berkeley.edu> }

