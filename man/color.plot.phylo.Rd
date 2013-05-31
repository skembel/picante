\name{color.plot.phylo}
\alias{color.plot.phylo}
\title{Color tip labels based on trait}
\description{
  Plots a phylogeny with tip labels colored to indicate continuous or discrete trait values}
\usage{
  color.plot.phylo(phylo, df, trait, taxa.names,
                    num.breaks = ifelse(is.factor(df[,trait]),
                        length(levels(df[,trait])), 12),
                    col.names = rainbow(ifelse(length(num.breaks) > 1,
                        length(num.breaks) - 1, num.breaks)),
                    cut.labs = NULL,
                    leg.title = NULL,
                    main = trait,
                    leg.cex = 1,
                    tip.labs = NULL,
                    ...)
        }
\arguments{
  \item{phylo}{An object of class \code{phylo}}
  \item{df}{A dataframe containing the traits to be plotted}
  \item{trait}{A string representing the name of column in the dataframe to be plotted}
  \item{taxa.names}{A string representing the name of column in the dataframe that contains the names of the taxa}
  \item{num.breaks}{For continuous traits, the number of bins to separate the data into}
  \item{col.names}{A vector of colors to use for tip labels}
  \item{leg.title}{A title for the tip color legend}
  \item{main}{A main title for the plot}
  \item{cut.labs}{A main title for the plot}
  \item{leg.cex}{A main title for the plot}
  \item{tip.labs}{A main title for the plot}
  \item{...}{Additional argument to pass to the \code{plot.phylo} function}
}
\details{
  If if \code{trait} is a factor then each level of the factor is plotted with the corresponding \code{col.names} value (if \code{length(num.breaks) > length(col.names)} colors are recycled.) If \code{trait} is not a factor then it is assumed to be continuous and \code{trait} is evenly divided into \code{num.breaks} levels.
}
\value{
  The command is invoked for its side effect, a plot of the \code{phylo} with tips colored based on \code{trait}
}

\author{ Peter Cowan <pdc@berkeley.edu> }
\keyword{color}

