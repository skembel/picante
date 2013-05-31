\name{evol.distinct}
\Rdversion{1.1}
\alias{evol.distinct}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Species' evolutionary distinctiveness
}
\description{
Calculates evolutionary distinctiveness measures for a suite of species by: a) equal splits (Redding and Mooers 2006) b) fair proportions (Isaac et al., 2007). Returns a datafram with species identifiers and species scores.
}
\usage{
evol.distinct(tree, type = c("equal.splits", "fair.proportion"),
    scale = FALSE, use.branch.lengths = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
an object of class phylo}
  \item{type}{
a) equal splits (Redding and Mooers 2006) or b) fair proportions (Isaac et al., 2007)
}
  \item{scale}{
The scale option refers to whether or not the phylogeny should be scaled to a depth of 1 or, in the case of an ultrametric tree,  scaled such that branch lengths are relative.
}
  \item{use.branch.lengths}{
If use.branch.lengths=FALSE, then all branch lengths are changed to 1.
}
}

\references{
Redding, D.W. and Mooers, A.O. (2006). Incorporating evolutionary measures into conservation prioritisation. Conservation Biology, 20, 1670-1678.

Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.

Mark Vellend, William K. Cornwell, Karen Magnuson-Ford, and Arne Mooers. In press. Measuring phylogenetic biodiversity. In: Biological diversity: frontiers in measurement and assessment.  Edited by Anne Magurran and Brian McGill.
}

\author{
Karen Magnuson-Ford, Will Cornwell, Arne Mooers, Mark Vellend
}

\note{
This function will return a vector of evolutionary distinctivenss for every species in the given tree.
If only a subset of values are needed there are two, concetually distinct options: either prune the tree first and then pass the tree in or subset the resulting vector.  These two options will provide very different outputs.  
}


