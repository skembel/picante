\name{tax.distinctiveness}
\Rdversion{1.1}
\alias{tax.distinctiveness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Taxonomic distinctiveness sensu Vane-Wright or May
}
\description{
Taxic diversity: Vane-Wright et al., 1991 and May 1990 which accounts for polytomies by counting the number of branches descending from each node that lies on the path from a spp tip to the root (not just counting the number of nodes).
}
\usage{
tax.distinctiveness(tree, type = c("Vane-Wright", "May"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
an object of class phylo
}
  \item{type}{
specify "Vane-Wright" or "May"
}
}


\references{
Vane-Wright, R.I., Humphries, C.J. and Williams, P.H. (1991). What to protect? - Systematics and the agony of choice. Biological Conservation, 55, 235-254.

May, R.M. (1990). Taxonomy as destiny. Nature, 347, 129-130.

Mark Vellend, William K. Cornwell, Karen Magnuson-Ford, and Arne Mooers. In press. Measuring phylogenetic biodiversity In: Biological diversity: frontiers in measurement and assessment. Edited by Anne Magurran and Brian McGill.
}

\author{
Karen Magnuson-Ford, Will Cornwell, Arne Mooers, Mark Vellend
}

