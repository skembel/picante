##' Color tip labels based on trait
##' 
##' Plots a phylogeny with tip labels colored to indicate continuous or
##' discrete trait values
##' 
##' If  \code{trait} is a factor then each level of the factor is plotted
##' with the corresponding \code{col.names} value (if \code{length(num.breaks)
##' > length(col.names)} colors are recycled.) If \code{trait} is not a factor
##' then it is assumed to be continuous and \code{trait} is evenly divided into
##' \code{num.breaks} levels.
##' 
##' @param phylo An object of class \code{phylo}
##' @param df A dataframe containing the traits to be plotted
##' @param trait A string representing the name of column in the dataframe to
##' be plotted
##' @param taxa.names A string representing the name of column in the dataframe
##' that contains the names of the taxa
##' @param num.breaks For continuous traits, the number of bins to separate the
##' data into
##' @param col.names A vector of colors to use for tip labels
##' @param leg.title A title for the tip color legend
##' @param main A main title for the plot
##' @param cut.labs A main title for the plot
##' @param leg.cex A main title for the plot
##' @param tip.labs A main title for the plot
##' @param ... Additional argument to pass to the \code{plot.phylo} function
##' @return The command is invoked for its side effect, a plot of the
##' \code{phylo} with tips colored based on \code{trait}
##' @author Peter Cowan <pdc@@berkeley.edu>
##' @keywords color
##' @export color.plot.phylo
`color.plot.phylo` <-
  function( # this function takes a phylogeny and a dataframe
           # with the trait of interest and the taxa columns
           # explicitly passed
           # The number of colors defaults to the number of
           # factors (this cannot be changed safely) for factors
           # and 12 for continous traits
           # A user supplied vector of colors can be passed to
           # col.names, be sure this equals num.breaks
           # leg.title provides a title for the legend
           phylo,
           df,
           trait,
           taxa.names,
           num.breaks = ifelse(is.factor(df[, trait]),
             length(levels(df[, trait])), 12
           ),
           col.names = rainbow(ifelse(length(num.breaks) > 1,
             length(num.breaks) - 1, num.breaks
           )),
           cut.labs = NULL,
           leg.title = NULL,
           main = trait,
           leg.cex = 1,
           tip.labs = NULL,
           ...) {
    # Get the initial par settings so they can be restored
    init.par <- par(mar = c(0, 0, 1, 0))
    # some data input error checking, all taxa in tree and df
    # no missing data values
    stopifnot(
      trait %in% names(df), taxa.names %in% names(df),
      class(df) == "data.frame", class(phylo) == "phylo"
    )
    len.tips <- length(phylo$tip.label)
    len.taxa <- length(df[, taxa.names])
    if (len.tips != len.taxa ||
      sum(phylo$tip.label %in% df[, taxa.names]) != len.taxa) {
      stop(
        "ERROR. Missing taxa in tree or data frame; # tips: ",
        len.tips, "# taxa: ", len.taxa, "# tips in df: ",
        sum(phylo$tip.label %in% df[, taxa.names])
      )
    }
    # ensure that the order of the data frame matches the tips
    order <- match(phylo$tip.label, df[, taxa.names])
    ordered.trait <- df[trait][order, ]
    # cut up the trait and assign a list of colors
    if (is.factor(ordered.trait)) {
      levs <- levels(ordered.trait)
      tip.color <- col.names[match(ordered.trait, levs)]
    } else {
      tip.color <- as.character(cut(
        ordered.trait,
        breaks = num.breaks,
        labels = col.names
      ))
      # levs gets used in legend, make one for continous data
      levs <- levels(cut(ordered.trait, breaks = num.breaks))
    }
    if (!is.null(tip.labs)) {
      phylo$tip.label <- df[tip.labs][order, ]
    }
    plot.phylo(
      phylo,
      cex = .8,
      tip.color = tip.color,
      main = main,
      ...
    )
    title(line = 0)

    if (is.null(cut.labs)) cut.labs <- levs
    legend(
      "bottomleft", cut.labs,
      fill = col.names,
      inset = 0.05, title = leg.title,
      cex = leg.cex
    )
    # restore settings ?
    on.exit(par(init.par))
  }
