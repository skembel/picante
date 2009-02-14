### Peter's tip coloring

`color.plot.phylo` <- 


function(
    # this function takes a phylogeny and a dataframe
    # with the trait of interest and the taxa columns
    # explicitly passed
    # The number of colors defaults to the number of
    # factors (this cannot be chnaged safely) for factors
    # and 12 for continous traits
    # A user supplied vector of colors can be passed to
    # col.names, be sure this equals num.breaks
    # leg.title provides a title for the legend
                    phylo, df, trait, taxa.names,
                    num.breaks = ifelse(is.factor(df[,trait]),
                        length(levels(df[,trait])), 12),
                    col.names = rainbow(ifelse(length(num.breaks) > 1, length(num.breaks) - 1, num.breaks)),
                    cut.labs = NULL,
                    leg.title = NULL,
                    main = trait,
                    leg.cex = 1,
                    tip.labs = NULL,
                    ...
)
{
    # Get the initial par settings so they can be restored
    init.par <- par(mar = c(0, 0, 1, 0))
    # some data input error checking, all taxa in tree and df
    # no missing data values
    stopifnot( trait %in% names(df), taxa.names %in% names(df),
        class(df) == "data.frame", class(phylo) == "phylo")
    len.tips <- length(phylo$tip.label)
    len.taxa <- length(df[,taxa.names])
    if (len.tips != len.taxa |
        sum(phylo$tip.label %in% df[,taxa.names]) != len.taxa) {
            stop("ERROR. Missing taxa in tree or data frame; # tips: ",
                len.tips, "# taxa: ", len.taxa, "# tips in df: ",
                sum(phylo$tip.label %in% df[,taxa.names]))
        }
    # ensure that the order of the data frame matches the tips
    order <- match(phylo$tip.label, df[,taxa.names])
    ordered.trait <- df[trait][order,]
    # cut up the trait and assign a list of colors
        if(is.factor(ordered.trait)){
            levs <- levels(ordered.trait)
            tip.color <- rep("black", times = len.taxa)
            tip.color <- col.names[match(ordered.trait, levs)]
        }else{
            tip.color = as.character(cut(
                ordered.trait,
                breaks = num.breaks,
                labels = col.names
            ))
            # levs gets used in legend, make one for continous data
            levs <- levels(cut(ordered.trait, breaks = num.breaks))
        }
    if(!is.null(tip.labs)) {
        phylo$tip.label <- df[tip.labs][order,]
    }
    plot.phylo(
        phylo,
        ## y.lim = c(0,80),
        cex = .8,
        tip.color = tip.color,
        main = main,
        ...
    )
    title(line = 0)
    
    if(is.null(cut.labs)) cut.labs <- levs
    legend(
        "bottomleft", cut.labs,
        fill = col.names,
        inset = 0.05, title = leg.title,
        cex = leg.cex
    )
    # restore settings ?
    on.exit(par(init.par))
}
