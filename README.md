PICANTE
================
Steven Kembel

=======

## News

There was a bug in the ses.pd function in picante versions \<=1.7.
Please see [this issue](https://github.com/skembel/picante/issues/17)
for more details.

## Overview

R tools for integrating phylogenies and ecology

The PICANTE package provides tools for **P**hylocom **I**ntegration,
**C**ommunity **A**nalyses, **N**ull-models, **T**raits and
**E**volution in R.

The package includes functions for analyzing the phylogenetic and trait
diversity of ecological communities, comparative analyses, and the
display and manipulation of phenotypic and phylogenetic data.

## Installation

To install the package

``` r
install.packages("picante")
```

To import the package

``` r
library('picante')
```

To import a dataset

``` r
data(phylocom)
```

## Exemple of datasets

### Community

Picante uses the same community data format as the vegan package - a
matrix or data.frame with sites/samples in the rows and taxa in the
columns. The elements of this data frame should be numeric values
indicating the abundance or presence/absence (0/1) of taxa in different
samples.

One important thing to note is that most functions in picante will use
the labels on columns in the community data set to match community and
phylogenetic data together. You need to make sure your column names are
present and match the tip labels of the phylogeny or your analysis may
not work.

Most functions in picante do basic error checking and will report when
there are mismatches between the data present in the community and
phylogenetic data sets. Similarly, your communities/sites/samples can be
given informative names, and these should be contained in the row
labels, not in a column of the data.frame.

<style>
body {
text-align: justify}
</style>

``` r
community <-phylocom$sample
community[,1:10]
```

    ##         sp1 sp10 sp11 sp12 sp13 sp14 sp15 sp17 sp18 sp19
    ## clump1    1    0    0    0    0    0    0    0    0    0
    ## clump2a   1    2    2    2    0    0    0    0    0    0
    ## clump2b   1    0    0    0    0    0    0    2    2    2
    ## clump4    1    1    0    0    0    0    0    2    2    0
    ## even      1    0    0    0    1    0    0    1    0    0
    ## random    0    0    0    1    0    4    2    3    0    0

### Phylogenies

Picante uses the phylo format implemented in the ape package to
represent phylogenetic relationships among taxa.

``` r
## in this exemple, edge have been shorted an transposed for better visibility
phy <- phylocom$phylo
phy$edge <- t(phy$edge[1:10,])
phy
```

    ## 
    ## Phylogenetic tree with 32 tips and 31 internal nodes.
    ## 
    ## Tip labels:
    ##   sp1, sp2, sp3, sp4, sp5, sp6, ...
    ## Node labels:
    ##   A, B, C, D, E, F, ...
    ## 
    ## Rooted; includes branch lengths.

***Important*** : in this example, edge have been shorted an transposed
for better visibility

### Trait

Trait data include any kind of data associated with the taxa present in
a phylogeny. Most functions in picante work with trait data represented
as a vector (for individual traits) or data.frame. The documentation for
individual functions will explain which data format is expected.

``` r
traits <- phylocom$traits
head(traits)
```

    ##     traitA traitB traitC traitD
    ## sp1      1      1      1      0
    ## sp2      1      1      2      0
    ## sp3      2      1      3      0
    ## sp4      2      1      4      0
    ## sp5      2      2      1      0
    ## sp6      2      2      2      0

## Usage

To calculates mean pairwise distance (MPD) separating taxa in a
community

``` r
picante::mpd(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
```

    ##   clump1  clump2a  clump2b   clump4     even   random 
    ## 4.250000 4.944444 5.833333 6.944444 7.750000 7.109375

To calculates mean nearest taxon distance (MNTD) for taxa in a community

``` r
mntd(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)
```

    ## [1] 2.000 2.000 2.000 2.000 6.000 4.875

## Authors and contributors

**Package maintainer**: Steven Kembel

**Developers**: Peter Cowan, Matthew Helmus, Steven Kembel

**Contributors**: David Ackerly, Simon Blomberg, Will Cornwell, Peter
Cowan, Matthew Helmus, Steven Kembel, Helene Morlon, Cam Webb

**Supporters** Development of picante has been supported by NSERC,
NESCent, the Google Summer of Code, and the Gordon and Betty Moore
Foundation.

Thanks to Jonathan Davies, Kyle Dexter, Catherine Graham, Nathaniel
Hallinan,Nick Matzke, Alain Paquette, Emmanuel Paradis, Juan Parra, Dan
Rabosky, and Marten Winter for feedback and bug reports.

Thanks to R-Forge for hosting the project up to version 1.5.

## Bug reports and feature requests

Please use the Github Issues page for picante to submit bug reports and
feature requests. Pull requests are welcomed but are more likely to be
accepted if you first submit them as an issue.
