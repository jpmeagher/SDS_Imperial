library(ancestralbat)
library(ape)

x11(height = 5)
plot.phylo(mexicantree)
axisPhylo()
title(main = list('Mexican Bat Phylogeny', font = 1, cex = 1.5), xlab = 'Time (millions of years)')
