library(phytools)
library(stringr)
library(ggstance)


tree = read.tree('results/rootedTree.newick')
mut = read.table('results/4foldSubsRatesBrLen02Residuals.txt', header = TRUE, sep='\t')
row.names(mut) = mut$Species

tree$tip.label = str_replace_all(tree$tip.label, fixed('__'), '_')

a = setdiff(tree$tip.label, mut$Species)
tree2 = drop.tip(tree, a)

write.tree(tree2, 'results/tree02BrLen_mod.newick')

setdiff(tree2$tip.label, mut$Species)
setdiff(mut$Species, tree2$tip.label)

pdf('results/heatmapTree.pdf', width = 500, height = 1000)

obj<-contMap(tree2, mut[, 2:13], plot=FALSE)

phylo.heatmap(tree, mut[, 2:13], standardize=TRUE)

plotTree.barplot(tree2, mut$C_T, xlim=c(0,100))


dev.off()


d3 <- data.frame(id = rep(tr$tip.label, each=2),
                 value = abs(rnorm(60, mean=100, sd=50)),
                 category = rep(LETTERS[1:2], 30))

p3 <- facet_plot(tree2, panel = 'Stacked Barplot', data = mut, 
                 geom = geom_barh, 
                 mapping = aes(x = C_T), 
                 stat='identity' ) 
