library(ape)

mut = read.table('results/4foldSubs.txt', header = TRUE, sep='\t')

tree = read.tree('results/phylogeny/iqtree_mito.treefile')

numberOfSpecies = length(tree$tip.label)

is.rooted(tree)
rootedTree = root(tree, 'Locusta_migratoria_1', resolve.root	= TRUE)
write.tree(rootedTree, 'results/rootedTree.newick')

a = as.data.frame(rootedTree$edge)
a = cbind(a, rootedTree$edge.length)
externalBranches = a[a$V2 <= numberOfSpecies,]
externalBranches = cbind(externalBranches, tree$tip.label)

names(externalBranches) = c('Node', 'Tip', 'BranchLength', 'Species')

brLenLess05 = externalBranches[externalBranches$BranchLength < 0.05,]
brLenMore05 = externalBranches[externalBranches$BranchLength > 0.05,]

brLenLess02 = externalBranches[externalBranches$BranchLength < 0.02,]
brLenMore02 = externalBranches[externalBranches$BranchLength > 0.02,]

treeMore05 = drop.tip(rootedTree, as.character(brLenLess05$Species))
length(treeMore05$tip.label)

treeMore02 = drop.tip(rootedTree, as.character(brLenLess02$Species))
length(treeMore02$tip.label)

########################

mut$Subs = as.character(mut$Subs)

speciesNumber = length(unique(mut$Species))
mutTable = setNames(data.frame(matrix(0, ncol = 12, nrow = speciesNumber)), unique(mut$Subs))
Species = unique(mut$Species)
mutTable = cbind(Species, mutTable)


for(i in 1:nrow(mut)){
  row = mut[i, ]
  mutTable[mutTable$Species == row$Species, row$Subs] = mutTable[mutTable$Species == row$Species, row$Subs] + 1
}

table(mut[mut$Species == "AUS49_Tumulitermes_sp._1",]$Subs)


data = merge(mutTable, externalBranches[, c('Species', 'BranchLength')], by='Species')

dataMore05 = data[data$BranchLength > 0.05, ]
dataMore02 = data[data$BranchLength > 0.02, ]

a = matrix(0, ncol=12, nrow=length(unique(dataMore02$Species)))

for(i in 1:12){
  a[, i] = dataMore02[, i+1] / dataMore02$BranchLength
}

rates = as.data.frame(a)
names(rates) = sub(' ', '', paste(names(mutTable[, -1]), '_rate'))

subsRates = cbind(dataMore02, rates)

write.table(subsRates, 'results/4foldSubsRatesBrLen02.txt', sep='\t', row.names = FALSE, quote = FALSE)
