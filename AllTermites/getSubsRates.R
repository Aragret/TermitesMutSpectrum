library(ape)

mut = read.table('results/4foldSubs.txt', header = TRUE, sep='\t')

tree = read.tree('results/phylogeny/iqtree_mito.treefile')

numberOfSpecies = length(tree$tip.label)

a = as.data.frame(tree$edge)
a = cbind(a, tree$edge.length)
externalBranches = a[a$V2 <= numberOfSpecies,]
externalBranches = cbind(externalBranches, tree$tip.label)

names(externalBranches) = c('Node', 'Tip', 'BranchLength', 'Species')

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

a = matrix(0, ncol=12, nrow=speciesNumber)

for(i in 1:12){
  a[, i] = data[, i+1] / data$BranchLength
}

rates = as.data.frame(a)
names(rates) = sub(' ', '', paste(names(mutTable[, -1]), '_rate'))

subsRates = cbind(data, rates)

write.table(subsRates, 'results/4foldSubsRates.txt', sep='\t', row.names = FALSE, quote = FALSE)
