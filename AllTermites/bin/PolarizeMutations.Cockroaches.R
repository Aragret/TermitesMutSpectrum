library(ape)
library(stringr)
library(seqinr)

tree = read.tree('../results/cockroaches11_19/phylogeny/iqtree_mito.treefile')

anc = read.table('../results/cockroaches11_19/phylogeny/iqtree_mito.state', header=TRUE)
# anc$Node = sub('Node', '', anc$Node)

agg = aggregate(State ~ Node, anc, paste, collapse = "")

table(anc$Node)

tree$edge
numberOfSpecies = length(tree$tip.label)

a = as.data.frame(tree$edge)
a = cbind(a, tree$edge.length)
externalBranches = a[a$V2 <= numberOfSpecies,]
externalBranches = cbind(externalBranches, tree$tip.label)

names(externalBranches) = c('Node', 'Tip', 'BranchLength', 'Species')

tree = root(tree, 'Locusta_migratoria_1')

labelsTree = as.data.frame(as.character(tree$node.label))
names(labelsTree) = 'NodeLabels'
#labelsTree$Node = lapply(as.character(labelsTree$NodeLabels), strsplit, '/')

a = do.call(rbind, str_split(labelsTree$NodeLabels, '/'))

labelsTree = cbind(labelsTree, a)
names(labelsTree) = c('NodeLabels', 'Node', 'Bootstrap', 'Bootstrap2')

ancSeqs = merge(agg, labelsTree, by='Node')

pdf('~/Desktop/aaa.pdf', width = 20, height = 100)
plot(tree); edgelabels(tree$edge[,1], width = 0.05, height = 0.05); nodelabels(tree$node.label, width = 0.05, height = 0.05, col = 'red')

dev.off()

### Just add 536 ???? (number of species)

ancSeqs$Node = as.numeric(sub('Node', '', ancSeqs$Node))
ancSeqs$Node = ancSeqs$Node + numberOfSpecies

data = merge(ancSeqs, externalBranches, by='Node')

aligns = read.fasta('../results/cockroaches11_19/align_cat.fasta', as.string = TRUE)

a = data.frame(Species=names(aligns), Seqs=unlist(getSequence(aligns, as.string=T)))

data$Species = as.character(data$Species)
a$Species = as.character(a$Species)

# a$Species = sapply(a$Species, function(x) sub('?', '', x))
a$Species = str_replace_all(a$Species, fixed('?'), '')
a$Species = str_replace_all(a$Species, fixed(';'), '')
a$Species = str_replace_all(a$Species, fixed('__'), '_')

data$Species = str_replace_all(data$Species, fixed('__'), '_')

anc_desSeqs = merge(data, a, by='Species')

setdiff(data$Species, anc_desSeqs$Species)

######## get codons

codonTable = as.data.frame(NULL)

for (y in 1:nrow(anc_desSeqs)){
  # y = 1
  mdns <- as.character(anc_desSeqs$State[y])
  msns <- toupper(as.character(anc_desSeqs$Seqs[y]))
  brLength = anc_desSeqs$BranchLength[y]
  s1 <- s2c(mdns)       #conversion of a string into a vector of chars
  s2 <-  s2c(msns)
  cdns1 <- splitseq(s1, 0, 3) #split a sequence into sub-sequences (codons)
  cdns2 <- splitseq(s2, 0, 3)
  
  codonposition <-c(1:length(cdns1))
    
  final <- c()
  #allframes <- c()
  for (x in codonposition){
    # x = 1
    if (cdns1[x] != cdns2[x]){
      taa1 <- translate(s2c(cdns1[x]), frame = 0, sens = "F", numcode = 5, NAstring = "X", ambiguous = FALSE)
      taa2 <- translate(s2c(cdns2[x]), frame = 0, sens = "F", numcode = 5, NAstring = "X", ambiguous = FALSE)
      species <- anc_desSeqs$Species[y]
      # print(c(species, x))
      if (x > 1 & x < length(codonposition))
      {
        final <- rbind(final, c(species, x, cdns1[x-1], cdns1[x], cdns1[x+1], cdns2[x-1], cdns2[x], cdns2[x+1], taa1, taa2, brLength))
      }
      if (x == 1)
      {
        final <- rbind(final, c(species, x, '---', cdns1[x], cdns1[x+1], '---', cdns2[x], cdns2[x+1], taa1, taa2, brLength))
      }
      if (x == length(codonposition))
      {
        final <- rbind(final, c(species, x, cdns1[x-1], cdns1[x], '---', cdns2[x-1], cdns2[x], '---', taa1, taa2, brLength))
      }
    }
  }
  da <- data.frame(final)
  colnames(da) <- c("Species", "CodonPosition", "PreviousAncCodon", "AncestorCodon", "NextAncCodon", "PreviousDesCodon", "DescendantCodon", "NextDesCodon", "AncestralAA", "DescendantAA", "BranchLength")
  codonTable = rbind(codonTable, da)
}

table(codonTable$Species)

setdiff(anc_desSeqs$Species, as.character(codonTable$Species)) # nothing !

write.table(codonTable, '../results/cockroaches11_19/PolarizeMutations.CodonsTable.Cockroaches.txt', sep = '\t', quote = FALSE, row.names = FALSE)
