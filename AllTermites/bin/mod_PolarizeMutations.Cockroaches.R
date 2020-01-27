rm(list = ls(all=TRUE))

library(ape)
library(seqinr)
library(stringr)

tree = read.tree('../results/rootedTreeCockroaches.newick')

numberOfSpecies = length(tree$tip.label)

a = as.data.frame(tree$edge)
a = cbind(a, tree$edge.length)
externalBranches = a[a$V2 <= numberOfSpecies,]
externalBranches = cbind(externalBranches, tree$tip.label)

names(externalBranches) = c('Node', 'Tip', 'BranchLength', 'Species')

tempTree = drop.tip(tree, 'ISO845_Zootermopsis_angusticollis_1')
tempExtBranches = a[a$V2 <= numberOfSpecies - 1,]
tempExtBranches = cbind(tempExtBranches, tempTree$tip.label)

anc = read.table('../results/cockroaches11_19/phylogeny/iqtree_mito.state', header=TRUE)
# anc$Node = sub('Node', '', anc$Node)

agg = aggregate(State ~ Node, anc, paste, collapse = "")

table(anc$Node)

#############################################################################

# tree = read.tree('../results/rootedTreeCockroaches.newick')

tree = read.tree('../results/cockroaches11_19/phylogeny/iqtree_mito.treefile')

a = as.data.frame(tree$edge)
a = cbind(a, tree$edge.length)

names(a) = c('Node', 'Tip', 'BranchLength')

externalBranches = a[a$Tip <= numberOfSpecies,]
externalBranches = cbind(externalBranches, tree$tip.label)

names(externalBranches) = c('Node', 'Tip', 'BranchLength', 'Species')



library(geiger)

max_node_number = max(tree$edge)
min_node_number = length(tree$tip.label) + 1

one_line = c()
for (i in min_node_number:max_node_number){
  descendants = tips(tree, i)
  if (length(descendants) == 2){
    one_line = rbind(one_line, c(i, descendants))
  }
}

sisters = as.data.frame(one_line)
names(sisters) = c('NodeNumber', 'Species1', 'Species2')

sistersBranchL = merge(sisters, externalBranches, by.x = 'Species1', by.y = 'Species')
sistersBranchL = merge(sistersBranchL, externalBranches[, c('Species', 'BranchLength')],
                       by.x = 'Species2', by.y = 'Species')

toRemove = c()
for(i in 1:nrow(sistersBranchL)){
  if(sistersBranchL[i,]$BranchLength.x >= 0.02 & sistersBranchL[i,]$BranchLength.y >= 0.02){
    toRemove = c(toRemove, i)
  }
}

smallBranches = sistersBranchL[-toRemove,] #110 rows

tipsToRemove = c()
nodesNeeded = c()

for(i in 1:nrow(smallBranches)){
  # i = 1
  if(smallBranches[i,]$BranchLength.x >= 0.02){
    tipsToRemove = c(tipsToRemove, as.character(smallBranches[i,]$Species2))
    next
  }
  if(smallBranches[i,]$BranchLength.y >= 0.02){
    tipsToRemove = c(tipsToRemove, as.character(smallBranches[i,]$Species1))
    next
  }
  tipsToRemove = c(tipsToRemove, as.character(smallBranches[i,]$Species2))
  lastAncNode = smallBranches[i,]$Node
  previousAnc = a[a$Tip == lastAncNode,]$Node
  brLenPrevAnc = a[a$Tip == previousAnc,]$BranchLength
  while(brLenPrevAnc <= 0.02){
    previousAnc = a[a$Tip == previousAnc,]$Node
    if(previousAnc == 680){
      print('680 ancestor')
      break
    }
    brLenPrevAnc = brLenPrevAnc + a[a$Tip == previousAnc,]$BranchLength
  }
  nodesNeeded = rbind(nodesNeeded, c(as.character(smallBranches[i,]$Species1), 
                                     previousAnc, brLenPrevAnc))
}

nodesNeeded = as.data.frame(nodesNeeded)
names(nodesNeeded) = c('Species', 'Node', 'BranchLength')

prunedTree = drop.tip(tree, tipsToRemove)

externalBranches$Species = as.character(externalBranches$Species)
nodesNeeded$Species = as.character(nodesNeeded$Species)

modExternalBranches = externalBranches[!(externalBranches$Species %in% c(tipsToRemove, nodesNeeded$Species)),]
modExternalBranches = modExternalBranches[, -2]
modExternalBranches = rbind(modExternalBranches, nodesNeeded)

modExternalBranches$Node = as.integer(modExternalBranches$Node)
modExternalBranches$BranchLength = as.numeric(modExternalBranches$BranchLength)

summary(modExternalBranches$BranchLength)

a = merge(externalBranches, modExternalBranches)

############################################################################
numberOfSpecies = length(tree$tip.label)

agg$Node = as.numeric(sub('Node', '', agg$Node))
agg$Node = agg$Node + numberOfSpecies

data = merge(agg, modExternalBranches, by='Node')

### adding aligns

aligns = read.fasta('../results/cockroaches11_19/aligns/align_cat.fasta', as.string = TRUE)

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

summary(data$BranchLength)
nrow(data[data$BranchLength < 0.02,])

lessThan02 = data[data$BranchLength < 0.02,]

write.table(codonTable, '../results/cockroaches11_19/mod_PolarizeMutations.CodonsTable.Cockroaches.txt', sep = '\t', quote = FALSE, row.names = FALSE)

old_codonTable = read.table('../results/cockroaches11_19/PolarizeMutations.CodonsTable.Cockroaches.txt', header = TRUE, sep='\t')
