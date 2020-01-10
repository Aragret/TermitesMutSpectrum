rm(list = ls(all=TRUE))

library(gdata)
library(ggplot2)
library(gridExtra)

mut4f = read.table('../results/cockroaches11_19/mod_4foldSubs.Cockroaches.txt', header=TRUE, sep='\t')
mutAll = read.table('../results/cockroaches11_19/mod_AllSubs.Cockroaches.txt')

families = read.xls('../results/4foldSubsBrLen02Fractions.xlsx')

mutBrLen02 = mut4f[mut4f$BranchLength >= 0.02,]

mutBrLen02$Subs = as.character(mutBrLen02$Subs)

speciesNumber = length(unique(mutBrLen02$Species))
mutTable = setNames(data.frame(matrix(0, ncol = 12, nrow = speciesNumber)), unique(mutBrLen02$Subs))
Species = unique(mutBrLen02$Species)
mutTable = cbind(Species, mutTable)

for(i in 1:nrow(mutBrLen02)){
  row = mutBrLen02[i, ]
  mutTable[mutTable$Species == row$Species, row$Subs] = mutTable[mutTable$Species == row$Species, row$Subs] + 1
}

table(mutBrLen02[mutBrLen02$Species == "AUS49_Tumulitermes_sp._1",]$Subs)
mutTable[mutTable$Species == 'AUS49_Tumulitermes_sp._1',]

brLen = aggregate(mutBrLen02$BranchLength, by=list(mutBrLen02$Species), unique)
names(brLen) = c('Species', 'BranchLength')

mutTableBrlen = merge(mutTable, brLen)

#############################################################################
# get fractions

a = matrix(0, ncol=12, nrow=nrow(mutTableBrlen))

mutTableBrlen$sumOfSubs = mutTableBrlen[, "C_T"] + mutTableBrlen[, "A_G"] + mutTableBrlen[, "T_C"] + mutTableBrlen[, "G_T"] + mutTableBrlen[, "A_C"] +
  mutTableBrlen[, "G_A"] + mutTableBrlen[, "A_T"] + mutTableBrlen[, "G_C"] + mutTableBrlen[, "C_A"] + mutTableBrlen[, "C_G"] + mutTableBrlen[, "T_A"] +
  mutTableBrlen[, "T_G"]

for(i in 1:12){
  a[, i] = mutTableBrlen[, i+1] / mutTableBrlen$sumOfSubs
}

fractions = as.data.frame(a)
names(fractions) = sub(' ', '', paste(names(mutTableBrlen[, 2:13]), '_fr'))

mutFr = cbind(mutTableBrlen, fractions)

summary(mutFr$A_T_fr)

pdf('../results/cockroaches11_19/histOfMutFreqs.pdf')
par(mfrow=c(3,4))
for(i in 16:27){
  # i = 16
  hist(mutFr[, i], main = names(mutFr)[i], xlim = c(0, 0.6))
}

dev.off()
