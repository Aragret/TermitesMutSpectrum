rm(list=ls(all=TRUE))

library(gdata)

mut = read.table('../results/cockroaches11_19/4foldSubs.Cockroaches.txt', header = TRUE, sep='\t')

families = read.xls('../results/4foldSubsBrLen02Fractions.xlsx')

mutBrLen02 = mut[mut$BranchLength >= 0.02,]

### get table one row - one species

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
mutTable[mutTable$Species == 'AUS49_Tumulitermes_sp._1',]

brLen = aggregate(mut$BranchLength, by=list(mut$Species), unique)
names(brLen) = c('Species', 'BranchLength')

mutTableBrlen = merge(mutTable, brLen)

###### subs fractions

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

cor.test(mutFr$BranchLength, mutFr$sumOfSubs)
plot(mutFr$BranchLength, mutFr$sumOfSubs)

cor.test(mutFr$BranchLength, mutFr$C_T_fr)
plot(mutFr$BranchLength, mutFr$C_T_fr)

summary(mutFr$BranchLength)

mutFrBrlen02 = mutFr[mutFr$BranchLength >= 0.02,]

pdf('../results/cockroaches11_19/CorrBrLenFractions.pdf')

cor.test(mutFrBrlen02$BranchLength, mutFrBrlen02$sumOfSubs)
plot(mutFrBrlen02$BranchLength, mutFrBrlen02$sumOfSubs)

cor.test(mutFrBrlen02$BranchLength, mutFrBrlen02$C_T_fr)
plot(mutFrBrlen02$BranchLength, mutFrBrlen02$C_T_fr)
plot(mutFrBrlen02$BranchLength, mutFrBrlen02$A_G_fr)

summary(mutFrBrlen02$BranchLength)

for(i in 16:27){
  print(names(mutFrBrlen02)[i])
  print(cor.test(mutFrBrlen02[, i], mutFrBrlen02$BranchLength))
}


################
# check different branch length

brLenLess1 = mutFrBrlen02[mutFrBrlen02$BranchLength <= 0.1,]
brLenMore1 = mutFrBrlen02[mutFrBrlen02$BranchLength > 0.1,]

plot(brLenLess1$C_T_fr, brLenLess1$BranchLength, main = 'brlen less than 0.1')
plot(brLenLess1$A_G_fr, brLenLess1$BranchLength, main = 'brlen less than 0.1')
cor.test(brLenLess1$C_T_fr, brLenLess1$BranchLength, method = 'spearman')

plot(brLenMore1$C_T_fr, brLenMore1$BranchLength, main = 'brlen more than 0.1')
plot(brLenMore1$A_G_fr, brLenMore1$BranchLength, main = 'brlen more than 0.1')
cor.test(brLenMore1$C_T_fr, brLenMore1$BranchLength, method = 'spearman')

for(i in 16:27){
  print(names(brLenLess1)[i])
  print(cor.test(brLenLess1[, i], brLenLess1$BranchLength, method = 'spearman'))
}

for(i in 16:27){
  print(names(brLenMore1)[i])
  print(cor.test(brLenMore1[, i], brLenMore1$BranchLength, method = 'spearman'))
}

dev.off()

####### add meta data

data = merge(mutFrBrlen02, families[, c('Species', 'Taxonomy', 'Worker', 'diet.wood.soil')], by='Species', all.x = TRUE)

write.table(data, '../results/cockroaches11_19/4foldSubsBrLen02Fractions.Cockroaches.txt', sep='\t', quote = FALSE, row.names = FALSE)
