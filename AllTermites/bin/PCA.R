
mut = read.table('results/4foldSubsRatesBrLen02Residuals.txt', header=TRUE, sep='\t')

a = matrix(0, ncol=12, nrow=nrow(mut))

mut$sumOfSubs = mut[, "C_T"] + mut[, "A_G"] + mut[, "T_C"] + mut[, "G_T"] + mut[, "A_C"] +
  mut[, "G_A"] + mut[, "A_T"] + mut[, "G_C"] + mut[, "C_A"] + mut[, "C_G"] + mut[, "T_A"] +
  mut[, "T_G"]

for(i in 1:12){
  a[, i] = mut[, i+1] / mut$sumOfSubs
}

fractions = as.data.frame(a)
names(fractions) = sub(' ', '', paste(names(mut[, 2:13]), '_fr'))

mutFr = cbind(mut, fractions)

write.table(mutFr[, c(1:14, 39:51)], 'results/4foldSubsBrLen02Fractions.txt', sep='\t',
            row.names = FALSE, quote = FALSE)

# MATRIX = mut[, c(27:38)] # residuals
# MATRIX = mut[, c(2:13)] # subs number
MATRIX = mutFr[, c(40:51)]

row.names(MATRIX)=mutFr$Species
matrix = MATRIX

PCA = prcomp(matrix, center = TRUE, scale = TRUE) #FALSE) # I don't scale because we analyze the same units (fraction from MutSpec) 
print(PCA)  
summary(PCA)
MATRIX$Pca1 = PCA$x[,1]
MATRIX$Pca2 = PCA$x[,2]
MATRIX$Pca3 = PCA$x[,3]
MATRIX$Pca4 = PCA$x[,4]

pdf('results/PCA.pdf')
biplot(PCA, choices=c(1,2), col = c('white','black'), cex = 0.8) #  biplot(princomp(USArrests),choices=c(1,3))
dev.off()
