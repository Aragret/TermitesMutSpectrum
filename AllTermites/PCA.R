
mut = read.table('results/4foldSubsRatesBrLen02.txt', header=TRUE, sep='\t')

MATRIX = mut[, c(27:38)] # residuals
MATRIX = mut[, c(2:13)] # subs number
row.names(MATRIX)=mut$Species
matrix = MATRIX

PCA = prcomp(matrix, center = TRUE, scale = TRUE) #FALSE) # I don't scale because we analyze the same units (fraction from MutSpec) 
print(PCA)  
summary(PCA)
MATRIX$Pca1 = PCA$x[,1]
MATRIX$Pca2 = PCA$x[,2]
MATRIX$Pca3 = PCA$x[,3]
MATRIX$Pca4 = PCA$x[,4]

MATRIX = cbind(MATRIX, MUT[, 'Family'])
names(MATRIX)[17] = 'Family'

pdf('results/PCA.pdf')
biplot(PCA, choices=c(1,2), col = c('white','black'), cex = 0.8) #  biplot(princomp(USArrests),choices=c(1,3))
dev.off()
