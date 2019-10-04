rm(list=ls(all=TRUE))

Mut = read.table('../../Body/3Results/Polymorphisms.MutSpecData.txt', header = TRUE)

## get ancestral nucleotide 
Mut$Subs = as.character(Mut$Subs)
TRIM<-function(x)	 {unlist(strsplit(x,'_'))[1]}
Mut$AncestralNuc = apply(as.matrix(Mut$Subs), 1 , FUN = TRIM)

## normalize by freq of ancestral nucleotides
Mut$Number = 1
A = Mut[Mut$AncestralNuc == 'A',]; A$Number = A$A/(A$A+A$T+A$G+A$C)
T = Mut[Mut$AncestralNuc == 'T',]; T$Number = T$T/(T$A+T$T+T$G+T$C)
G = Mut[Mut$AncestralNuc == 'G',]; G$Number = G$G/(G$A+G$T+G$G+G$C)
C = Mut[Mut$AncestralNuc == 'C',]; C$Number = C$C/(C$A+C$T+C$G+C$C)

Mut = rbind(A,T); Mut = rbind(Mut,G); Mut = rbind(Mut,C);
Mut$Number = 1/Mut$Number
agg = aggregate(Mut$Number, by = list(Mut$Species,Mut$Subs), FUN = sum)
names(agg)=c('Species','Subs','Freq')

## make vector of 12 Subs for each species
Template = data.frame(unique(Mut$Subs)); names(Template) = c('Subs'); Template$Freq = 0;
VecOfSpecies = unique(agg$Species)
for (i in 1:length(VecOfSpecies))
{ # i = 2
  Temp = agg[agg$Species == VecOfSpecies[i],]
  Template$Species = VecOfSpecies[i]
  Temp = merge(Temp,Template, by = c('Subs','Species'), all = TRUE)
  Temp[is.na(Temp)] <- 0
  Temp$Freq = Temp$Freq.x + Temp$Freq.y
  ALL = sum(Temp$Freq)
  
  Line = data.frame(VecOfSpecies[i], Temp[Temp$Subs == 'A_T',]$Freq/ALL, Temp[Temp$Subs == 'A_G',]$Freq/ALL, Temp[Temp$Subs == 'A_C',]$Freq/ALL, Temp[Temp$Subs == 'T_A',]$Freq/ALL, Temp[Temp$Subs == 'T_G',]$Freq/ALL, Temp[Temp$Subs == 'T_C',]$Freq/ALL, Temp[Temp$Subs == 'G_A',]$Freq/ALL, Temp[Temp$Subs == 'G_T',]$Freq/ALL, Temp[Temp$Subs == 'G_C',]$Freq/ALL, Temp[Temp$Subs == 'C_A',]$Freq/ALL, Temp[Temp$Subs == 'C_T',]$Freq/ALL, Temp[Temp$Subs == 'C_G',]$Freq/ALL)
  names(Line)=c('Species','A_T','A_G','A_C','T_A','T_G','T_C','G_A','G_T','G_C','C_A','C_T','C_G')
  if (i == 1) {Final = Line}
  if (i >  1) {Final = rbind(Final,Line)}
}

write.table(Final, '../../Body/3Results/Polymorphisms.MutSpecData.SynNorm.OneSpOneRow.txt', quote = FALSE, row.names = FALSE)

#####################################################################################

MUT = read.table('../../Body/3Results/Polymorphisms.MutSpecData.SynNorm.OneSpOneRow.txt', header = TRUE)

taxon = read.csv('../../Body/2Derived/TermitesFamilies.csv')
MUT = merge(MUT, taxon[, 1:2])

MUT$TsTv = (MUT$T_C + MUT$C_T + MUT$G_A + MUT$A_G) / (MUT$T_A + MUT$A_T + MUT$G_C + MUT$C_G + MUT$G_T + MUT$T_G + MUT$C_A + MUT$A_C)
summary(MUT$TsTv)
MutTsTv = MUT[MUT$TsTv < Inf,]

summary(MutTsTv[MutTsTv$Family == 'Termitidae',]$TsTv)
summary(MutTsTv[MutTsTv$Family == 'Rhinotermitidae',]$TsTv)

# too low number of mutations

###########################################################################
pdf('../../Body/4Figures/Polymorphisms.PCA.pdf')

MATRIX = MUT[, c(2:13)]
row.names(MATRIX)=MUT$Species
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

biplot(PCA, choices=c(1,2), col = c('white','black'), cex = 0.8) #  biplot(princomp(USArrests),choices=c(1,3))

###########################################################################
### only with two families

MATRIX = MUT[MUT$Family == 'Termitidae' | MUT$Family == 'Rhinotermitidae', c(2:13)]
row.names(MATRIX)=MUT[MUT$Family == 'Termitidae' | MUT$Family == 'Rhinotermitidae', 'Species']
matrix = MATRIX

PCA = prcomp(matrix, center = TRUE, scale = TRUE) #FALSE) # I don't scale because we analyze the same units (fraction from MutSpec) 
print(PCA)  
summary(PCA)
MATRIX$Pca1 = PCA$x[,1]
MATRIX$Pca2 = PCA$x[,2]
MATRIX$Pca3 = PCA$x[,3]
MATRIX$Pca4 = PCA$x[,4]

MATRIX = cbind(MATRIX, MUT[MUT$Family == 'Termitidae' | MUT$Family == 'Rhinotermitidae', 'Family'])
names(MATRIX)[17] = 'Family'

biplot(PCA, choices=c(1,2), col = c('white','black'), cex = 0.8) #  biplot(princomp(USArrests),choices=c(1,3))

for(i in 1:nrow(MATRIX)){
  if(MATRIX[i, 'Family'] == 'Termitidae'){
    MATRIX$SIZE[i] = 0
  }
  if(MATRIX[i, 'Family'] == 'Rhinotermitidae'){
    MATRIX$SIZE[i] = 1
  }
}
MATRIX$SIZE 
summary(MATRIX$SIZE)

plot(MATRIX$Pca1,MATRIX$Pca2, pch = 16, cex = 1.5, col = rgb(0,MATRIX$SIZE,0,1), xlim = c(-11,3), ylim = c(-3,4.5)); 
biplot(PCA, choices=c(1,2), scale = 0.5, col = c('grey','black'), cex = 0.5) # , xlim = c(-11,3), ylim = c(-3,4.5))

dev.off()
