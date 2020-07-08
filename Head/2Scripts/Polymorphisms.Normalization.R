rm(list=ls(all=TRUE))

VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')
length(unique(VecOfSynFourFoldDegenerateSites)) # 32

unzip("../../Body/2Derived/POLARIZEDBR_DATA_ML.zip")
List = list.files("../../Body/2Derived/POLARIZEDBR_DATA_ML/")

for (i in 1:length(List)){
  # i = 1
  infile = paste("../../Body/2Derived/POLARIZEDBR_DATA_ML/",as.character(List[i]),sep='') 
  if (length(grep('POLARISED',infile)) > 0)
  {
    Species = gsub('\\..*','',as.character(List[i]))
    Gene = gsub(Species,'',as.character(List[i])); Gene = gsub('.POLARISED.txt','',Gene); Gene = gsub('\\.POLARISED.txt','',Gene); Gene = gsub('\\.','',Gene); 
    GeneSpecies = read.table(infile, header = TRUE)
    GeneSpecies = GeneSpecies[GeneSpecies$BranchPosition == 'External',]
    ExternalSeqsTogether = paste(GeneSpecies$MoreShallowNodeSeq,collapse = '')
    ExternalSeqsTogether = unlist(strsplit(ExternalSeqsTogether,'')) # 5700/3
    CodonsVec = c(); StartNuc = 1
    if (length(ExternalSeqsTogether)/3 == round(length(ExternalSeqsTogether)/3))  # if divide by 3 without the rest
    {
      for (j in 1:(length(ExternalSeqsTogether)/3))
      {
        CodonsVec = c(CodonsVec,paste(ExternalSeqsTogether[StartNuc : (StartNuc+2)],collapse = ''))
        StartNuc = StartNuc+3
      }
      AllCodons = length(CodonsVec)        # 1021
      CodonsVecNeutral = CodonsVec[CodonsVec %in% VecOfSynFourFoldDegenerateSites]
      NeutralCodons = length(CodonsVecNeutral) # 1900
      data.frame(table(CodonsVecNeutral))
      
      CodonsVecNeutral = gsub("CTA|GTA|TCA|CCA|ACA|GCA|CGA|GGA",'A',CodonsVecNeutral)
      CodonsVecNeutral = gsub("CTT|GTT|TCT|CCT|ACT|GCT|CGT|GGT",'T',CodonsVecNeutral)
      CodonsVecNeutral = gsub("CTG|GTG|TCG|CCG|ACG|GCG|CGG|GGG",'G',CodonsVecNeutral)
      CodonsVecNeutral = gsub("CTC|GTC|TCC|CCC|ACC|GCC|CGC|GGC",'C',CodonsVecNeutral)
      
      Line=c(Species,Gene,length(CodonsVecNeutral[CodonsVecNeutral == 'A']),length(CodonsVecNeutral[CodonsVecNeutral == 'T']),length(CodonsVecNeutral[CodonsVecNeutral == 'G']),length(CodonsVecNeutral[CodonsVecNeutral == 'C']), AllCodons, NeutralCodons)
      if (i == 1) {Final = Line}
      if (i >  1) {Final = rbind(Final,Line)}
    }
  }
}

Final = as.data.frame(Final); names(Final)=c('Species','Gene','A','T','G','C','NumberOfAllCodons',"NumberOfFourFoldDegenCodons")
write.table(Final, "../../Body/3Results/Polymorphisms.Normalization.NeutralATGC.txt", quote = FALSE, row.names = FALSE, sep='\t')

## delete all unziped files
files <- list.files("../../Body/2Derived/POLARIZEDBR_DATA_ML/")
for (i in 1:length(files)) 
{ # i = 1
  file = paste('../../Body/2Derived/POLARIZEDBR_DATA_ML/',files[i],sep='')
  if (file.exists(file)) file.remove(file)
}

################### 
Taxa = read.csv("../../Body/2Derived/TermitesFamilies.csv") 


Final = merge(Final,Taxa)
str(Final)
Final[,3:6] <- sapply(Final[,3:6], function(x) as.numeric(as.character(x)))
Agg1 = aggregate(list(Final$A,Final$T,Final$G,Final$C), by = list(Final$Species, Final$Family), FUN = sum); names(Agg1) = c('Species','Family','A','T','G','C')

Agg1$FrA = Agg1$A/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)
Agg1$FrT = Agg1$T/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)
Agg1$FrG = Agg1$G/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)
Agg1$FrC = Agg1$C/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)

Agg2 = aggregate(list(Agg1$FrA,Agg1$FrT,Agg1$FrG,Agg1$FrC), by = list(Agg1$Family), FUN = mean); names(Agg2) = c('Family','FrA','FrT','FrG','FrC')
Agg2$Sum = Agg2$FrA + Agg2$FrT + Agg2$FrG + Agg2$FrC # 1 
# Agg2 = Agg2[Agg2$Class %in% c('Actinopterygii','Amphibia','Reptilia','Mammalia','Aves'),]

row.names(Agg2)=Agg2$Family
BARPLOT = t(as.matrix(Agg2[,c(2:5)]))
BARPLOT = BARPLOT[,c(1,5,2,4,3)] # column order
BARPLOT = BARPLOT[c(3,2,4,1),] # row order

# Get the stacked barplot
ColG = rgb(0.1,0.1,0.1,0.5)
ColT = rgb(0.1,0.1,1,0.5)
ColC = rgb(0.1,1,0.1,0.5)
ColA = rgb(1,0.1,0.1,0.5)

pdf("../../Body/4Figures/Polymorphisms.Normalization.R.01.pdf", width = 30, height = 20)

barplot(BARPLOT, col = c(ColG,ColT,ColC,ColA), border="white", space=0.04, font.axis=2, xlab="")

PieChartTable = read.table("../../Body/3Results/Polymorphisms.BetweenFamiliesWithoutNormalization.PieChartTable.txt")
PieChartTable$AncestralNuc = NA
for (i in 1:nrow(PieChartTable)) {PieChartTable$AncestralNuc[i] = unlist(strsplit(as.character(PieChartTable$Subs[i]),split = '_'))[1]}

PieChartTable = merge(PieChartTable,Agg2[,c(1:5)])
PieChartTable$Normalised1Number = 0
for (i in 1:nrow(PieChartTable))
{
  if (PieChartTable$AncestralNuc[i] == 'A') {PieChartTable$Normalised1Number[i] = PieChartTable$Number[i] / PieChartTable$FrA[i]}
  if (PieChartTable$AncestralNuc[i] == 'T') {PieChartTable$Normalised1Number[i] = PieChartTable$Number[i] / PieChartTable$FrT[i]}
  if (PieChartTable$AncestralNuc[i] == 'G') {PieChartTable$Normalised1Number[i] = PieChartTable$Number[i] / PieChartTable$FrG[i]}
  if (PieChartTable$AncestralNuc[i] == 'C') {PieChartTable$Normalised1Number[i] = PieChartTable$Number[i] / PieChartTable$FrC[i]}
}

Agg = aggregate(PieChartTable$Normalised1Number, by = list(PieChartTable$Family), FUN = sum); names(Agg) = c('Family','Total')
PieChartTable = merge(PieChartTable,Agg)
PieChartTable$Normalised2Number = PieChartTable$Normalised1Number/PieChartTable$Total

write.table(PieChartTable,"../../Body/3Results/Polymorphisms.Normalization.Normalized12Fractions.txt", quote = FALSE, row.names = FALSE)

par(mfrow=c(1,5))
VecOfFamilies = as.character(unique(Agg$Family))
for (i in 1:length(VecOfFamilies))
{
  pie(PieChartTable[PieChartTable$Family == VecOfFamilies[i],]$Normalised2Number, labels = PieChartTable[PieChartTable$Family == VecOfFamilies[i],]$Subs, main = VecOfFamilies[i], col=rainbow(12))
}

dev.off()

##########barplots########
library("ggpubr")

NormFrac = read.table("../../Body/3Results/Polymorphisms.Normalization.Normalized12Fractions.txt", header = TRUE)


pdf("../../Body/4Figures/Polymorphisms.Normalization.R.03.Bars.pdf", width = 15, height = 10)

a = ggbarplot(NormFrac, 'Subs', 'Normalised2Number', xlab="Substitution types", title = 'all',
              fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE)

b = ggbarplot(NormFrac[NormFrac$Family == 'Termitidae',], 'Subs', 'Normalised2Number', xlab="Substitution types", title = 'Termitidae',
              fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE)

c = ggbarplot(NormFrac[NormFrac$Family == 'Rhinotermitidae',], 'Subs', 'Normalised2Number', xlab="Substitution types", title = 'Rhinotermitidae',
              fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE)

lower = NormFrac[NormFrac$Family %in% c('Rhinotermitidae', 'Hodotermitidae', 'Termopsidae'),]

all_lower = ggbarplot(lower, 
                      'Subs', 'Normalised2Number', xlab="Substitution types", title = 'Lower termites',
fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE,
add = 'mean') + ylim(c(0, max(lower$Number + 0.05)))

higher = ggbarplot(NormFrac[NormFrac$Family == 'Termitidae',], 'Subs', 'Normalised2Number', xlab="Substitution types", title = 'Higher termites',
                   fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE) +
ylim(c(0, max(lower$Number + 0.05)))


ggarrange(a,                                                 # First row with scatter plot
          ggarrange(b, c, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A"                                        # Labels of the scatter plot
) 

ggarrange(all_lower, higher,
          ncol = 2, nrow = 2)

dev.off()
