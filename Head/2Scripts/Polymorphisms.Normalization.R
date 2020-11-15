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

unzip("../../Body/2Derived/POLARIZEDBR_DATA_ML.zip", exdir = '../../Body/2Derived/')
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
Final = read.table("../../Body/3Results/Polymorphisms.Normalization.NeutralATGC.txt",
                   header = TRUE, sep = '\t')

Final = merge(Final,Taxa)
str(Final)
Final[,3:6] <- sapply(Final[,3:6], function(x) as.numeric(as.character(x)))
Agg1 = aggregate(list(Final$A,Final$T,Final$G,Final$C), by = list(Final$Species, Final$Family), FUN = sum); names(Agg1) = c('Species','Family','A','T','G','C')

Agg1$FrA = Agg1$A/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)
Agg1$FrT = Agg1$T/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)
Agg1$FrG = Agg1$G/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)
Agg1$FrC = Agg1$C/(Agg1$A+Agg1$T+Agg1$G+Agg1$C)

Agg2 = aggregate(list(Agg1$FrA,Agg1$FrT,Agg1$FrG,Agg1$FrC), by = list(Agg1$Family), FUN = median); names(Agg2) = c('Family','FrA','FrT','FrG','FrC')
Agg2$Sum = Agg2$FrA + Agg2$FrT + Agg2$FrG + Agg2$FrC # 1 

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
library(cowplot)

NormFrac = read.table("../../Body/3Results/Polymorphisms.Normalization.Normalized12Fractions.txt", header = TRUE)


pdf("../../Body/4Figures/Polymorphisms.Normalization.R.03.Bars.pdf", width = 15, height = 10)

a = ggbarplot(NormFrac, 'Subs', 'Normalised2Number', xlab="Substitution types", 
              ylab = 'Normalized frequencies', title = 'all',
              fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE) +
  scale_x_discrete(labels = sort(unique(NormFrac$Subs), decreasing = TRUE))

b = ggbarplot(NormFrac[NormFrac$Family == 'Termitidae',], 'Subs', 'Normalised2Number', 
              xlab="Substitution types", ylab = 'Normalized frequencies', title = 'Termitidae',
              fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE) +
  scale_x_discrete(labels = sort(unique(NormFrac$Subs), decreasing = TRUE))

rhino = NormFrac[NormFrac$Family == 'Rhinotermitidae',]
c = ggbarplot(rhino, 'Subs', 'Normalised2Number', xlab="Substitution types", 
              ylab = 'Normalized frequencies', title = 'Rhinotermitidae',
              fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", 
                                                         "#bdbdbd", "#feb24c", "#f03b20", 
                                                         "#bdbdbd", "#bdbdbd", "#bdbdbd",
                                                         "#2c7fb8", "#bdbdbd"), 
              combine = TRUE) +
  scale_x_discrete(labels = c('T_G', 'T_C', 'T_A', 'G_T', 'G_A', 'C_T', 'C_G', 
                              'C_A', 'A_T', 'A_G', 'A_C'))

lower = NormFrac[NormFrac$Family %in% c('Rhinotermitidae', 'Hodotermitidae', 'Termopsidae'),]

all_lower = ggbarplot(lower, 
                      'Subs', 'Normalised2Number', xlab="Substitution types", title = 'Lower termites',
fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE,
add = 'mean') + ylim(c(0, max(lower$Number + 0.05))) +
  scale_x_discrete(labels = c('T_G', 'T_C', 'T_A', 'G_T', 'G_A', 'C_T', 'C_G', 
                              'C_A', 'A_T', 'A_G', 'A_C'))

higher = ggbarplot(NormFrac[NormFrac$Family == 'Termitidae',], 'Subs', 'Normalised2Number', xlab="Substitution types", title = 'Higher termites',
                   fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", 
                                                              "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", 
                                                              "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE) +
ylim(c(0, max(lower$Number + 0.05))) + 
  scale_x_discrete(labels = sort(unique(NormFrac$Subs), decreasing = TRUE))


lessSocialPlot = ggbarplot(NormFrac[NormFrac$Family %in% lessSocial,], 
                       'Subs', 'Normalised2Number', xlab="Substitution types", title = 'Less social termites',
                       fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#2c7fb8"), combine = TRUE,
                       add = 'mean') + ylim(c(0, 0.45)) +
  scale_x_discrete(labels = c('T_G', 'T_C', 'T_A', 'G_T', 'G_A', 'C_T', 'A_T', 'A_G'))

moreSocialPlot = ggbarplot(NormFrac[NormFrac$Family %in% moreSocial,], 
                       'Subs', 'Normalised2Number', xlab="Substitution types", title = 'More social termites',
                       fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", 
                                                                  "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", 
                                                                  "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE,
                       add = 'mean') + ylim(c(0, 0.45)) +
  scale_x_discrete(labels = sort(unique(NormFrac$Subs), decreasing = TRUE))


ggarrange(a,                                                 # First row with scatter plot
          ggarrange(b, c, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A"                                        # Labels of the scatter plot
) 

ggarrange(all_lower, higher,
          ncol = 2, nrow = 2)

dev.off()

plots = plot_grid(lessSocialPlot, moreSocialPlot, nrow=1)

save_plot('~/Alina/other_projects/TermitesMutSpectrum/Body/4Figures/Polymorphisms.MoreLessSocial.pdf', 
          plots, base_height = 7)

### compare An>Gn in higher and lower termites

library(dplyr)

mut = read.table('../../Body/3Results/Polymorphisms.MutSpecData.txt', header = TRUE)

# External Internal 
# 662     1013

# mann-whitney test for normilized mutspec

normMut = mut %>%
  group_by(Species) %>% 
  mutate(
    AllSubs = length(Subs)
  ) %>%
  filter(Subs == 'T_C') %>% 
  mutate(
    T_C = length(Subs) / AllSubs
  ) %>% 
  distinct(Species, .keep_all = TRUE) %>%
  mutate(
    T_fr = T / (A + T + G + C),
    T_C_norm = T_C / T_fr
  )

t.test(normMut[normMut$Family == 'Termitidae',]$T_C_norm, normMut[normMut$Family != 'Termitidae',]$T_C_norm)

# 20 vs 

# t = 3.664, df = 30.091, p-value = 0.0009496
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.9909639 3.4861512
# sample estimates:
#   mean of x mean of y 
# 4.973166  2.734609


wilcox.test(normMut[normMut$Family == 'Termitidae',]$T_C_norm, normMut[normMut$Family != 'Termitidae',]$T_C_norm,
            conf.int = TRUE)

# W = 235, p-value = 0.00387
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   0.7241917 3.4648685
# sample estimates:
#   difference in location 
# 1.855148 

# fisher test

a = nrow(mut[mut$Subs == 'T_C' & mut$Family != 'Termitidae',]) # 250
b = nrow(mut[mut$Subs == 'T_C' & mut$Family == 'Termitidae',]) # 267
c = nrow(mut[mut$Subs != 'T_C' & mut$Family != 'Termitidae',]) # 659
d = nrow(mut[mut$Subs != 'T_C' & mut$Family == 'Termitidae',]) # 499

fisher.test(matrix(c(a, b, c, d), ncol = 2))

# p-value = 0.001208
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5726741 0.8778123
# sample estimates:
#   odds ratio 
# 0.7091567

###################################
### A>G in more and less social

unique(normMut$Family)

# why only 35 species ?

# more social versus less social

moreSocial = c('Termitidae', 'Rhinotermitidae', 'Mastotermitidae')
lessSocial = c('Hodotermitidae', 'Termopsidae')

t.test(normMut[normMut$Family %in% moreSocial,]$T_C_norm,
       normMut[normMut$Family %in% lessSocial,]$T_C_norm)

# t = 0.51448, df = 2.8361, p-value = 0.6443
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -2.790630  3.825537
# sample estimates:
#   mean of x mean of y 
# 4.058137  3.540684

# fisher

a = nrow(mut[mut$Subs == 'T_C' & mut$Family %in% moreSocial,]) # 439
b = nrow(mut[mut$Subs == 'T_C' & mut$Family %in% lessSocial,]) # 78
c = nrow(mut[mut$Subs != 'T_C' & mut$Family %in% moreSocial,]) # 1022
d = nrow(mut[mut$Subs != 'T_C' & mut$Family %in% lessSocial,]) # 136

fisher.test(matrix(c(a, b, c, d), ncol = 2))

# p-value = 0.06816
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5498975 1.0255416
# sample estimates:
#   odds ratio 
# 0.749069 

### all transitions in more and less social

# normMutAll = mut %>%
#   group_by(Species) %>% 
#   mutate(
#     AllSubs = length(Subs)
#   ) %>%
#   mutate(
#     T_C = length(Subs[Subs == 'T_C']) / AllSubs,
#     A_G = length(Subs[Subs == 'A_G']) / AllSubs,
#     C_T = length(Subs[Subs == 'C_T']) / AllSubs,
#     G_A = length(Subs[Subs == 'G_A']) / AllSubs
#   ) %>% 
#   distinct(Species, .keep_all = TRUE) %>%
#   mutate(
#     T_fr = T / (A + T + G + C),
#     A_fr = A / (A + T + G + C),
#     C_fr = C / (A + T + G + C),
#     G_fr = G / (A + T + G + C),
#     T_C_norm = T_C / T_fr,
#     A_G_norm = A_G / A_fr,
#     C_T_norm = C_T / C_fr,
#     G_A_norm = G_A / G_fr
#   )


a = nrow(mut[mut$Subs == 'G_A' & mut$Family %in% moreSocial,]) # 329
b = nrow(mut[mut$Subs == 'G_A' & mut$Family %in% lessSocial,]) # 25
c = nrow(mut[mut$Subs != 'G_A' & mut$Family %in% moreSocial,]) # 1132
d = nrow(mut[mut$Subs != 'G_A' & mut$Family %in% lessSocial,]) # 189

fisher.test(matrix(c(a, b, c, d), ncol = 2))

# p-value = 0.0001573
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.412716 3.543389
# sample estimates:
#   odds ratio 
# 2.196349 

a = nrow(mut[mut$Subs == 'C_T' & mut$Family %in% moreSocial,]) # 411
b = nrow(mut[mut$Subs == 'C_T' & mut$Family %in% lessSocial,]) # 76
c = nrow(mut[mut$Subs != 'C_T' & mut$Family %in% moreSocial,]) # 1050
d = nrow(mut[mut$Subs != 'C_T' & mut$Family %in% lessSocial,]) # 138

fisher.test(matrix(c(a, b, c, d), ncol = 2))

# p-value = 0.02949
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.5206664 0.9759042
# sample estimates:
#   odds ratio 
# 0.7108976 

a = nrow(mut[mut$Subs == 'A_G' & mut$Family %in% moreSocial,]) # 177
b = nrow(mut[mut$Subs == 'A_G' & mut$Family %in% lessSocial,]) # 18
c = nrow(mut[mut$Subs != 'A_G' & mut$Family %in% moreSocial,]) # 1284
d = nrow(mut[mut$Subs != 'A_G' & mut$Family %in% lessSocial,]) # 196

fisher.test(matrix(c(a, b, c, d), ncol = 2))

# p-value = 0.1372
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.8968913 2.6513302
# sample estimates:
#   odds ratio 
# 1.500698 

##################################################################
### fisher with external branches

mutExt = mut[mut$Branch == 'External',]

a = nrow(mutExt[mutExt$Subs == 'T_C' & mutExt$Family %in% moreSocial,]) # 177
b = nrow(mutExt[mutExt$Subs == 'T_C' & mutExt$Family %in% lessSocial,]) # 37
c = nrow(mutExt[mutExt$Subs != 'T_C' & mutExt$Family %in% moreSocial,]) # 394
d = nrow(mutExt[mutExt$Subs != 'T_C' & mutExt$Family %in% lessSocial,]) # 57

fisher.test(matrix(c(a, b, c, d), ncol = 2))

# p-value = 0.09572
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.424520 1.101423
# sample estimates:
#   odds ratio 
# 0.6807404 

a = nrow(mutExt[mutExt$Subs == 'G_A' & mutExt$Family %in% moreSocial,]) # 147
b = nrow(mutExt[mutExt$Subs == 'G_A' & mutExt$Family %in% lessSocial,]) # 13
c = nrow(mutExt[mutExt$Subs != 'G_A' & mutExt$Family %in% moreSocial,]) # 421
d = nrow(mutExt[mutExt$Subs != 'G_A' & mutExt$Family %in% lessSocial,]) # 81

fisher.test(matrix(c(a, b, c, d), ncol = 2))

# p-value = 0.01294
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.158943 4.385800
# sample estimates:
#   odds ratio 
# 2.173407 

a = nrow(mutExt[mutExt$Subs == 'C_T' & mutExt$Family %in% moreSocial,]) # 161
b = nrow(mutExt[mutExt$Subs == 'C_T' & mutExt$Family %in% lessSocial,]) # 35
c = nrow(mutExt[mutExt$Subs != 'C_T' & mutExt$Family %in% moreSocial,]) # 407
d = nrow(mutExt[mutExt$Subs != 'C_T' & mutExt$Family %in% lessSocial,]) # 59

fisher.test(matrix(c(a, b, c, d), ncol = 2))

# p-value = 0.0881
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.4138786 1.0875116
# sample estimates:
#   odds ratio 
# 0.6672598 

a = nrow(mutExt[mutExt$Subs == 'A_G' & mutExt$Family %in% moreSocial,]) # 46
b = nrow(mutExt[mutExt$Subs == 'A_G' & mutExt$Family %in% lessSocial,]) # 1
c = nrow(mutExt[mutExt$Subs != 'A_G' & mutExt$Family %in% moreSocial,]) # 522
d = nrow(mutExt[mutExt$Subs != 'A_G' & mutExt$Family %in% lessSocial,]) # 93

fisher.test(matrix(c(a, b, c, d), ncol = 2))

# p-value = 0.00854
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.360119 333.723434
# sample estimates:
#   odds ratio 
# 8.182109 


