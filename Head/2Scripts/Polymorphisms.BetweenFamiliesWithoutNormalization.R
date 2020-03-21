rm(list=ls(all=TRUE))

library("ggpubr")

pdf("../../Body/4Figures/Polymorphisms.BetweenFamiliesWithoutNormalization.R.01.pdf", width = 30, height = 20)

ColG = rgb(0.1,0.1,0.1,0.5)
ColT = rgb(0.1,0.1,1,0.5)
ColC = rgb(0.1,1,0.1,0.5)
ColA = rgb(1,0.1,0.1,0.5)

MUT = read.table("../../Body/3Results/Mutational_spectra_in_Termites.txt", header = TRUE)
length(unique(MUT$Species)) # 37

Taxa = read.csv("../../Body/2Derived/TermitesFamilies.csv") 

#########################
VecOfNormalSubstitutions <- c('A_C','C_A',
                              'A_G','G_A',
                              'C_G','G_C',
                              'C_T','T_C',
                              'G_T','T_G',
                              'T_A','A_T')
nrow(MUT)

SP = data.frame(table(MUT$Species)); names(SP) = c('Species','NumberOfAllSubst')
SPN = data.frame(table(MUT[MUT$Subs %in% VecOfNormalSubstitutions,]$Species)); names(SPN) = c('Species','NumberOfNormalSubst')
SP = merge(SP,SPN); SP$FractionOfNormal = SP$NumberOfNormalSubst/SP$NumberOfAllSubst
hist(SP$FractionOfNormal)
summary(SP$FractionOfNormal)
SpeciesToDelete = as.character(SP[SP$FractionOfNormal <=0.95,]$Species); length(SpeciesToDelete) # 1
# MUT = MUT[!MUT$Species %in% SpeciesToDelete,]
# MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]

##### FILTER 2: Synonymous Substitutions
nrow(MUT) # 1884
MUT = MUT[as.character(MUT$AncestralAA) == as.character(MUT$DescendantAA),]; nrow(MUT) # 1675

##### FILTER 3: fourfold degenerate sites:
VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG')
length(unique(VecOfSynFourFoldDegenerateSites)) # 32
nrow(MUT) # 1675
MUT4f = MUT[MUT$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & MUT$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]; nrow(MUT4f) # 209120
# 611 subs

MUT = merge(MUT, Taxa[, c('Species', 'Family')], all.x = TRUE)  ##### NOT ALL SPECIES HAVE TAXONOMY!!!!


PieChartTable = c()
Equil = c()
# par(mfrow = c(1,5))
VecOfFamilies = as.character(unique(MUT$Family))
for (i in 1:length(VecOfFamilies))
{ # i = 1
  Temp = MUT[MUT$Family == VecOfFamilies[i],]
  NumberOfSpecies = length(unique(Temp$Species))
  title = paste(VecOfFamilies[i], ', N = ', NumberOfSpecies, sep = '')
  Temp$number = 1
  Agg = aggregate(Temp$number, by = list(Temp$Subs), FUN = sum); names(Agg) = c('Subs','Number')
  Agg$Number = Agg$Number/sum(Agg$Number)
  sum(Agg$Number) # 1 - 100%
  # pie(Agg$Number, labels = Agg$Subs, main = title, col=rainbow(length(Agg$Subs)))
  
  Agg$Family = VecOfFamilies[i];
  PieChartTable = rbind(PieChartTable,Agg)
  
  ToGFromG = sum(Agg[Agg$Subs %in% c('A_G','T_G','C_G'),]$Number) /  sum(Agg[Agg$Subs %in% c('G_A','G_T','G_C'),]$Number)
  ToTFromT = sum(Agg[Agg$Subs %in% c('A_T','G_T','C_T'),]$Number) /  sum(Agg[Agg$Subs %in% c('T_A','T_G','T_C'),]$Number)
  ToAFromA = sum(Agg[Agg$Subs %in% c('T_A','G_A','C_A'),]$Number) /  sum(Agg[Agg$Subs %in% c('A_T','A_G','A_C'),]$Number)
  ToCFromC =  sum(Agg[Agg$Subs %in% c('T_C','G_C','A_C'),]$Number) / sum(Agg[Agg$Subs %in% c('C_T','C_G','C_A'),]$Number)
  
  Equil = rbind(Equil,c(i,VecOfFamilies[i],ToGFromG,ToTFromT,ToAFromA,ToCFromC))
}

PieChartTable = data.frame(PieChartTable)
write.table(PieChartTable,"../../Body/3Results/Polymorphisms.BetweenFamiliesWithoutNormalization.PieChartTable.txt")
Equil = data.frame(Equil); names(Equil) = c('i','Classes','G','T','C','A')

###################################

a = ggbarplot(PieChartTable, 'Subs', 'Number', xlab="Substitution types", title = 'all',
          fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE)

b = ggbarplot(PieChartTable[PieChartTable$Family == 'Termitidae',], 'Subs', 'Number', xlab="Substitution types", title = 'Termitidae',
          fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE)

c = ggbarplot(PieChartTable[PieChartTable$Family == 'Rhinotermitidae',], 'Subs', 'Number', xlab="Substitution types", title = 'Rhinotermitidae',
          fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE)

ggarrange(a,                                                 # First row with scatter plot
          ggarrange(b, c, ncol = 2, labels = c("B", "C")), # Second row with box and dot plots
          nrow = 2, 
          labels = "A"                                        # Labels of the scatter plot
) 

# ggbarplot(PieChartTable[PieChartTable$Family == 'Hodotermitidae',], 'Subs', 'Number', xlab="Substitution types", title = 'Hodotermitidae',
#           fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"))

# ggbarplot(PieChartTable[PieChartTable$Family == 'Termopsidae',], 'Subs', 'Number', xlab="Substitution types", title = 'Termopsidae',
#           fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"))


### equilibrium for each species:
Equilibrium = c()
MUT$number = 1
AGG = aggregate(MUT$number, by = list(MUT$Subs,MUT$Species,MUT$Family), FUN = sum); names(AGG) = c('Subs','Species','Family','Number')
VecOfSpecies = unique(MUT$Species)
for (i in 1:length(VecOfSpecies))
{ # i = 1
  Agg = AGG[AGG$Species == VecOfSpecies[i],]
  if (sum(Agg$Number) >= 15)
  {
    ToGFromG = log2(sum(Agg[Agg$Subs %in% c('A_G','T_G','C_G'),]$Number) /  sum(Agg[Agg$Subs %in% c('G_A','G_T','G_C'),]$Number))
    ToTFromT = log2(sum(Agg[Agg$Subs %in% c('A_T','G_T','C_T'),]$Number) /  sum(Agg[Agg$Subs %in% c('T_A','T_G','T_C'),]$Number))
    ToAFromA = log2(sum(Agg[Agg$Subs %in% c('T_A','G_A','C_A'),]$Number) /  sum(Agg[Agg$Subs %in% c('A_T','A_G','A_C'),]$Number))
    ToCFromC =  log2(sum(Agg[Agg$Subs %in% c('T_C','G_C','A_C'),]$Number) / sum(Agg[Agg$Subs %in% c('C_T','C_G','C_A'),]$Number))
    
    Equilibrium = rbind(Equilibrium, c(as.character(VecOfSpecies[i]), as.character(Agg$Family[1]), ToGFromG, ToTFromT, ToAFromA, ToCFromC))
  }
}
Equilibrium = data.frame(Equilibrium); names(Equilibrium)=c('Species','Family','G','T','A','C')
Equilibrium[,3:6] <- sapply(Equilibrium[,3:6], function(x) as.numeric(as.character(x)))
write.table(Equilibrium,"../../Body/3Results/Polymorphisms.BetweenFamiliesWithoutNormalization.EquilibriumLog2ToFrom.txt", quote = FALSE, row.names = FALSE)

par(mfrow=c(1,1))
Equilibrium = Equilibrium[Equilibrium$Family %in% VecOfFamilies,]
boxplot(
  Equilibrium[Equilibrium$Family == 'Termitidae',]$G, Equilibrium[Equilibrium$Family == 'Termitidae',]$T, Equilibrium[Equilibrium$Family == 'Termitidae',]$C, Equilibrium[Equilibrium$Family == 'Termitidae',]$A, 
  Equilibrium[Equilibrium$Family == 'Rhinotermitidae',]$G, Equilibrium[Equilibrium$Family == 'Rhinotermitidae',]$T, Equilibrium[Equilibrium$Family == 'Rhinotermitidae',]$C, Equilibrium[Equilibrium$Family == 'Rhinotermitidae',]$A, 
  Equilibrium[Equilibrium$Family == 'Hodotermitidae',]$G, Equilibrium[Equilibrium$Family == 'Hodotermitidae',]$T, Equilibrium[Equilibrium$Family == 'Hodotermitidae',]$C, Equilibrium[Equilibrium$Family == 'Hodotermitidae',]$A, 
  Equilibrium[Equilibrium$Family == 'Termopsidae',]$G, Equilibrium[Equilibrium$Family == 'Termopsidae',]$T, Equilibrium[Equilibrium$Family == 'Termopsidae',]$C, Equilibrium[Equilibrium$Family == 'Termopsidae',]$A, 
  at = c(1,2,3,4, 6,7,8,9, 11,12,13,14, 16,17,18,19),
  outline = FALSE, notch = FALSE, col = c(ColG,ColT,ColC,ColA), names = rep(c('G','T','C','A'),4))
title(xlab = 'Termitidae, Rhinotermitidae, Hodotermitidae, Termopsidae')
abline(h = 0, col = 'red')

dev.off()

#################################################################################
### Ts/Tv 

MUT4f = merge(MUT4f, Taxa[, c('Species', 'Family')], all.x = TRUE)
MUT4f$Ts = replicate(nrow(MUT4f), 0)
MUT4f$Tv = replicate(nrow(MUT4f), 0)

for(i in 1:nrow(MUT4f)){
  # i = 1
  if(MUT4f[i, 'Subs'] %in% c('A_G', 'G_A', 'T_C', 'C_T')){
    MUT4f[i, 'Ts'] = 1
  }
  else(MUT4f[i, 'Tv'] = 1)
}

agg = aggregate(MUT4f[, c('Ts', 'Tv')], by=list(MUT4f$Family), sum)
names(agg) = c('Family', 'Ts', 'Tv')

