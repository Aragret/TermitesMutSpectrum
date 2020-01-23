rm(list=ls(all=TRUE))

library(ggplot2)
library(ggpubr)

codonTable = read.table('../results/nd6_22_01//mod_PolarizeMutations.CodonsTable.nd6.txt', sep='\t', header=TRUE)

withoutGapsCodonTable = codonTable[as.character(codonTable$DescendantCodon) != '---',]
# 982204 subs

length(unique(codonTable$Species)) # 571
length(unique(withoutGapsCodonTable$Species)) # 571

# setdiff(codonTable$Species, withoutGapsCodonTable$Species)

FirstCodon = withoutGapsCodonTable$AncestorCodon
SecondCodon = withoutGapsCodonTable$DescendantCodon

data = data.frame(FirstCodon, SecondCodon)
data$FirstSecond = paste(data$FirstCodon, data$SecondCodon,sep = '_')

COMPAR<-function(x)	{
  cod1 <- unlist(strsplit(x,'_'))[1]
  cod2 <- unlist(strsplit(x,'_'))[2]
  #table <- c("FirstCOD", "SecondCOD")
  NumSub <- 0
  unlist(strsplit(cod1,''))[1]
  if (unlist(strsplit(cod1,''))[1] != unlist(strsplit(cod2,''))[1]){
    NumSub = 1 + NumSub
    FirstC <- unlist(strsplit(cod1,''))[1]
    SecC <- unlist(strsplit(cod2,''))[1]
  } 
  if (unlist(strsplit(cod1,''))[2] != unlist(strsplit(cod2,''))[2]){
    NumSub = 1 + NumSub
    FirstC <- unlist(strsplit(cod1,''))[2]
    SecC <- unlist(strsplit(cod2,''))[2]
  }
  if (unlist(strsplit(cod1,''))[3] != unlist(strsplit(cod2,''))[3])
  {
    NumSub = 1 + NumSub
    FirstC <- unlist(strsplit(cod1,''))[3]
    SecC <- unlist(strsplit(cod2,''))[3]
  }
  if (NumSub == 1) {nucleotides <- paste(FirstC,SecC,sep='_')}
  if (NumSub > 1)  {nucleotides  = 'MoreThanOne_SUBST'}
  return(nucleotides);
}

data$Subst = apply(as.matrix(data$FirstSecond), 1, COMPAR)  
withoutGapsCodonTable$Subs = data$Subst

tableSubs = withoutGapsCodonTable[withoutGapsCodonTable$Subs != 'MoreThanOne_SUBST',]

length(unique(tableSubs$Species)) #571

#### all subs

tableSubs$Subs = as.character(tableSubs$Subs)
mut = tableSubs[which(!(tableSubs$Subs %in% c('A_N', 'A_R', 'C_N', 'G_N', 'T_K', 'T_N', 'T_Y',
                                              '-_A', '-_C', 'C_Y', '-_G', 'G_R', '-_T', 'A_Y'))),]
# 749328

length(unique(mut$Species)) # 571

pdf('../results/nd6_22_01//mod_HistOfSubs.pdf')

ggplot(mut, aes(x = Subs)) +
  geom_histogram(stat = 'count')

Temp = mut
Temp$number = 1
Agg = aggregate(Temp$number, by = list(Temp$Subs), FUN = sum); names(Agg) = c('Subs','Number')
Agg$Number = Agg$Number/sum(Agg$Number)

ggbarplot(Agg, 'Subs', 'Number', xlab="Substitution types", title = 'all', 
          fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE)


# dev.off()


table(mut$Subs)

#### syn

syn = mut[as.character(mut$AncestralAA) == as.character(mut$DescendantAA),]
ggplot(syn, aes(x = Subs)) +
  geom_histogram(stat = 'count') + ggtitle('Synonymous')

####

write.table(mut, '../results/nd6_22_01/mod_AllSubs.nd6.txt', sep='\t', row.names = FALSE, quote = FALSE)
write.table(syn, '../results/nd6_22_01/mod_SynSubs.nd6.txt', sep='\t', row.names = FALSE, quote = FALSE)

################################################################
### 4fold deg

VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG',
                                     'AGT', 'AGC', 'AGA', 'AGG')
length(unique(VecOfSynFourFoldDegenerateSites)) # 36

mut4f = mut[mut$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & mut$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]; nrow(mut4f) # 209120
# ~300k subs

length(unique(mut4f$Species)) #571

Temp = mut4f
Temp$number = 1
Agg = aggregate(Temp$number, by = list(Temp$Subs), FUN = sum); names(Agg) = c('Subs','Number')
Agg$Number = Agg$Number/sum(Agg$Number)

ggbarplot(Agg, 'Subs', 'Number', xlab="Substitution types", title = '4-fold degenerate sites', 
          fill = 'Subs', color = 'Subs', palette = c("#bdbdbd", "#7fcdbb", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#feb24c", "#f03b20", "#bdbdbd", "#bdbdbd", "#bdbdbd", "#2c7fb8", "#bdbdbd"), combine = TRUE)

write.table(mut4f, '../results/nd6_22_01/mod_4foldSubs.nd6.txt', sep='\t', row.names = FALSE, quote = FALSE)

table(mut4f$Species)

dev.off()
