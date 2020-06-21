rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

if(!require(dplyr)){install.packages('dplyr')}
if(!require(tidyr)){install.packages('tidyr')}
if(!require(ggplot2)){install.packages('ggplot2')}
if(!require(ggpubr)){install.packages('ggpubr')}

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

mut = read.table('../results/nd6_22_01/mod_PolarizeMutations.CodonsTable.nd6.txt',
                 header=TRUE, sep='\t')

syn = mut[mut$AncestralAA == mut$DescendantAA,]

VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG',
                                     'AGT', 'AGC', 'AGA', 'AGG')

majorGenes = c('ATP6', 'ATP8', 'COI', 'COII', 'COIII', 'CyB', 'ND3', 'ND6', 'ND2')
minorGenes = c('ND1', 'ND4', 'ND4L', 'ND5')

cockroaches = c('Ectobiidae1', 'Tryonicidae', 'Blaberidae', 'Corydiidae', 'Ectobiidae2',
                'Lamproblattidae', 'Anaplectidae', 'Blattidae', 'Cryptocercidae', 'Ectobiidae3',
                'Nocticolidae')

lessSocial = c('Kalotermitidae', 'Stolotermitidae', 'Serritermitidae',
               'Archotermopsidae', 'Termitogetoninae', 'Prorhinotermitinae')

moreSocial = c('Coptotermitinae', 'Nasutitermitinae', 'Amitermes-group',
               'Macrotermitinae', 'Apicotermitinae', 'Rhinotermitinae',
               'Rhinotermitidae', 'Microcerotermes', 'Termes-group', 'Promiro',
               'Cubitermitinae', 'Foraminitermitinae', 'Syntermitinae',
               'Pericapritermes-group', 'Sphaerotermitinae', 'Neocapri-group',
               'Cephalo-group', 'Mastotermitidae', 'Prorhinotermitinae',
               'pericapritermes-group')


# ATP6.fasta.align.fna = 1-687
# ATP8.fasta.align.fna = 688-867
# COI.fasta.align.fna = 868-2433
# COII.fasta.align.fna = 2434-3120
# COIII.fasta.align.fna = 3121-3918
# CyB.fasta.align.fna = 3919-5091
# ND1.fasta.align.fna = 5092-6051
# ND2.fasta.align.fna = 6052-7128
# ND3.fasta.align.fna = 7129-7488
# ND4.fasta.align.fna = 7489-8871
# ND4L.fasta.align.fna = 8872-9171
# ND5.fasta.align.fna = 9172-10968
# ND6.fasta.align.fna = 10969-11556

syn$NuclPosition = syn$CodonPosition*2 + syn$CodonPosition - 2
syn = syn %>% mutate(MajorPosition = case_when(.$NuclPosition < 5091 | 
                                                 (.$NuclPosition > 6051 & .$NuclPosition < 7480) |
                                                 .$NuclPosition > 10968 ~ 1))

syn$MajorPosition[is.na(syn$MajorPosition)] <- 0

fourFold = syn[syn$AncestorCodon %in% VecOfSynFourFoldDegenerateSites & 
                 syn$DescendantCodon %in% VecOfSynFourFoldDegenerateSites,]

thirdPosAnc = substring(fourFold$AncestorCodon, first = 3, last = 3)
thirdPosDes = substring(fourFold$DescendantCodon, first = 3, last = 3)

fourFold$Subs = paste(thirdPosAnc, thirdPosDes, sep='_')

table(fourFold$Subs)

fourFold = fourFold[!(fourFold$Subs %in% c('A_A', 'C_C', 'T_T')),]

major = fourFold[fourFold$MajorPosition == 1,]
minor = fourFold[fourFold$MajorPosition == 0,]

### calculate mut spectrum for major genes
speciesNumber = length(unique(major$Species))
mutTableMajor = setNames(data.frame(matrix(0, ncol = 12, nrow = speciesNumber)), unique(major$Subs))
Species = unique(major$Species)
mutTableMajor = cbind(Species, mutTableMajor)

for(i in 1:nrow(major)){
  row = major[i, ]
  mutTableMajor[mutTableMajor$Species == row$Species, row$Subs] = mutTableMajor[mutTableMajor$Species == row$Species, row$Subs] + 1
}

table(major[major$Species == "AUS49_Tumulitermes_sp._1",]$Subs)
mutTableMajor[mutTableMajor$Species == 'AUS49_Tumulitermes_sp._1',]

### calculate mut spectrum for minor genes
speciesNumber = length(unique(minor$Species))
mutTableMinor = setNames(data.frame(matrix(0, ncol = 12, nrow = speciesNumber)), unique(minor$Subs))
Species = unique(minor$Species)
mutTableMinor = cbind(Species, mutTableMinor)

for(i in 1:nrow(minor)){
  row = minor[i, ]
  mutTableMinor[mutTableMinor$Species == row$Species, row$Subs] = mutTableMinor[mutTableMinor$Species == row$Species, row$Subs] + 1
}

table(minor[minor$Species == "AUS49_Tumulitermes_sp._1",]$Subs)
mutTableMinor[mutTableMinor$Species == 'AUS49_Tumulitermes_sp._1',]

### get fractions of mutations for major genes

a = matrix(0, ncol=12, nrow=nrow(mutTableMajor))

mutTableMajor$sumOfSubs = mutTableMajor[, "C_T"] + mutTableMajor[, "A_G"] + mutTableMajor[, "T_C"] + mutTableMajor[, "G_T"] + mutTableMajor[, "A_C"] +
  mutTableMajor[, "G_A"] + mutTableMajor[, "A_T"] + mutTableMajor[, "G_C"] + mutTableMajor[, "C_A"] + mutTableMajor[, "C_G"] + mutTableMajor[, "T_A"] +
  mutTableMajor[, "T_G"]

for(i in 1:12){
  a[, i] = mutTableMajor[, i+1] / mutTableMajor$sumOfSubs
}

fractionsMajor = as.data.frame(a)
names(fractionsMajor) = sub(' ', '', paste(names(mutTableMajor[, 2:13]), '_fr'))

mutFrMajor = cbind(mutTableMajor, fractionsMajor)


### get fractions of mutations for minor genes

a = matrix(0, ncol=12, nrow=nrow(mutTableMinor))

mutTableMinor$sumOfSubs = mutTableMinor[, "C_T"] + mutTableMinor[, "A_G"] + mutTableMinor[, "T_C"] + mutTableMinor[, "G_T"] + mutTableMinor[, "A_C"] +
  mutTableMinor[, "G_A"] + mutTableMinor[, "A_T"] + mutTableMinor[, "G_C"] + mutTableMinor[, "C_A"] + mutTableMinor[, "C_G"] + mutTableMinor[, "T_A"] +
  mutTableMinor[, "T_G"]

for(i in 1:12){
  a[, i] = mutTableMinor[, i+1] / mutTableMinor$sumOfSubs
}

fractionsMinor = as.data.frame(a)
names(fractionsMinor) = sub(' ', '', paste(names(mutTableMinor[, 2:13]), '_fr'))

mutFrMinor = cbind(mutTableMinor, fractionsMinor)

####################################################################################
### get normalized frequencies

count = read.table('../results/nucleotide_content06_20/ATGCforEachGene4fold.txt', 
                   header = TRUE, sep='\t')

count$Species = sub(';', '', count$Species)
count$Species = unlist(lapply(count$Species, paste, '_1', sep=''))
count$Species = sub('__', '_', count$Species)
count$Species = sub('?', '', count$Species, fixed = TRUE)
count$Species = sub('___', '__', count$Species)


countMajor = count[count$MajorStrand == 1,]
countMinor = count[count$MajorStrand == 0,]

### get all nucleotide counts for each species

aggMajor = countMajor %>% group_by(Species) %>% summarise_at(vars(A, T, G, C), sum)

# sum(countMajor[countMajor$Species == '3.10.3_Amitermes_capito',]$A) # 796

aggMinor = countMinor %>% group_by(Species) %>% summarise_at(vars(A, T, G, C), sum)

### get nucleotide frequencies

aggMajorFreq = aggMajor %>%
  mutate(
    A_fr = A / (A + T + G + C),
    T_fr = T / (A + T + G + C),
    G_fr = G / (A + T + G + C),
    C_fr = C / (A + T + G + C))

aggMinorFreq = aggMinor %>%
  mutate(
    A_fr = A / (A + T + G + C),
    T_fr = T / (A + T + G + C),
    G_fr = G / (A + T + G + C),
    C_fr = C / (A + T + G + C))

mutWithNuclMajor = merge(mutFrMajor, aggMajorFreq, by='Species')
mutWithNuclMinor = merge(mutFrMinor, aggMinorFreq, by='Species')

setdiff(mutFrMajor$Species, mutWithNuclMajor$Species) # 0
setdiff(mutFrMinor$Species, mutWithNuclMinor$Species) # 0

### normalization

mutWithNuclMinorNormalized = mutWithNuclMinor %>%
  mutate(
    A_T_norm = A_T_fr / A_fr,
    A_G_norm = A_G_fr / A_fr,
    A_C_norm = A_C_fr / A_fr,
    T_A_norm = T_A_fr / T_fr,
    T_G_norm = T_G_fr / T_fr,
    T_C_norm = T_C_fr / T_fr,
    G_A_norm = G_A_fr / G_fr,
    G_T_norm = G_T_fr / G_fr,
    G_C_norm = G_C_fr / G_fr,
    C_A_norm = C_A_fr / C_fr,
    C_T_norm = C_T_fr / C_fr,
    C_G_norm = C_G_fr / C_fr
  )

mutWithNuclMajorNormalized = mutWithNuclMajor %>%
  mutate(
    A_T_norm = A_T_fr / A_fr,
    A_G_norm = A_G_fr / A_fr,
    A_C_norm = A_C_fr / A_fr,
    T_A_norm = T_A_fr / T_fr,
    T_G_norm = T_G_fr / T_fr,
    T_C_norm = T_C_fr / T_fr,
    G_A_norm = G_A_fr / G_fr,
    G_T_norm = G_T_fr / G_fr,
    G_C_norm = G_C_fr / G_fr,
    C_A_norm = C_A_fr / C_fr,
    C_T_norm = C_T_fr / C_fr,
    C_G_norm = C_G_fr / C_fr
  )

### merge with families 
families = read.table('../results/nd6_22_01/mutSpectrumFractions.txt', header = TRUE, sep='\t')
families = families[, c('Species', 'Taxonomy')]

mutWithNuclMinorNormalized = merge(mutWithNuclMinorNormalized, families, by='Species', all.x = TRUE)
mutWithNuclMajorNormalized = merge(mutWithNuclMajorNormalized, families, by='Species', all.x = TRUE)

mutWithNuclMinorNormalized = mutWithNuclMinorNormalized %>%
  filter(!is.na(Taxonomy)) %>%
  mutate(
    Cockroaches = case_when(.$Taxonomy %in% cockroaches ~ 1,
                           !(.$Taxonomy %in% cockroaches) ~ 0)
  )

mutWithNuclMajorNormalized = mutWithNuclMajorNormalized %>%
  filter(!is.na(Taxonomy)) %>%
  mutate(
    Cockroaches = case_when(.$Taxonomy %in% cockroaches ~ 1,
                            !(.$Taxonomy %in% cockroaches) ~ 0)
  )

DFtallMinor <- mutWithNuclMinorNormalized %>% 
  select(Species, Cockroaches, A_T_norm:C_G_norm) %>%
  gather(key = Subs, value = Value, A_T_norm:C_G_norm)

DFtallMajor <- mutWithNuclMajorNormalized %>% 
  select(Species, Cockroaches, A_T_norm:C_G_norm) %>%
  gather(key = Subs, value = Value, A_T_norm:C_G_norm)

DFtallMinor$Cockroaches = as.factor(DFtallMinor$Cockroaches)
DFtallMajor$Cockroaches = as.factor(DFtallMajor$Cockroaches)

minor_cockroaches = 
  ggbarplot(DFtallMinor, 'Subs', 'Value', xlab="Substitution types",
          fill = 'Cockroaches',
          position = position_dodge(),
          add = 'mean_se',
          title = 'Minor genes') + 
  scale_fill_discrete(name = '', labels = c("Termites", "Cockroaches")) +
  scale_x_discrete(labels = sub('_norm', '', unique(DFtallMinor$Subs)))

major_cockroaches = 
  ggbarplot(DFtallMajor, 'Subs', 'Value', xlab="Substitution types",
          fill = 'Cockroaches',
          position = position_dodge(),
          add = 'mean_se',
          title = 'Major genes') + 
  scale_fill_discrete(name = '', labels = c("Termites", "Cockroaches")) + 
  scale_x_discrete(labels = sub('_norm', '', unique(DFtallMajor$Subs)))


DFtallMajor[DFtallMajor$Value == max(DFtallMajor$Value),] # Reticulitermes_flavipes
DFtallMinor[DFtallMinor$Value == max(DFtallMinor$Value),] # Cryptocercus_meridianus

##############################################################################
### barplots for cockroaches, less social and more social termites

mutWithNuclMinorNormalized = mutWithNuclMinorNormalized %>%
  filter(!is.na(Taxonomy)) %>%
  mutate(
    Sociality = case_when(
      .$Cockroaches == 1 ~ 0,
      .$Taxonomy %in% lessSocial ~ 1,
      .$Taxonomy %in% moreSocial ~ 2
    )
  )

mutWithNuclMajorNormalized = mutWithNuclMajorNormalized %>%
  filter(!is.na(Taxonomy)) %>%
  mutate(
    Sociality = case_when(
      .$Cockroaches == 1 ~ 0,
      .$Taxonomy %in% lessSocial ~ 1,
      .$Taxonomy %in% moreSocial ~ 2
    )
  )

DFtallMinor <- mutWithNuclMinorNormalized %>% 
  select(Species, Sociality, A_T_norm:C_G_norm) %>%
  gather(key = Subs, value = Value, A_T_norm:C_G_norm)

DFtallMajor <- mutWithNuclMajorNormalized %>% 
  select(Species, Sociality, A_T_norm:C_G_norm) %>%
  gather(key = Subs, value = Value, A_T_norm:C_G_norm)

DFtallMinor$Sociality = as.factor(DFtallMinor$Sociality)
DFtallMajor$Sociality = as.factor(DFtallMajor$Sociality)

DFtallMajor = DFtallMajor[!is.na(DFtallMajor$Sociality),]
DFtallMinor = DFtallMinor[!is.na(DFtallMinor$Sociality),]

minor = ggbarplot(DFtallMinor, 'Subs', 'Value', xlab="Substitution types",
          fill = 'Sociality',
          position = position_dodge(),
          add = 'mean_se',
          title = 'Minor genes') + 
  scale_fill_manual(name = '', labels = c("Cockroaches", "Less sociale termites", 'More Social termites'),
                      values = RColorBrewer::brewer.pal(n = 3, name = "Purples")) +
  # scale_fill_brewer(palette = "Purples", breaks=c("Cockroaches", "Less sociale termites", 'More Social termites')) +
  scale_x_discrete(labels = sub('_norm', '', unique(DFtallMinor$Subs))) 

major = ggbarplot(DFtallMajor, 'Subs', 'Value', xlab="Substitution types",
          fill = 'Sociality',
          position = position_dodge(),
          add = 'mean_se',
          title = 'Major genes') + 
  scale_fill_manual(name = '', labels = c("Cockroaches", "Less sociale termites", 'More Social termites'),
                    values = RColorBrewer::brewer.pal(n = 3, name = "Purples")) +
  # scale_fill_brewer(palette = "Purples", breaks=c("Cockroaches", "Less sociale termites", 'More Social termites')) +
  scale_x_discrete(labels = sub('_norm', '', unique(DFtallMajor$Subs))) 

plots = plot_grid(major_cockroaches, minor_cockroaches, major, minor, nrow = 2)

save_plot('../results/nucleotide_content06_20/normMutSpecCockroachesTermites.pdf', plots, base_height = 10)

write.table(mutWithNuclMajorNormalized, '../results/nucleotide_content06_20/normMutSpecMajorStrand.txt', 
            sep='\t', row.names = FALSE, quote = FALSE)
write.table(mutWithNuclMinorNormalized, '../results/nucleotide_content06_20/normMutSpecMinorStrand.txt', 
            sep='\t', row.names = FALSE, quote = FALSE)

write.table(DFtallMajor, '../results/nucleotide_content06_20/normMutSpecMajorStrandLong.txt', 
            sep='\t', row.names = FALSE, quote = FALSE)
write.table(DFtallMinor, '../results/nucleotide_content06_20/normMutSpecMinorStrandLong.txt', 
            sep='\t', row.names = FALSE, quote = FALSE)
