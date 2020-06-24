rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

if(!require(dplyr)){install.packages('dplyr')}
if(!require(tidyr)){install.packages('tidyr')}
if(!require(ggpubr)){install.packages('ggpubr')}
if(!require(cowplot)){install.packages('cowplot')}
if(!require(seqinr)){install.packages('seqinr')}

library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
library(seqinr)

mutMajor = read.table('../results/nucleotide_content06_20/normMutSpecMajorStrand.txt',
                      header = TRUE, sep='\t')

mutMinor = read.table('../results/nucleotide_content06_20/normMutSpecMinorStrand.txt',
                      header = TRUE, sep='\t')

mutMajorTsTv = mutMajor %>%
  mutate(
    A_tstv = .$A_G / (.$A_C + .$A_G + .$A_T),
    T_tstv = .$T_C / (.$T_C + .$T_A + .$T_G),
    G_tstv = .$G_A / (.$G_A + .$G_T + .$G_C),
    C_tstv = .$C_T / (.$C_T + .$C_A + .$C_G)
  )

mutMinorTsTv = mutMinor %>%
  mutate(
    A_tstv = .$A_G / (.$A_C + .$A_G + .$A_T),
    T_tstv = .$T_C / (.$T_C + .$T_A + .$T_G),
    G_tstv = .$G_A / (.$G_A + .$G_T + .$G_C),
    C_tstv = .$C_T / (.$C_T + .$C_A + .$C_G)
  )

summary(mutMajorTsTv$A_tstv)
summary(mutMajorTsTv$T_tstv)
summary(mutMajorTsTv$G_tstv)
summary(mutMajorTsTv$C_tstv)

mutMajorTsTv = mutMajorTsTv %>% 
  mutate(
    T_tstv = replace_na(mutMajorTsTv$T_tstv, 0),
    G_tstv = replace_na(mutMajorTsTv$G_tstv, 0)
  )

mutMinorTsTv = mutMinorTsTv %>% 
  mutate(
    A_tstv = replace_na(mutMinorTsTv$A_tstv, 0),
    G_tstv = replace_na(mutMinorTsTv$G_tstv, 0),
    C_tstv = replace_na(mutMinorTsTv$C_tstv, 0)
  )


DFtallMinor <- mutMinorTsTv %>% 
  select(Species, Sociality, A_tstv:C_tstv) %>%
  gather(key = Tstv, value = Value, A_tstv:C_tstv)

DFtallMajor <- mutMajorTsTv %>% 
  select(Species, Sociality, A_tstv:C_tstv) %>%
  gather(key = Tstv, value = Value, A_tstv:C_tstv)


complSubs = unlist(lapply(DFtallMinor$Tstv, 
                          function(x){
                            y = paste(toupper(comp(s2c(x)[1])), collapse = '_')
                            y
                          }
))

DFtallMinor$ComplSubs = complSubs

DFtallMajor$Tstv = sub('_tstv', '', DFtallMajor$Tstv)

AllStrands = merge(DFtallMajor, DFtallMinor[, c('Species', 'ComplSubs', 'Value')], by.x = c('Species', 'Tstv'),
                   by.y = c('Species', 'ComplSubs'))

AllStrands$Value = AllStrands$Value.x + AllStrands$Value.y

DFtallMinor$Sociality = as.factor(DFtallMinor$Sociality)
DFtallMajor$Sociality = as.factor(DFtallMajor$Sociality)
AllStrands$Sociality = as.factor(AllStrands$Sociality)

DFtallMajor = DFtallMajor[!is.na(DFtallMajor$Sociality),]
DFtallMinor = DFtallMinor[!is.na(DFtallMinor$Sociality),]
AllStrands = AllStrands[!is.na(AllStrands$Sociality),]

minor = ggboxplot(DFtallMinor, 'Tstv', 'Value', xlab="Ts fractions",
                  fill = 'Sociality',
                  notch = TRUE,
                  # position = position_dodge(),
                  title = 'Minor genes') + 
  scale_fill_manual(name = '', labels = c("Cockroaches", "Less social termites", 'More social termites'),
                    values = RColorBrewer::brewer.pal(n = 3, name = "Purples")) +
  # scale_fill_brewer(palette = "Purples", breaks=c("Cockroaches", "Less sociale termites", 'More Social termites')) +
  scale_x_discrete(labels = c('A', 'T', 'G', 'C'))

major = ggboxplot(DFtallMajor, 'Tstv', 'Value', xlab="Ts fractions",
                  fill = 'Sociality',
                  notch = TRUE,
                  # position = position_dodge(),
                  title = 'Major genes') + 
  scale_fill_manual(name = '', labels = c("Cockroaches", "Less social termites", 'More social termites'),
                    values = RColorBrewer::brewer.pal(n = 3, name = "Purples")) +
  # scale_fill_brewer(palette = "Purples", breaks=c("Cockroaches", "Less sociale termites", 'More Social termites')) +
  scale_x_discrete(labels = c('A', 'T', 'G', 'C'))

all = ggboxplot(AllStrands, 'Tstv', 'Value', xlab="Ts fractions",
                fill = 'Sociality',
                notch = TRUE,
                # position = position_dodge(),
                title = 'All genes') + 
  scale_fill_manual(name = '', labels = c("Cockroaches", "Less social termites", 'More social termites'),
                    values = RColorBrewer::brewer.pal(n = 3, name = "Purples")) +
  # scale_fill_brewer(palette = "Purples", breaks=c("Cockroaches", "Less sociale termites", 'More Social termites')) +
  scale_x_discrete(labels = c('A', 'T', 'G', 'C'))

plots = plot_grid(minor, major, all, nrow = 3)

save_plot('../results/nucleotide_content06_20/NucleotideSpecificTsTvRatio.pdf', plots, base_height = 12)
