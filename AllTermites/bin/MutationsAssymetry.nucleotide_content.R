rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(cowplot)

mutMajor = read.table('../results/nucleotide_content06_20/normMutSpecMajorStrand.txt',
                      header = TRUE, sep='\t')

mutMinor = read.table('../results/nucleotide_content06_20/normMutSpecMinorStrand.txt',
                      header = TRUE, sep='\t')

major = mutMajor %>%
  select(Species, Cockroaches, Sociality, A_T_norm:C_G_norm) %>%
  mutate(
    from_A = A_T_norm + A_G_norm + A_C_norm,
    to_A = T_A_norm + G_A_norm + C_A_norm,
    from_T = T_A_norm + T_C_norm + T_G_norm,
    to_T = A_T_norm + C_T_norm + G_T_norm,
    from_C = C_A_norm + C_T_norm + C_G_norm,
    to_C = A_C_norm + T_C_norm + G_C_norm,
    from_G = G_A_norm + G_T_norm + G_C_norm,
    to_G = A_G_norm + T_G_norm + C_G_norm
  )

minor = mutMinor %>%
  select(Species, Cockroaches, Sociality, A_T_norm:C_G_norm) %>%
  mutate(
    from_A = T_A_norm + T_C_norm + T_G_norm, # complementary to match major strand
    to_A = A_T_norm + C_T_norm + G_T_norm,
    from_T = A_T_norm + A_G_norm + A_C_norm,
    to_T = T_A_norm + G_A_norm + C_A_norm,
    from_C = G_A_norm + G_T_norm + G_C_norm,
    to_C = A_G_norm + T_G_norm + C_G_norm,
    from_G = C_A_norm + C_T_norm + C_G_norm,
    to_G = A_C_norm + T_C_norm + G_C_norm
  )

DFtallMinor <- minor %>% 
  select(Species, Cockroaches, Sociality, from_A:to_G) %>%
  gather(key = Subs, value = Value, from_A:to_G)

DFtallMajor <- major %>% 
  select(Species, Cockroaches, Sociality, from_A:to_G) %>%
  gather(key = Subs, value = Value, from_A:to_G)

AllStrands = merge(DFtallMajor, DFtallMinor[, c('Species', 'Subs', 'Value')], by.x = c('Species', 'Subs'),
                   by.y = c('Species', 'Subs'))

AllStrands$Value = AllStrands$Value.x + AllStrands$Value.y

AllStrands$Cockroaches = as.factor(AllStrands$Cockroaches)
DFtallMinor$Cockroaches = as.factor(DFtallMinor$Cockroaches)
DFtallMajor$Cockroaches = as.factor(DFtallMajor$Cockroaches)

AllStrands$Sociality = as.factor(AllStrands$Sociality)
DFtallMinor$Sociality = as.factor(DFtallMinor$Sociality)
DFtallMajor$Sociality = as.factor(DFtallMajor$Sociality)

minor_cockroaches = ggbarplot(DFtallMinor, 'Subs', 'Value', xlab="Substitution types",
          fill = 'Cockroaches',
          position = position_dodge(),
          add = 'mean_se',
          title = 'Minor genes') + 
  scale_fill_discrete(name = '', labels = c("Termites", "Cockroaches"))

major_cockroaches = ggbarplot(DFtallMajor, 'Subs', 'Value', xlab="Substitution types",
          fill = 'Cockroaches',
          position = position_dodge(),
          add = 'mean_se',
          title = 'Major genes') + 
  scale_fill_discrete(name = '', labels = c("Termites", "Cockroaches")) 

all_cockroaches = ggbarplot(AllStrands, 'Subs', 'Value', xlab="Substitution types",
          fill = 'Cockroaches',
          position = position_dodge(),
          add = 'mean_se',
          title = 'All genes',
          order = unique(DFtallMajor$Subs)) + 
  scale_fill_discrete(name = '', labels = c("Termites", "Cockroaches"))


minor = ggbarplot(DFtallMinor[!is.na(DFtallMinor$Sociality),], 'Subs', 'Value', xlab="Substitution types",
          fill = 'Sociality',
          position = position_dodge(),
          add = 'mean_se',
          title = 'Minor genes') + 
  scale_fill_manual(name = '', labels = c("Cockroaches", "Less social termites", 'More social termites'),
                    values = RColorBrewer::brewer.pal(n = 3, name = "Purples")) +
  # scale_fill_brewer(palette = "Purples", breaks=c("Cockroaches", "Less sociale termites", 'More Social termites')) +
  scale_x_discrete(labels = sub('_norm', '', unique(DFtallMinor$Subs))) 

major = ggbarplot(DFtallMajor[!is.na(DFtallMajor$Sociality),], 'Subs', 'Value', xlab="Substitution types",
                  fill = 'Sociality',
                  position = position_dodge(),
                  add = 'mean_se',
                  title = 'Major genes') + 
  scale_fill_manual(name = '', labels = c("Cockroaches", "Less social termites", 'More social termites'),
                    values = RColorBrewer::brewer.pal(n = 3, name = "Purples")) +
  # scale_fill_brewer(palette = "Purples", breaks=c("Cockroaches", "Less sociale termites", 'More Social termites')) +
  scale_x_discrete(labels = sub('_norm', '', unique(DFtallMajor$Subs))) 

allStrands = ggbarplot(AllStrands[!is.na(AllStrands$Sociality),], 'Subs', 'Value', xlab="Substitution types",
                       fill = 'Sociality',
                       position = position_dodge(),
                       add = 'mean_se',
                       title = 'All genes',
                       order = unique(DFtallMajor$Subs)) + 
  scale_fill_manual(name = '', labels = c("Cockroaches", "Less social termites", 'More social termites'),
                    values = RColorBrewer::brewer.pal(n = 3, name = "Purples"))
allStrands

plots1 = plot_grid(major_cockroaches, minor_cockroaches, all_cockroaches, nrow = 2)
plots2 = plot_grid(major, minor, allStrands, nrow = 2)

save_plot('../results/nucleotide_content06_20/normMutAsymmetryRoaches.pdf', plots1,
          base_height = 8)
save_plot('../results/nucleotide_content06_20/normMutAsymmetrySociality.pdf', plots2,
          base_height = 8)

