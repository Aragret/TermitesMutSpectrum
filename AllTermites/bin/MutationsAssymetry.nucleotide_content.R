rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

if(!require(dplyr)){install.packages('dplyr')}
if(!require(tidyr)){install.packages('tidyr')}
if(!require(ggplot2)){install.packages('ggplot2')}
if(!require(ggpubr)){install.packages('ggpubr')}
if(!require(cowplot)){install.packages('cowplot')}
if(!require(RColorBrewer)){install.packages('RColorBrewer')}

library(dplyr)
library(tidyr)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)


mutMajor = read.table('../results/nucleotide_content06_20/normMutSpecMajorStrand.txt',
                      header = TRUE, sep='\t')

mutMinor = read.table('../results/nucleotide_content06_20/normMutSpecMinorStrand.txt',
                      header = TRUE, sep='\t')

major = mutMajor %>%
  select(Species, Termites, Sociality, A_T_norm:C_G_norm) %>%
  mutate(
    from_A = A_T_norm + A_G_norm + A_C_norm,
    to_A = T_A_norm + G_A_norm + C_A_norm,
    from_T = T_A_norm + T_C_norm + T_G_norm,
    to_T = A_T_norm + C_T_norm + G_T_norm,
    from_C = C_A_norm + C_T_norm + C_G_norm,
    to_C = A_C_norm + T_C_norm + G_C_norm,
    from_G = G_A_norm + G_T_norm + G_C_norm,
    to_G = A_G_norm + T_G_norm + C_G_norm
  ) %>%
  mutate(
    A_asym = from_A / to_A, 
    T_asym = from_T / to_T,
    G_asym = from_G / to_G,
    C_asym = from_C / to_C
  )

minor = mutMinor %>%
  select(Species, Termites, Sociality, A_T_norm:C_G_norm) %>%
  mutate(
    from_A = T_A_norm + T_C_norm + T_G_norm, # complementary to match major strand
    to_A = A_T_norm + C_T_norm + G_T_norm,
    from_T = A_T_norm + A_G_norm + A_C_norm,
    to_T = T_A_norm + G_A_norm + C_A_norm,
    from_C = G_A_norm + G_T_norm + G_C_norm,
    to_C = A_G_norm + T_G_norm + C_G_norm,
    from_G = C_A_norm + C_T_norm + C_G_norm,
    to_G = A_C_norm + T_C_norm + G_C_norm
  ) %>% mutate(
    A_asym = from_A / to_A, 
    T_asym = from_T / to_T,
    G_asym = from_G / to_G,
    C_asym = from_C / to_C
  )

DFtallMinor <- minor %>% 
  select(Species, Termites, Sociality, from_A:to_G) %>%
  gather(key = Subs, value = Value, from_A:to_G)

DFtallMajor <- major %>% 
  select(Species, Termites, Sociality, from_A:to_G) %>%
  gather(key = Subs, value = Value, from_A:to_G)

DFtallMinorAsymm <- minor %>% 
  select(Species, Termites, Sociality, A_asym:C_asym) %>%
  gather(key = Subs, value = Value, A_asym:C_asym)

DFtallMajorAsymm <- major %>% 
  select(Species, Termites, Sociality, A_asym:C_asym) %>%
  gather(key = Subs, value = Value, A_asym:C_asym)


AllStrands = merge(DFtallMajor, DFtallMinor[, c('Species', 'Subs', 'Value')], by.x = c('Species', 'Subs'),
                   by.y = c('Species', 'Subs'))

AllStrandsAsymm = merge(DFtallMajorAsymm, DFtallMinorAsymm[, c('Species', 'Subs', 'Value')], by.x = c('Species', 'Subs'),
                        by.y = c('Species', 'Subs'))
                   
AllStrands$Value = AllStrands$Value.x + AllStrands$Value.y
AllStrandsAsymm$Value = AllStrandsAsymm$Value.x + AllStrandsAsymm$Value.y

AllStrands$Termites = as.factor(AllStrands$Termites)
DFtallMinor$Termites = as.factor(DFtallMinor$Termites)
DFtallMajor$Termites = as.factor(DFtallMajor$Termites)
AllStrandsAsymm$Termites = as.factor(AllStrandsAsymm$Termites)

AllStrands$Sociality = as.factor(AllStrands$Sociality)
DFtallMinor$Sociality = as.factor(DFtallMinor$Sociality)
DFtallMajor$Sociality = as.factor(DFtallMajor$Sociality)

my_colors <- RColorBrewer::brewer.pal(3, "Purples")[c(1,3)]

minor_cockroaches = ggbarplot(DFtallMinor, 'Subs', 'Value', xlab="Substitution types",
          fill = 'Termites',
          position = position_dodge(),
          add = 'mean_se',
          title = 'Minor genes') + 
  scale_fill_discrete(name = '', labels = c("Cockroaches", "Termites"))

major_cockroaches = ggbarplot(DFtallMajor, 'Subs', 'Value', xlab="Substitution types",
          fill = 'Termites',
          position = position_dodge(),
          add = 'mean_se',
          title = 'Major genes') + 
  scale_fill_discrete(name = '', labels = c("Cockroaches", "Termites")) 

all_cockroaches = ggbarplot(AllStrands, 'Subs', 'Value', xlab="Substitution types",
          fill = 'Termites',
          position = position_dodge(),
          add = 'mean_se',
          title = 'All genes',
          order = unique(DFtallMajor$Subs)) + 
  scale_fill_manual(values = my_colors, '', labels = c('Cockroaches', 'Termites'))

ggbarplot(AllStrandsAsymm[AllStrandsAsymm$Value != 'Inf',], 'Subs', 'Value', xlab="Substitution types",
          fill = 'Termites',
          position = position_dodge(),
          add = 'mean_se',
          title = 'All genes') +
          # order = unique(DFtallMajor$Subs)) + 
  scale_fill_discrete(name = '', labels = c("Cockroaches", "Termites"))

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

