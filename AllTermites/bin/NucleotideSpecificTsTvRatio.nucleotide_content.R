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

mutMinorComp = mutMinor

complSubs = unlist(lapply(names(mutMinorComp)[2:13], 
                          function(x){
                            y = paste(comp(s2c(x)[c(1, 3)], forceToLower = FALSE), 
                                      collapse = '_')
                            y
                          }
))

names(mutMinorComp)[2:13] = complSubs

allStrands = bind_rows(mutMajor %>% add_rownames(), 
                       mutMinorComp %>% add_rownames()) %>% 
  group_by(Species) %>% 
  select(Species:sumOfSubs) %>%
  summarise_all(sum)

allStrandsTsTv = allStrands %>%
  mutate(
    A_tstv = .$A_G / (.$A_C + .$A_G + .$A_T),
    T_tstv = .$T_C / (.$T_C + .$T_A + .$T_G),
    G_tstv = .$G_A / (.$G_A + .$G_T + .$G_C),
    C_tstv = .$C_T / (.$C_T + .$C_A + .$C_G)
  )

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

allStrandsTsTv = allStrandsTsTv %>% 
  mutate(
    G_tstv = replace_na(allStrandsTsTv$G_tstv, 0)
  )

DFtallMinor <- mutMinorTsTv %>% 
  select(Species, Sociality, A_tstv:C_tstv) %>%
  gather(key = Tstv, value = Value, A_tstv:C_tstv)

DFtallMajor <- mutMajorTsTv %>% 
  select(Species, Sociality, A_tstv:C_tstv) %>%
  gather(key = Tstv, value = Value, A_tstv:C_tstv)

DFtallAll = left_join(allStrandsTsTv %>% select(Species, A_tstv:C_tstv), 
                      mutMajor %>% select(Species, Taxonomy:Sociality)) %>%
  gather(key = Tstv, value = Value, A_tstv:C_tstv)


DFtallMinor$Sociality = as.factor(DFtallMinor$Sociality)
DFtallMajor$Sociality = as.factor(DFtallMajor$Sociality)
DFtallAll$Sociality = as.factor(DFtallAll$Sociality)

DFtallMajor = DFtallMajor[!is.na(DFtallMajor$Sociality),]
DFtallMinor = DFtallMinor[!is.na(DFtallMinor$Sociality),]
DFtallAll = DFtallAll[!is.na(DFtallAll$Sociality),]

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

all = ggboxplot(DFtallAll, 'Tstv', 'Value', xlab="Ts fractions",
                fill = 'Sociality',
                notch = TRUE,
                # position = position_dodge(),
                title = 'All genes') + 
  scale_fill_manual(name = '', labels = c("Cockroaches", "Less social termites", 'More social termites'),
                    values = RColorBrewer::brewer.pal(n = 3, name = "Purples")) +
  # scale_fill_brewer(palette = "Purples", breaks=c("Cockroaches", "Less sociale termites", 'More Social termites')) +
  scale_x_discrete(labels = c('A', 'T', 'G', 'C'), limits = c('T_tstv', 'A_tstv',
                                                              'C_tstv', 'G_tstv'))

plots = plot_grid(minor, major, all, nrow = 3)

save_plot('../results/nucleotide_content06_20/NucleotideSpecificTsTvRatio.pdf', plots, base_height = 12)

#################################################################################
### statistics for cockroaches, kalo and termitidae

data = left_join(allStrandsTsTv %>% select(Species, A_tstv:C_tstv), 
                 mutMajor %>% select(Species, Taxonomy:Sociality))

t.test(data[data$Termites == 0,]$A_tstv, data[data$Termites == 1,]$A_tstv)

# t = -12.03, df = 132.39, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.3104319 -0.2227632
# sample estimates:
#   mean of x mean of y 
# 0.3133276 0.5799251 

t.test(data[data$Termites == 0,]$T_tstv, data[data$Termites == 1,]$T_tstv)

# t = -6.2081, df = 144.88, p-value = 5.349e-09
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.17065962 -0.08823504
# sample estimates:
#   mean of x mean of y 
# 0.5278543 0.6573016

t.test(data[data$Termites == 0,]$G_tstv, data[data$Termites == 1,]$G_tstv)

# t = -4.9035, df = 114.61, p-value = 3.132e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.3080039 -0.1307575
# sample estimates:
#   mean of x mean of y 
# 0.5332236 0.7526043 

t.test(data[data$Termites == 0,]$C_tstv, data[data$Termites == 1,]$C_tstv)

# t = -0.65735, df = 125.09, p-value = 0.5122
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.04617999  0.02315170
# sample estimates:
#   mean of x mean of y 
# 0.6557153 0.6672295 

higher_termites = c("Apicotermitinae", "Cephalo-group", "Microcerotermes", 
                    "Termes-group", "Nasutitermitinae", "Amitermes-group", 
                    "Promiro", "Macrotermitinae", "Cubitermitinae", "Foraminitermitinae", 
                    "Syntermitinae", "Sphaerotermitinae", "Neocapri-group", 
                    'Pericapritermes-group', "pericapritermes-group")

t.test(data[data$Taxonomy == 'Kalotermitidae',]$A_tstv, 
       data[data$Taxonomy %in% higher_termites,]$A_tstv)

# t = -5.5234, df = 75.168, p-value = 4.561e-07
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.20171763 -0.09478359
# sample estimates:
#   mean of x mean of y 
# 0.4441776 0.5924282

t.test(data[data$Taxonomy == 'Kalotermitidae',]$T_tstv, 
       data[data$Taxonomy %in% higher_termites,]$T_tstv)

# t = -2.6071, df = 81.588, p-value = 0.01085
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.12131989 -0.01630111
# sample estimates:
#   mean of x mean of y 
# 0.5831763 0.6519868 

t.test(data[data$Taxonomy == 'Kalotermitidae',]$G_tstv, 
       data[data$Taxonomy %in% higher_termites,]$G_tstv)

# t = -1.3991, df = 77.032, p-value = 0.1658
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1453222  0.0253827
# sample estimates:
#   mean of x mean of y 
# 0.7014977 0.7614675 

t.test(data[data$Taxonomy == 'Kalotermitidae',]$C_tstv, 
       data[data$Taxonomy %in% higher_termites,]$C_tstv)

# t = -0.14621, df = 78.149, p-value = 0.8841
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.04225638  0.03647408
# sample estimates:
#   mean of x mean of y 
# 0.6669161 0.6698073 

# add lower and higher termites

higher_termites = c("Apicotermitinae", "Cephalo-group", "Microcerotermes", 
                    "Termes-group", "Nasutitermitinae", "Amitermes-group", 
                    "Promiro", "Macrotermitinae", "Cubitermitinae", "Foraminitermitinae", 
                    "Syntermitinae", "Sphaerotermitinae", "Neocapri-group", 
                    'Pericapritermes-group', "pericapritermes-group")

lowerTermites = data[data$Termites == 1 & !(data$Taxonomy %in% higher_termites),]
higherTermites = data[data$Taxonomy %in% higher_termites,]

t.test(lowerTermites$A_tstv, higherTermites$A_tstv)

# t = -1.8413, df = 181.28, p-value = 0.06722
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.080095907  0.002769091
# sample estimates:
#   mean of x mean of y 
# 0.5537648 0.5924282 

t.test(lowerTermites$T_tstv, higherTermites$T_tstv)

# t = 0.79555, df = 203.86, p-value = 0.4272
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.02468517  0.05808021
# sample estimates:
#   mean of x mean of y 
# 0.6686844 0.6519868 

t.test(lowerTermites$G_tstv, higherTermites$G_tstv)

# t = -0.93722, df = 213.21, p-value = 0.3497
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.08505146  0.03023613
# sample estimates:
#   mean of x mean of y 
# 0.7340598 0.7614675 

t.test(lowerTermites$C_tstv, higherTermites$C_tstv)

# t = -0.56255, df = 208.86, p-value = 0.5743
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.03590626  0.01996348
# sample estimates:
#   mean of x mean of y 
# 0.6618359 0.6698073 


# plots

to_plot = DFtallAll %>%
  filter(Termites == 1) %>%
  mutate(
    Higher = as.factor(case_when(.$Taxonomy %in% higherTermites$Taxonomy ~ 1,
                       .$Taxonomy %in% lowerTermites$Taxonomy ~ 0))
  )


all_higher = ggboxplot(to_plot, 'Tstv', 'Value', xlab="Ts fractions",
                fill = 'Higher',
                notch = TRUE,
                # position = position_dodge(),
                title = 'All genes') + 
  scale_fill_manual(name = '', labels = c("Lower termites", "Higher termites"),
                    values = RColorBrewer::brewer.pal(n = 3, name = "Purples")) +
  # scale_fill_brewer(palette = "Purples", breaks=c("Cockroaches", "Less sociale termites", 'More Social termites')) +
  scale_x_discrete(labels = c('A', 'T', 'G', 'C'), limits = c('T_tstv', 'A_tstv',
                                                              'C_tstv', 'G_tstv'))

all_higher
