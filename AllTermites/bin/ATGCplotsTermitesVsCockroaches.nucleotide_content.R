rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

if(!require(dplyr)){install.packages('dplyr')}
if(!require(dplyr)){install.packages('ggplot2')}
if(!require(dplyr)){install.packages('cowplot')}

library(dplyr)
library(ggplot2)
library(cowplot)


cont = read.table('../results/nucleotide_content06_20/ATGCforEachGene4fold.txt', 
                  header = TRUE, sep='\t')
cont$Species = sub(';', '', cont$Species)
cont$Species = unlist(lapply(cont$Species, paste, '_1', sep=''))
cont$Species = sub('__', '_', cont$Species)
cont$Species = sub('?', '', cont$Species, fixed = TRUE)
cont$Species = sub('___', '__', cont$Species)


mut = read.table('../results/nd6_22_01/mutSpectrumFractions.txt', header = TRUE, sep='\t')
families = mut[, c('Species', 'Taxonomy')]

data = merge(cont, families, by='Species')

setdiff(families$Species, data$Species)

boxplot(data$A, data$T, data$G, data$C)

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

data = data[!is.na(data$Taxonomy),] # remove Mantids

data$A_fr = data$A / (data$A + data$T + data$G + data$C)
data$T_fr = data$T / (data$A + data$T + data$G + data$C)
data$G_fr = data$G / (data$A + data$T + data$G + data$C)
data$C_fr = data$C / (data$A + data$T + data$G + data$C)


for(i in 1:nrow(data)){
  if(data$Taxonomy[i] %in% cockroaches){
    data$Cockroaches[i] = 1
  }
  if(!(data$Taxonomy[i] %in% cockroaches))
  {data$Cockroaches[i] = 0}
}

data$Cockroaches = as.factor(data$Cockroaches)

major_strand = data[data$MajorStrand == 1,]

minor_strand = data[data$MajorStrand == 0,]

a_major = ggplot(major_strand, aes(Cockroaches, A_fr, fill = Cockroaches)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Purples") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Termites', 'Cockroaches')) +
  xlab('') + ylab('') + ggtitle('A in major genes')

a_major

t_major = ggplot(major_strand, aes(Cockroaches, T_fr, fill = Cockroaches)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Termites', 'Cockroaches')) +
  xlab('') + ylab('') + ggtitle('T in major genes')

t_major

g_major = ggplot(major_strand, aes(Cockroaches, G_fr, fill = Cockroaches)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Termites', 'Cockroaches')) +
  xlab('') + ylab('') + ggtitle('G in major genes')

g_major

c_major = ggplot(major_strand, aes(Cockroaches, C_fr, fill = Cockroaches)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Termites', 'Cockroaches')) +
  xlab('') + ylab('') + ggtitle('C in major genes')

c_major

#################################
# minor genes

a_minor = ggplot(minor_strand, aes(Cockroaches, A_fr, fill = Cockroaches)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Termites', 'Cockroaches')) +
  xlab('') + ylab('') + ggtitle('A in minor genes')

a_minor

t_minor = ggplot(minor_strand, aes(Cockroaches, T_fr, fill = Cockroaches)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Termites', 'Cockroaches')) +
  xlab('') + ylab('') + ggtitle('T in minor genes')

t_minor

g_minor = ggplot(minor_strand, aes(Cockroaches, G_fr, fill = Cockroaches)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Termites', 'Cockroaches')) +
  xlab('') + ylab('') + ggtitle('G in minor genes')

g_minor

c_minor = ggplot(minor_strand, aes(Cockroaches, C_fr, fill = Cockroaches)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Termites', 'Cockroaches')) +
  xlab('') + ylab('') + ggtitle('C in minor genes')

c_minor

plots = plot_grid(a_major, t_major, g_major, c_major, a_minor, t_minor, g_minor, c_minor, nrow = 2)

save_plot('../results/nucleotide_content06_20/ATGCplotsTermitesVsCockroaches.pdf', plots, base_height = 7)

#################################################################################
### Compare cockroaches, less and more social termites

data = data %>% 
  mutate(
    Sociality = case_when(
      .$Cockroaches == 1 ~ 0,
      .$Taxonomy %in% lessSocial ~ 1,
      .$Taxonomy %in% moreSocial ~ 2
    )
  )

data$Sociality = as.factor(data$Sociality)

data = data[!is.na(data$Sociality),]

strand.labs <- c("Minor strand", "Major strand")
names(strand.labs) <- c("0", "1")

a = ggplot(data, aes(Sociality, A_fr, fill = Sociality)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Less social termites', 'More social termites')) +
  xlab('') + ylab('') + ggtitle('A fraction') +
  facet_grid(MajorStrand ~ ., labeller = labeller(MajorStrand = strand.labs)) +
  scale_fill_brewer(palette="Purples")

t = ggplot(data, aes(Sociality, T_fr, fill = Sociality)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Less social termites', 'More social termites')) +
  xlab('') + ylab('') + ggtitle('T fraction') +
  facet_grid(MajorStrand ~ ., labeller = labeller(MajorStrand = strand.labs)) + 
  scale_fill_brewer(palette="Purples")

c = ggplot(data, aes(Sociality, C_fr, fill = Sociality)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Less social termites', 'More social termites')) +
  xlab('') + ylab('') + ggtitle('C fraction') +
  facet_grid(MajorStrand ~ ., labeller = labeller(MajorStrand = strand.labs)) + 
  scale_fill_brewer(palette="Purples")

g = ggplot(data, aes(Sociality, G_fr, fill = Sociality)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Less social termites', 'More social termites')) +
  xlab('') + ylab('') + ggtitle('G fraction') +
  facet_grid(MajorStrand ~ ., labeller = labeller(MajorStrand = strand.labs)) +
  scale_fill_brewer(palette="Purples")

plots = plot_grid(a, t, g, c, nrow = 1)

save_plot('../results/nucleotide_content06_20/ATGCplotsSociality.pdf', plots, base_height = 15)
