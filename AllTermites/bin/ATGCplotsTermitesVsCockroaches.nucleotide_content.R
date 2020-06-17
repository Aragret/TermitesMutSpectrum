rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

library(ggplot2)
library(cowplot)

cont = read.table('../results/nucleotide_content06_20/ATGCforEachGene4fold.txt', 
                  header = TRUE, sep='\t')
cont$Species = sub(';', '', cont$Species)
cont$Species = unlist(lapply(cont$Species, paste, '_1', sep=''))
cont$Species = sub('__', '_', cont$Species)
cont$Species = sub('?', '', cont$Species, fixed = TRUE)


mut = read.table('../results/nd6_22_01/mutSpectrumFractions.txt', header = TRUE, sep='\t')
families = mut[, c('Species', 'Taxonomy')]

data = merge(cont, families, by='Species')

setdiff(families$Species, a$Species)

boxplot(data$A, data$T, data$G, data$C)

cockroaches = c('Ectobiidae1', 'Tryonicidae', 'Blaberidae', 'Corydiidae', 'Ectobiidae2',
                'Lamproblattidae', 'Anaplectidae', 'Blattidae', 'Cryptocercidae', 'Ectobiidae3',
                'Nocticolidae')

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
  # scale_fill_brewer(palette="Set3") +
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
