rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

if(!require(dplyr)){install.packages('dplyr')}
if(!require(dplyr)){install.packages('ggplot2')}
if(!require(dplyr)){install.packages('cowplot')}

library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)


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

higher_termites = c("Apicotermitinae", "Cephalo-group", "Microcerotermes", 
                    "Termes-group", "Nasutitermitinae", "Amitermes-group", 
                    "Promiro", "Macrotermitinae", "Cubitermitinae", "Foraminitermitinae", 
                    "Syntermitinae", "Sphaerotermitinae", "Neocapri-group", 
                    'Pericapritermes-group', "pericapritermes-group")

data = data[!is.na(data$Taxonomy),] # remove Mantids

data$A_fr = data$A / (data$A + data$T + data$G + data$C)
data$T_fr = data$T / (data$A + data$T + data$G + data$C)
data$G_fr = data$G / (data$A + data$T + data$G + data$C)
data$C_fr = data$C / (data$A + data$T + data$G + data$C)


for(i in 1:nrow(data)){
  if(data$Taxonomy[i] %in% cockroaches){
    data$Termites[i] = 0
  }
  if(!(data$Taxonomy[i] %in% cockroaches))
  {data$Termites[i] = 1}
}

data$Termites = as.factor(data$Termites)

major_strand = data[data$MajorStrand == 1,]

minor_strand = data[data$MajorStrand == 0,]

my_colors <- RColorBrewer::brewer.pal(2, "Purples")[c(1,3)]

a_major = ggplot(major_strand, aes(Termites, A_fr, fill = Termites)) +
  geom_violin() + 
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('A in major genes')

a_major

t_major = ggplot(major_strand, aes(Termites, T_fr, fill = Termites)) +
  geom_violin() + 
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('T in major genes')

t_major

g_major = ggplot(major_strand, aes(Termites, G_fr, fill = Termites)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('G in major genes')

g_major

c_major = ggplot(major_strand, aes(Termites, C_fr, fill = Termites)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('C in major genes')

c_major

#################################
# minor genes

a_minor = ggplot(minor_strand, aes(Termites, A_fr, fill = Termites)) +
  geom_violin() + 
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('A in minor genes')

a_minor

t_minor = ggplot(minor_strand, aes(Termites, T_fr, fill = Termites)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('T in minor genes')

t_minor

g_minor = ggplot(minor_strand, aes(Termites, G_fr, fill = Termites)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('G in minor genes')

g_minor

c_minor = ggplot(minor_strand, aes(Termites, C_fr, fill = Termites)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('C in minor genes')

c_minor

plots = plot_grid(a_major, t_major, g_major, c_major, a_minor, t_minor, g_minor, c_minor, nrow = 2)

save_plot('../results/nucleotide_content06_20/ATGCplotsTermitesVsCockroaches.pdf', plots, base_height = 7)

#################################################################################
### Compare cockroaches, less and more social termites

data = data %>% 
  mutate(
    Sociality = case_when(
      .$Termites == 0 ~ 0,
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

a

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

t

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

c

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

g

plots = plot_grid(a, t, g, c, nrow = 1)

save_plot('../results/nucleotide_content06_20/ATGCplotsSociality.pdf', plots, base_height = 15)

############################################################
## merge both strands

names(minor_strand) = c('Species', 'Gene', 'T', 'A', 'C', 'G', 'MajorStrand',
                        'Taxonomy', 'T_fr', 'A_fr', 'C_fr', 'G_fr', 'Termites')

all_strands = full_join(major_strand, minor_strand)

all_nucl = all_strands %>%
  group_by(Species) %>%
  summarise(A_both = sum(A), T_both = sum(T), G_both = sum(G), C_both = sum(C))

# sum(all_strands[all_strands$Species == "B044_Panchlora_nivea_1",]$G)

all_fr = all_nucl %>%
  mutate(
    A_fr = A_both / (A_both + T_both + C_both + G_both),
    T_fr = T_both / (A_both + T_both + C_both + G_both),
    G_fr = G_both / (A_both + T_both + C_both + G_both),
    C_fr = C_both / (A_both + T_both + C_both + G_both)
  )

all_fr_families = left_join(all_fr, families)

all_fr_sociality = all_fr_families %>% 
  mutate(
    Sociality = case_when(
      .$Taxonomy %in% cockroaches ~ 0,
      .$Taxonomy %in% lessSocial ~ 1,
      .$Taxonomy %in% moreSocial ~ 2
    ),
    Termites = case_when(
      .$Taxonomy %in% cockroaches ~ 0,
      .$Taxonomy %in% c(lessSocial, moreSocial) ~ 1
    )
  )


all_fr_sociality$Sociality = as.factor(all_fr_sociality$Sociality)
all_fr_sociality$Termites = as.factor(all_fr_sociality$Termites)

all_fr_sociality = all_fr_sociality[!is.na(all_fr_sociality$Sociality),]

  
my_colors <- RColorBrewer::brewer.pal(2, "Purples")[c(1,3)]

t_all = ggplot(all_fr_sociality, aes(Sociality, A_fr, fill = Sociality)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Less\n social\n termites', 'More\n social\n termites')) +
  xlab('') + ylab('') + ggtitle('T fraction') +
  scale_fill_brewer(palette="Purples")

t_all

a_all = ggplot(all_fr_sociality, aes(Sociality, T_fr, fill = Sociality)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Less\n social\n termites', 'More\n social\n termites')) +
  xlab('') + ylab('') + ggtitle('A fraction') +
  scale_fill_brewer(palette="Purples")

a_all

g_all = ggplot(all_fr_sociality, aes(Sociality, C_fr, fill = Sociality)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Less\n social\n termites', 'More\n social\n termites')) +
  xlab('') + ylab('') + ggtitle('G fraction') +
  scale_fill_brewer(palette="Purples")

c_all = ggplot(all_fr_sociality, aes(Sociality, G_fr, fill = Sociality)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Less\n social\n termites', 'More\n social\n termites')) +
  xlab('') + ylab('') + ggtitle('C fraction') +
  scale_fill_brewer(palette="Purples")

c_all

# cockroaches

t_all_roaches = ggplot(all_fr_sociality, aes(Termites, A_fr, fill = Termites)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('T fraction')

a_all_roaches = ggplot(all_fr_sociality, aes(Termites, T_fr, fill = Termites)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('A fraction')
  

g_all_roaches = ggplot(all_fr_sociality, aes(Termites, C_fr, fill = Termites)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('G fraction')

c_all_roaches = ggplot(all_fr_sociality, aes(Termites, G_fr, fill = Termites)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Termites')) +
  scale_fill_manual(values = my_colors) +
  xlab('') + ylab('') + ggtitle('C fraction')

plots3 = plot_grid(a_all, t_all, g_all, c_all, a_all_roaches, t_all_roaches, 
                   g_all_roaches, c_all_roaches, nrow = 2)

save_plot('../results/nucleotide_content06_20/ATGCplotMergedStrands.pdf',
          plots3, base_height = 10)

write.table(all_fr_sociality, '../results/nucleotide_content06_20/ATGCplotMergedStrands.txt',
            sep='\t', quote = FALSE, row.names = FALSE)

################################################################################
### pairwise comparisons

# cockroaches vs termites 

t.test(all_fr_sociality[all_fr_sociality$Termites == 0,]$A_fr,
       all_fr_sociality[all_fr_sociality$Termites == 1,]$A_fr)

# t = -17.632, df = 109.79, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.1505364 -0.1201153
# sample estimates:
#   mean of x mean of y 
# 0.5599491 0.6952750 

t.test(all_fr_sociality[all_fr_sociality$Termites == 0,]$T_fr,
       all_fr_sociality[all_fr_sociality$Termites == 1,]$T_fr)

# t = 32.474, df = 108.19, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.1780313 0.2011776
# sample estimates:
#   mean of x mean of y 
# 0.3021161 0.1125116 

t.test(all_fr_sociality[all_fr_sociality$Termites == 0,]$G_fr,
       all_fr_sociality[all_fr_sociality$Termites == 1,]$G_fr)

# t = -14.484, df = 180.91, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.03270325 -0.02486117
# sample estimates:
#   mean of x  mean of y 
# 0.03536801 0.06415021 

t.test(all_fr_sociality[all_fr_sociality$Termites == 0,]$C_fr,
       all_fr_sociality[all_fr_sociality$Termites == 1,]$C_fr)

# t = -6.6622, df = 120.21, p-value = 8.585e-10
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.03307351 -0.01791924
# sample estimates:
#   mean of x mean of y 
# 0.1025668 0.1280632 

######### kalo vs termitidae

nrow(all_fr_sociality[all_fr_sociality$Taxonomy == 'Kalotermitidae',]) #62
nrow(all_fr_sociality[all_fr_sociality$Taxonomy %in% higher_termites,]) #272

t.test(all_fr_sociality[all_fr_sociality$Taxonomy == 'Kalotermitidae',]$A_fr,
       all_fr_sociality[all_fr_sociality$Taxonomy %in% higher_termites,]$A_fr)

# t = -5.0685, df = 69.964, p-value = 3.138e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.04811842 -0.02094295
# sample estimates:
#   mean of x mean of y 
# 0.6649318 0.6994625 

t.test(all_fr_sociality[all_fr_sociality$Taxonomy == 'Kalotermitidae',]$T_fr,
       all_fr_sociality[all_fr_sociality$Taxonomy %in% higher_termites,]$T_fr)

# t = 5.2714, df = 95.34, p-value = 8.416e-07
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.01089211 0.02405170
# sample estimates:
#   mean of x mean of y 
# 0.1307202 0.1132483 

t.test(all_fr_sociality[all_fr_sociality$Taxonomy == 'Kalotermitidae',]$G_fr,
       all_fr_sociality[all_fr_sociality$Taxonomy %in% higher_termites,]$G_fr)

# t = 0.94508, df = 75.503, p-value = 0.3476
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.003768564  0.010573215
# sample estimates:
#   mean of x  mean of y 
# 0.06667796 0.06327563 

t.test(all_fr_sociality[all_fr_sociality$Taxonomy == 'Kalotermitidae',]$C_fr,
       all_fr_sociality[all_fr_sociality$Taxonomy %in% higher_termites,]$C_fr)

# t = 3.1812, df = 75.103, p-value = 0.002133
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.005104899 0.022208020
# sample estimates:
#   mean of x mean of y 
# 0.1376700 0.1240136 


# add lower and higher termites

lowerTermites = all_fr_sociality[all_fr_sociality$Termites == 1 & !(all_fr_sociality$Taxonomy %in% higher_termites),]
higherTermites = all_fr_sociality[all_fr_sociality$Taxonomy %in% higher_termites,]

# add Stylotermitidae to lower termites

t.test(lowerTermites$A_fr, higherTermites$A_fr)

# t = -2.7213, df = 167.49, p-value = 0.00719
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.022579745 -0.003592539
# sample estimates:
#   mean of x mean of y 
# 0.6863764 0.6994625 

t.test(lowerTermites$T_fr, higherTermites$T_fr)

# t = -0.73293, df = 205.4, p-value = 0.4644
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.008494537  0.003890433
# sample estimates:
#   mean of x mean of y 
# 0.1109462 0.1132483 

t.test(lowerTermites$G_fr, higherTermites$G_fr)

# t = 1.1012, df = 199.35, p-value = 0.2721
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.002161036  0.007627186
# sample estimates:
#   mean of x  mean of y 
# 0.06600871 0.06327563 

t.test(lowerTermites$C_fr, higherTermites$C_fr)

# t = 4.6339, df = 213.29, p-value = 6.243e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   0.007271959 0.018038279
# sample estimates:
#   mean of x mean of y 
# 0.1366687 0.1240136 

all_fr_higher = all_fr_families %>%
  mutate(
    Higher = case_when(.$Taxonomy %in% higher_termites ~ 2,
                       .$Taxonomy %in% c(unique(as.character(lowerTermites$Taxonomy)), 
                                         'Stylotermitidae') ~ 1,
                       .$Taxonomy %in% cockroaches ~ 0)
  )

all_fr_higher$Higher = as.factor(all_fr_higher$Higher)

t_all_higher = ggplot(all_fr_higher, aes(Higher, A_fr, fill = Higher)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Lower\n termites', 'Higher\n termites')) +
  xlab('') + ylab('') + ggtitle('T fraction') +
  scale_fill_brewer(palette="Purples")

t_all_higher

# all_fr_higher[is.na(all_fr_higher$Higher),]

summary(all_fr_higher[all_fr_higher$Taxonomy == 'Coptotermitinae',]$A_fr)
summary(all_fr_higher[all_fr_higher$Taxonomy == 'Rhinotermitinae',]$A_fr)
summary(all_fr_higher[all_fr_higher$Taxonomy == 'Kalotermitidae',]$A_fr)

a_all_higher = ggplot(all_fr_higher, aes(Higher, T_fr, fill = Higher)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Lower\n termites', 'Higher\n termites')) +
  xlab('') + ylab('') + ggtitle('A fraction') +
  scale_fill_brewer(palette="Purples")

a_all_higher


g_all_higher = ggplot(all_fr_higher, aes(Higher, C_fr, fill = Higher)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Lower\n termites', 'Higher\n termites')) +
  xlab('') + ylab('') + ggtitle('G fraction') +
  scale_fill_brewer(palette="Purples")

g_all_higher

c_all_higher = ggplot(all_fr_higher, aes(Higher, G_fr, fill = Higher)) +
  geom_violin() + 
  # scale_fill_brewer(palette="Set3") +
  theme_minimal() + 
  stat_summary(fun.y = 'median', geom = 'point', shape = 20, size = 3, col = 'midnightblue') +
  theme(axis.text.x = element_text(color = "black", size = 12), 
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none") + 
  scale_x_discrete(labels=c('Cockroaches', 'Lower\n termites', 'Higher\n termites')) +
  xlab('') + ylab('') + ggtitle('C fraction') +
  scale_fill_brewer(palette="Purples")

c_all_higher


plots4 = plot_grid(a_all_higher, t_all_higher, g_all_higher, c_all_higher, nrow = 2)

save_plot('../results/nucleotide_content06_20/ATGCplotHigherLowerTermites.pdf',
          plots4, base_height = 10)
