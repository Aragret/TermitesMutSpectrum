rm(list = ls(all=TRUE))

library(phytools)
library(ggplot2)

mut = read.table('../results/nd6_22_01/mutSpectrumFractions.txt', header = TRUE, sep='\t')
tree = read.tree('../results/nd6_22_01/phylogeny_rooted/iqtree_mito.treefile')

mut$TsTv = (mut$C_T + mut$A_G + mut$T_C + mut$G_A) / (mut$G_T + mut$A_T + mut$G_C +
                                                        mut$A_C + mut$C_A + mut$C_G +
                                                        mut$T_A + mut$T_G)

mut$Ts = mut$C_T + mut$A_G + mut$T_C + mut$G_A
mut$Tv = mut$G_T + mut$A_T + mut$G_C +  mut$A_C + mut$C_A + mut$C_G + mut$T_A + mut$T_G

summary(mut$TsTv)

###################################################################################

higher_termites = c("Apicotermitinae", "Cephalo-group", "Microcerotermes", 
                    "Termes-group", "Nasutitermitinae", "Amitermes-group", 
                    "Promiro", "Macrotermitinae", "Cubitermitinae", "Foraminitermitinae", 
                    "Syntermitinae", "Sphaerotermitinae", "Neocapri-group", 
                    'Pericapritermes-group', "pericapritermes-group")

for(i in 1:nrow(mut)){
  # i = 1
  if(mut[i, 'Taxonomy'] %in% higher_termites){
    mut[i, 'HigherTermites'] = 1
  }
  else{mut[i, 'HigherTermites'] = 0}
}

mut$HigherTermites = as.factor(mut$HigherTermites)

cockroaches = c('Ectobiidae1', 'Tryonicidae', 'Blaberidae', 'Corydiidae', 'Ectobiidae2',
                'Lamproblattidae', 'Anaplectidae', 'Blattidae', 'Cryptocercidae', 'Ectobiidae3',
                'Nocticolidae')

mut = mut[!is.na(mut$Taxonomy),]

for(i in 1:nrow(mut)){
  if(mut$Taxonomy[i] %in% cockroaches){
    mut$Cockroaches[i] = 1
  }
  if(!(mut$Taxonomy[i] %in% cockroaches))
  {mut$Cockroaches[i] = 0}
}

mut$Cockroaches = as.factor(mut$Cockroaches)

summary(mut$Cockroaches)

ggplot(mut, aes(TsTv, fill = mut$Cockroaches)) +
  # geom_histogram(aes(fill = filter_workers$Worker), alpha = 0.4) +
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'dodge') +
  scale_fill_manual(values=c("#404080", "#69b3a2"))

#################################################################################

data = mut[which(as.character(mut$Species) %in% tree$tip.label),]
df_vec <- as.character(data$Species)
tree_vec <- tree$tip.label
a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)
row.names(data) = data$Species
tree2 <- drop.tip(tree, b)


phylANOVA(tree2, data$Cockroaches, data$TsTv)
