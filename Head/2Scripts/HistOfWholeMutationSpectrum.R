library(ggplot2)

mut = read.table('../../Body/3Results/Mutational_spectra_in_Termites.txt', header=TRUE)
taxons = read.csv('../../Body/2Derived/TermitesFamilies.csv')

table(mut$Subs)

length(table(mut[as.character(mut$AncestralAA) == as.character(mut$DescendantAA),]$Subs))

mut$Subs = as.character(mut$Subs)

mut = mut[which(!(mut$Subs %in% c('A_N', 'A_R', 'C_N', 'G_N', 'T_K', 'T_N', 'T_Y'))),]

data = merge(mut, taxons, by='Species', all.x = TRUE)

plotHist <- function(data, taxon) {
  temp_data = data[data$Family == taxon,]
  graph = ggplot(temp_data, aes(x = Subs)) +
    geom_histogram(stat = 'count') +
    ggtitle(taxon)
  return(graph)
}

pdf('../../Body/4Figures/HistOfWholeMutationSpectrum.R.pdf')

ggplot(data = mut, aes(x = Subs)) +
  geom_histogram(stat = 'count')

for(taxon in unique(as.character(data$Family))){
  print(plotHist(data, taxon))
}

dev.off()

nrow(mut[mut$Subs == 'C_T' | mut$Subs == 'G_A',]) # 904
nrow(mut[mut$Subs == 'T_C' | mut$Subs == 'A_G',]) # 756
