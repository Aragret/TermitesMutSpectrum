library(ggplot2)

mut = read.table('../../Body/3Results/Mutational_spectra_in_Termites.txt', header=TRUE)

mut = mut[which(!(mut$Subs %in% c('A_N', 'A_R', 'C_N', 'G_N', 'T_K', 'T_N', 'T_Y'))),]

mut$Subs = as.character(mut$Subs)

# table(mut$Subs)

taxons = read.csv('../../Body/2Derived/TermitesFamilies.csv')
data = merge(mut, taxons, by='Species', all.x = TRUE)

mutFraction = function(data){
  a = as.data.frame(table(data$Subs))
  all_mut = sum(a$Freq)
  names(a) = c('Subs', 'Number')
  a$SubsFreq = a$Number / all_mut
  return(a)
}

higher = mutFraction(data[data$Family == 'Termitidae' | data$Family == 'Rhinotermitidae',])
lower = mutFraction(data[data$Family != 'Termitidae' & data$Family != 'Rhinotermitidae',])

ggplot(higher, aes(Subs, SubsFreq)) +
  geom_bar(stat = "identity")

ggplot(lower, aes(Subs, SubsFreq)) +
  geom_bar(stat = "identity")

# I see less G>A in "lower" termites

# Try clusterisation ?
