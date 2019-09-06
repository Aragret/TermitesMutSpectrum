library(ggplot2)

mut = read.table('../../Body/3Results/Mutational_spectra_in_Termites.txt', header=TRUE)

table(mut$Subs)

length(table(mut[as.character(mut$AncestralAA) == as.character(mut$DescendantAA),]$Subs))

mut$Subs = as.character(mut$Subs)

mut = mut[which(!(mut$Subs %in% c('A_N', 'A_R', 'C_N', 'G_N', 'T_K', 'T_N', 'T_Y'))),]

ggplot(data = mut, aes(x = Subs)) +
  geom_histogram(stat = 'count')

