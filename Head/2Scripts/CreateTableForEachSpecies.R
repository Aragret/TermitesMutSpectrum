mut = read.table('../../Body/3Results/Mutational_spectra_in_Termites.txt', header=TRUE)

mut = mut[which(!(mut$Subs %in% c('A_N', 'A_R', 'C_N', 'G_N', 'T_K', 'T_N', 'T_Y'))),]

mut$Subs = as.character(mut$Subs)

mutTable = setNames(data.frame(matrix(0, ncol = 12, nrow = 37)), unique(mut$Subs))
Species = unique(mut$Species)
mutTable = cbind(Species, mutTable)


for(i in 1:nrow(mut)){
  row = mut[i, ]
  mutTable[mutTable$Species == row$Species, row$Subs] = mutTable[mutTable$Species == row$Species, row$Subs] + 1
}

table(mut[mut$Species == "Macrotermes_gilvus",]$Subs)

write.table(mutTable, '../../Body/2Derived/MutNumbersForEachSpecies.txt', sep='\t', 
            row.names = FALSE, quote = FALSE)
