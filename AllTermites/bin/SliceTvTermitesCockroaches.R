rm(list = ls(all=TRUE))

mut = read.table('../results/nd6_22_01/mutSpectrumFractions.txt', header = TRUE, sep='\t')

mut$TsTv = (mut$C_T + mut$A_G + mut$T_C + mut$G_A) / (mut$G_T + mut$A_T + mut$G_C +
                                                        mut$A_C + mut$C_A + mut$C_G +
                                                        mut$T_A + mut$T_G)

mut$Ts = mut$C_T + mut$A_G + mut$T_C + mut$G_A
mut$Tv = mut$G_T + mut$A_T + mut$G_C +  mut$A_C + mut$C_A + mut$C_G + mut$T_A + mut$T_G

############################################################################

workers = mut[mut$Worker == 1,]
workers = workers[!is.na(workers$Worker),]
withoutWorkers = mut[(mut$Worker == 0) & !is.na(mut$Worker),]

filter_workers = rbind(workers, withoutWorkers)

mut$Soldier = as.factor(mut$Soldier)
filter_soldiers = mut[!is.na(mut$Soldier),]

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

###############################################################################

summary(mut$Tv)

quantile(mut$Tv, probs = seq(0, 1, 1/5))[2]

nrow(mut[mut$Tv == quantile(mut$Tv, probs = seq(0, 1, 1/5))[2],])

one_line = c()
for(i in 1:5){
  # i = 3
  start = quantile(mut$Tv, probs = seq(0, 1, 1/5))[i]
  end = quantile(mut$Tv, probs = seq(0, 1, 1/5))[i + 1]
  temp_mut = mut[mut$Tv >= start & mut$Tv < end,]
  if(i == 4){
    temp_mut = mut[mut$Tv >= start & mut$Tv <= end,]
  }
  cockroaches = temp_mut[temp_mut$Cockroaches == 1,]
  termites = temp_mut[temp_mut$Cockroaches == 0,]
  result = wilcox.test(cockroaches$Ts, termites$Ts)
  one_line = rbind(one_line, c(nrow(cockroaches), nrow(termites), end, result$statistic, result$p.value))
}

wilcoxTableQuantiles = as.data.frame(one_line)
names(wilcoxTableQuantiles) = c('N_cockroaches', 'N_termites', 'End', 'W', 'Pvalue')

write.table(wilcoxTableQuantiles, '../results/nd6_22_01/wilcoxTableQuantiles.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)

##########################

vecTv = seq(0, max(mut$Tv), 100)

one_line = c()
for(i in 1:(length(vecTv) - 1)){
  # i = 1
  start = vecTv[i]; end = vecTv[i + 1]
  temp_mut = mut[mut$Tv >= start & mut$Tv < end,]
  cockroaches = temp_mut[temp_mut$Cockroaches == 1,]
  termites = temp_mut[temp_mut$Cockroaches == 0,]
  result = wilcox.test(cockroaches$Ts, termites$Tv)
  one_line = rbind(one_line, c(nrow(cockroaches), nrow(termites), result$statistic, result$p.value))
}

# not enough observations after 200 transitions
