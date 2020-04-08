rm(list = ls(all=TRUE))

mut = read.table('../results/nd6_22_01/mutSpectrumFractions.txt', header = TRUE, sep='\t')

mut$TsTv = (mut$C_T + mut$A_G + mut$T_C + mut$G_A) / (mut$G_T + mut$A_T + mut$G_C +
                                                        mut$A_C + mut$C_A + mut$C_G +
                                                        mut$T_A + mut$T_G)

mut$Ts = mut$C_T + mut$A_G + mut$T_C + mut$G_A
mut$Tv = mut$G_T + mut$A_T + mut$G_C +  mut$A_C + mut$C_A + mut$C_G + mut$T_A + mut$T_G

############################################################################

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

summary(mut$Ts)

quantile(mut$Ts, probs = seq(0, 1, 1/5))

nrow(mut[mut$Ts == quantile(mut$Ts, probs = seq(0, 1, 1/5))[2],])

one_line = c()
for(i in 1:5){
  # i = 3
  start = quantile(mut$Ts, probs = seq(0, 1, 1/5))[i]
  end = quantile(mut$Ts, probs = seq(0, 1, 1/5))[i + 1]
  temp_mut = mut[mut$Ts >= start & mut$Ts < end,]
  if(i == 4){
    temp_mut = mut[mut$Ts >= start & mut$Ts <= end,]
  }
  cockroaches = temp_mut[temp_mut$Cockroaches == 1,]
  termites = temp_mut[temp_mut$Cockroaches == 0,]
  cockroachesTv = median(cockroaches$Tv)
  termitesTv = median(termites$Tv)
  result = wilcox.test(cockroaches$Tv, termites$Tv)
  one_line = rbind(one_line, c(nrow(cockroaches), nrow(termites), end, cockroachesTv, termitesTv, result$statistic, result$p.value))
}

wilcoxTableQuantiles = as.data.frame(one_line)
names(wilcoxTableQuantiles) = c('N_cockroaches', 'N_termites', 'End', 'medianTvCockroaches', 'medianTvTermites', 'W', 'Pvalue')

write.table(wilcoxTableQuantiles, '../results/nd6_22_01/wilcoxTableQuantiles.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)

##########################

vecTs = seq(0, max(mut$Ts), 50)

one_line = c()
for(i in 1:(length(vecTs) - 1)){
  # i = 1
  start = vecTs[i]; end = vecTs[i + 1]
  temp_mut = mut[mut$Ts >= start & mut$Ts < end,]
  if(i == (length(vecTs) - 1)){
    end = max(mut$Ts)
    temp_mut = mut[mut$Ts >= start & mut$Ts <= end,]
  }
  cockroaches = temp_mut[temp_mut$Cockroaches == 1,]
  termites = temp_mut[temp_mut$Cockroaches == 0,]
  cockroachesTv = median(cockroaches$Tv)
  termitesTv = median(termites$Tv)
  result = wilcox.test(cockroaches$Tv, termites$Tv)
  one_line = rbind(one_line, c(nrow(cockroaches), nrow(termites), end, cockroachesTv, termitesTv, result$statistic, result$p.value))
}

wilcoxTableSegments = as.data.frame(one_line)
names(wilcoxTableSegments) = c('N_cockroaches', 'N_termites', 'End', 'medianTvCockroaches', 'medianTvTermites', 'W', 'Pvalue')

write.table(wilcoxTableSegments, '../results/nd6_22_01/wilcoxTableSegments.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)

# not enough observations after 150 transitions
