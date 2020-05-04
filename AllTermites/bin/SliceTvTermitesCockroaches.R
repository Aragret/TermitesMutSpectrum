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

finiteMut = mut[mut$TsTv != 'Inf',]

###############################################################################

summary(mut$Ts)

mut = mut[mut$Ts <= 150,] # not enough observations for cockroaches after 150 Ts

quantile(mut$Ts, probs = seq(0, 1, 1/5)) # 5 segments of transitions

one_line = c()
for(i in 1:5){
  # i = 3
  start = quantile(mut$Ts, probs = seq(0, 1, 1/5))[i]
  end = quantile(mut$Ts, probs = seq(0, 1, 1/5))[i + 1]
  temp_mut = mut[mut$Ts >= start & mut$Ts < end,]
  if(i == 5){
    temp_mut = mut[mut$Ts >= start & mut$Ts <= end,]
  }
  cockroaches = temp_mut[temp_mut$Cockroaches == 1,]
  termites = temp_mut[temp_mut$Cockroaches == 0,]
  cockroachesTv = median(cockroaches$Tv)
  termitesTv = median(termites$Tv)
  result = t.test(cockroaches$Tv, termites$Tv)
  one_line = rbind(one_line, c(nrow(cockroaches), nrow(termites), end, cockroachesTv, termitesTv, result$statistic, result$p.value))
}

TtestTableQuantiles = as.data.frame(one_line)
names(TtestTableQuantiles) = c('N_cockroaches', 'N_termites', 'End', 'medianTvCockroaches', 'medianTvTermites', 't', 'Pvalue')

write.table(TtestTableQuantiles, '../results/nd6_22_01/tTestTableQuantiles.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)


one_line = c()
for(i in 1:5){
  # i = 3
  start = quantile(finiteMut$Ts, probs = seq(0, 1, 1/5))[i]
  end = quantile(finiteMut$Ts, probs = seq(0, 1, 1/5))[i + 1]
  temp_mut = finiteMut[finiteMut$Ts >= start & finiteMut$Ts < end,]
  if(i == 5){
    temp_mut = finiteMut[finiteMut$Ts >= start & finiteMut$Ts <= end,]
  }
  cockroaches = temp_mut[temp_mut$Cockroaches == 1,]
  termites = temp_mut[temp_mut$Cockroaches == 0,]
  cockroachesTsTv = median(cockroaches$TsTv)
  termitesTsTv = median(termites$TsTv)
  result = t.test(cockroaches$TsTv, termites$TsTv)
  one_line = rbind(one_line, c(nrow(cockroaches), nrow(termites), end, cockroachesTsTv, termitesTsTv, result$statistic, result$p.value))
}

tTestTableQuantilesTsTv = as.data.frame(one_line)
names(tTestTableQuantilesTsTv) = c('N_cockroaches', 'N_termites', 'End', 'medianTsTvCockroaches', 'medianTsTvTermites', 't', 'Pvalue')

write.table(tTestTableQuantilesTsTv, '../results/nd6_22_01/tTestTableQuantilesTsTv.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)
