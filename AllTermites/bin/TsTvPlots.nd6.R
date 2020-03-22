rm(list = ls(all=TRUE))

library(ggplot2)

mut = read.table('../results/nd6_22_01/mutSpectrumFractions.txt', header = TRUE, sep='\t')

mut$TsTv = (mut$C_T + mut$A_G + mut$T_C + mut$G_A) / (mut$G_T + mut$A_T + mut$G_C +
                                                        mut$A_C + mut$C_A + mut$C_G +
                                                        mut$T_A + mut$T_G)

mut$Ts = mut$C_T + mut$A_G + mut$T_C + mut$G_A
mut$Tv = mut$G_T + mut$A_T + mut$G_C +  mut$A_C + mut$C_A + mut$C_G + mut$T_A + mut$T_G


##################################################################################
### Workers

pdf('../results/nd6_22_01/TsTvPlots.nd6.pdf')

workers = mut[mut$Worker == 1,]
workers = workers[!is.na(workers$Worker),]
withoutWorkers = mut[(mut$Worker == 0) & !is.na(mut$Worker),]

filter_workers = rbind(workers, withoutWorkers)

summary(filter_workers$TsTv)

ggplot(filter_workers, aes(TsTv, fill = filter_workers$Worker)) +
  # geom_histogram(aes(fill = filter_workers$Worker), alpha = 0.4) +
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'dodge') +
  scale_fill_manual(values=c("#404080", "#69b3a2"))

mut$Soldier = as.factor(mut$Soldier)
filter_soldiers = mut[!is.na(mut$Soldier),]

ggplot(filter_soldiers, aes(TsTv, fill = filter_soldiers$Soldier)) +
  # geom_histogram(aes(fill = filter_workers$Worker), alpha = 0.4) +
  geom_histogram( color="#e9ecef", alpha=0.4, position = 'dodge') +
  scale_fill_manual(values=c("#404080", "#69b3a2"))

##################################################################################

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

plot(mut$sumOfSubs, mut$TsTv)
plot(mut$sumOfSubs, mut$Ts)
plot(mut$sumOfSubs, mut$Tv)
plot(mut$Ts, mut$Tv)

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

mut$Cockroaches = as.factor(mut$Cockroaches)

ggplot(mut, aes(Ts, Tv, col = Cockroaches)) + 
  geom_point() +
  geom_abline()

ggplot(mut, aes(Ts, Tv, col = Cockroaches)) + 
  geom_point(aes(size = HigherTermites)) +
  geom_abline()

#################################################################################

ggplot(mut[mut$Cockroaches == 0,], aes(Ts, Tv)) +
  geom_point() +
  geom_abline() +
  xlim(min(mut$Ts), max(mut$Ts)) +
  ylim(min(mut$Tv), max(mut$Tv)) +
  facet_wrap(~ Taxonomy)

ggplot(mut[mut$Cockroaches == 1,], aes(Ts, Tv)) +
  geom_point() +
  geom_abline() + xlim(min(mut$Ts), max(mut$Ts)) +
  ylim(min(mut$Tv), max(mut$Tv)) +
  facet_wrap(~ Taxonomy)

#####################################################################
mut$HigherTermites = as.numeric(mut$HigherTermites)
mut$Cockroaches = as.numeric(mut$Cockroaches)

mut$TermitesCockroaches = mut$HigherTermites + 1
mut$TermitesCockroaches = mut$TermitesCockroaches - mut$Cockroaches
mut$TermitesCockroaches = as.factor(mut$TermitesCockroaches)

summary(mut$TermitesCockroaches)

ggplot(mut, aes(x=TermitesCockroaches, y=TsTv, fill=mut$TermitesCockroaches)) +
  geom_violin() +
  #scale_fill_discrete(name = "", labels = c("Cockroaches (n=98)", "Lower termites (n=130)", 'Higher termites (n=272)')) +
  xlab('') + stat_summary(fun.y="median", geom="point",size=2) +
  scale_fill_grey(name = "", labels = c("Cockroaches (n=98)", "Lower termites (n=130)", 'Higher termites (n=272)'),
                  start = 1, end = 0.3)
  

dev.off()


#####################################################################################
mut$HigherTermites = as.factor(mut$HigherTermites)
mut$Cockroaches = as.factor(mut$Cockroaches)

wilcox.test(mut[mut$TermitesCockroaches == 0,]$TsTv, mut[mut$TermitesCockroaches == 1,]$TsTv)
# W = 3167, p-value = 1.852e-10

wilcox.test(mut[mut$TermitesCockroaches == 1,]$TsTv, mut[mut$TermitesCockroaches == 2,]$TsTv)
# W = 15564, p-value = 0.08751

wilcox.test(mut[mut$TermitesCockroaches == 0,]$TsTv, mut[mut$TermitesCockroaches == 2,]$TsTv)
# W = 5408, p-value < 2.2e-16


mut$HigherTermites = as.factor(as.numeric(mut$HigherTermites) - 1)
mut$Cockroaches = as.factor(as.numeric(mut$Cockroaches) - 1)

finiteMut = mut[mut$TsTv != 'Inf',]
summary(lm(scale(TsTv) ~ HigherTermites + Cockroaches + scale(BranchLength), data = finiteMut))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          0.14592    0.08252   1.768   0.0776 .  
# HigherTermites1     -0.07631    0.10072  -0.758   0.4490    
# Cockroaches1        -0.52973    0.12877  -4.114 4.56e-05 ***
#   scale(BranchLength) -0.25638    0.04469  -5.736 1.69e-08 ***

pdf('../results/nd6_22_01/TsTvPlotsBranchLength.nd6.pdf')

ggplot(finiteMut, aes(y = TsTv, x = BranchLength, col = Cockroaches)) +
  geom_point()

ggplot(finiteMut, aes(y = TsTv, x = sumOfSubs, col = Cockroaches)) +
  geom_point()

dev.off()
