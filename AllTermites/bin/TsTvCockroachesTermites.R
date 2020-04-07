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
summary(mut$Cockroaches)

mut$HigherTermites = as.numeric(mut$HigherTermites)
mut$Cockroaches = as.numeric(mut$Cockroaches)

mut$TermitesCockroaches = mut$HigherTermites + 1
mut$TermitesCockroaches = mut$TermitesCockroaches - mut$Cockroaches
mut$TermitesCockroaches = as.factor(mut$TermitesCockroaches)

summary(mut$TermitesCockroaches)

mut$HigherTermites = as.factor(as.numeric(mut$HigherTermites) - 1)
mut$Cockroaches = as.factor(as.numeric(mut$Cockroaches) - 1)


finiteMut = mut[mut$TsTv != 'Inf',]


# termites vs cockroaches

result1 <- aov(Ts ~ Tv * TermitesCockroaches, data = finiteMut)
print(summary(result1))

result2 <- aov(Ts ~ Tv + TermitesCockroaches, data = finiteMut)
print(summary(result2))

print(anova(result1,result2))


# cockroaches vs higher termites vs lower termites

result1 <- aov(Ts ~ Tv * Cockroaches, data = finiteMut)
print(summary(result1))

result2 <- aov(Ts ~ Tv + Cockroaches, data = finiteMut)
print(summary(result2))

print(anova(result1,result2))

#########################################################################

pdf('../results/nd6_22_01/TsTvCockroachesTermites.R.pdf')

hist(mut$Ts, breaks = 50)
hist(mut$Tv, breaks = 50)

Result <- glm(Ts ~ Tv + Cockroaches, data = mut, family = poisson)
summary(Result)

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   4.592e+00  6.928e-03  662.81   <2e-16 ***
#   Tv            2.336e-03  5.122e-05   45.60   <2e-16 ***
#   Cockroaches1 -4.308e-01  1.202e-02  -35.84   <2e-16 ***

Result <- glm(Ts ~ Tv * Cockroaches, data = mut, family = poisson)
summary(Result)

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      4.536e+00  7.641e-03 593.622  < 2e-16 ***
#   Tv               2.870e-03  5.826e-05  49.257  < 2e-16 ***
#   Cockroaches1    -9.446e-02  2.154e-02  -4.386 1.16e-05 ***
#   Tv:Cockroaches1 -2.095e-03  1.165e-04 -17.975  < 2e-16 ***

Result <- glm(Ts ~ Tv + Cockroaches + BranchLength, data = mut, family = poisson)
summary(Result)

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   4.592e+00  6.969e-03 659.013   <2e-16 ***
#   Tv            2.364e-03  5.921e-05  39.924   <2e-16 ***
#   Cockroaches1 -4.285e-01  1.225e-02 -34.980   <2e-16 ***
#   BranchLength -5.509e-02  5.822e-02  -0.946    0.344    

Result <- glm(Ts ~ (Tv + Cockroaches) * BranchLength, data = mut, family = poisson)
summary(Result)

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                4.455e+00  1.199e-02 371.509   <2e-16 ***
#   Tv                         2.376e-03  7.844e-05  30.294   <2e-16 ***
#   Cockroaches1              -3.919e-01  2.205e-02 -17.775   <2e-16 ***
#   BranchLength               3.101e+00  2.251e-01  13.778   <2e-16 ***
#   Tv:BranchLength           -8.118e-03  8.570e-04  -9.473   <2e-16 ***
#   Cockroaches1:BranchLength -4.691e-01  1.874e-01  -2.503   0.0123 *  


dev.off()
