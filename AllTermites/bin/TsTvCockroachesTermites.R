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

summary(mut$Cockroaches)

finiteMut = mut[mut$TsTv != 'Inf',]


# termites vs cockroaches

result1 <- aov(Ts ~ Tv * Cockroaches, data = finiteMut)
print(summary(result1))

result2 <- aov(Ts ~ Tv + Cockroaches, data = finiteMut)
print(summary(result2))

print(anova(result1,result2))

#########################################################################

pdf('../results/nd6_22_01/TsTvCockroachesTermites.R.pdf')

hist(mut$Ts, breaks = 50)
hist(mut$Tv, breaks = 50)


plot(mut$Tv, mut$Ts)
plot(mut$Tv, log(mut$Ts))
plot(log(mut$Tv), mut$Ts)
plot(log(mut$Tv), log(mut$Ts))


### residuals

Lm1 <- lm(Ts ~ Tv, mut)

mut$res <- residuals(Lm1)

hist(mut$res, main = 'Ts ~ Tv res')
plot(mut$Tv, mut$res, main = 'Ts ~ Tv res')


Lm1 <- lm(log(Ts) ~ Tv, mut)

mut$resLog <- residuals(Lm1)

hist(mut$resLog, breaks = 20, main='log(Ts) ~ Tv res')
plot(mut$Tv, mut$resLog, main='log(Ts) ~ Tv res')


Lm1 <- lm(log(Tv) ~ Ts, finiteMut)

finiteMut$resLog <- residuals(Lm1)

hist(finiteMut$resLog, breaks = 20, main='log(Tv) ~ Ts res')
plot(finiteMut$Ts, finiteMut$resLog, main='log(Ts) ~ Tv res')
plot(log(finiteMut$Ts), finiteMut$resLog, main='log(Ts) ~ Tv res')


#######################################################################

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
