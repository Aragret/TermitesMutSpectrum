rm(list = ls(all=TRUE))

mut = read.table('../results/nd6_22_01/mutSpectrumFractions.txt', header = TRUE, sep='\t')

mut$TsTv = (mut$C_T + mut$A_G + mut$T_C + mut$G_A) / (mut$G_T + mut$A_T + mut$G_C +
                                                        mut$A_C + mut$C_A + mut$C_G +
                                                        mut$T_A + mut$T_G)

mut$Ts = mut$C_T + mut$A_G + mut$T_C + mut$G_A
mut$Tv = mut$G_T + mut$A_T + mut$G_C +  mut$A_C + mut$C_A + mut$C_G + mut$T_A + mut$T_G

############################################################################

# families of cockroaches 
cockroaches = c('Ectobiidae1', 'Tryonicidae', 'Blaberidae', 'Corydiidae', 'Ectobiidae2',
                'Lamproblattidae', 'Anaplectidae', 'Blattidae', 'Cryptocercidae', 'Ectobiidae3',
                'Nocticolidae')

mut = mut[!is.na(mut$Taxonomy),] # remove Mantids

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

#########################################################################

pdf('../results/nd6_22_01/TsTvCockroachesTermites.R.pdf')

hist(mut$Ts, breaks = 50)
hist(mut$Tv, breaks = 50)


plot(mut$Tv, mut$Ts)
plot(mut$Tv, log(mut$Ts))
plot(log(mut$Tv), mut$Ts)
plot(log(mut$Tv), log(mut$Ts))


### residuals from models of Ts vs Tv

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
# 1 paragraph of results

Result <- glm(Tv ~ Ts + Cockroaches, data = mut, family = poisson)
summary(Result)

# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  3.433e+00  1.437e-02  238.94   <2e-16 ***
#   Ts           8.171e-03  9.522e-05   85.81   <2e-16 ***
#   Cockroaches1 8.225e-01  9.905e-03   83.04   <2e-16 ***



dev.off()
