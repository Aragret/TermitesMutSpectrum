library(dplyr)

mutCock = read.table('~/Documents/lab/TermitesMutSpectrum/Body/3Results/Cockroaches.MutSpecData.txt', header = TRUE)
mutTer = read.table('~/Documents/lab/TermitesMutSpectrum/Body/3Results/Termites.MutSpecData.txt', header = TRUE)

# mann-whitney test for normilized mutspec
test_prep <- function(mut) {
  normMut = mut %>%
    group_by(Species) %>% 
    mutate(
      AllSubs = length(Subs)
    ) %>%
    filter(Subs == 'T_C') %>% 
    mutate(
      T_C = length(Subs) / AllSubs
    ) %>% 
    distinct(Species, .keep_all = TRUE) %>%
    mutate(
      T_fr = T / (A + T + G + C),
      T_C_norm = T_C / T_fr
    )
}
normMutCock <- test_prep(mutCock)
normMutTer <- test_prep(mutTer)
t.test(normMutTer$T_C_norm, normMutCock$T_C_norm)



# t = 7.621, df = 40.123, p-value = 2.549e-09
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# 2.200041 3.787887
# sample estimates:
#   mean of x mean of y 
# 4.013784  1.019820  


wilcox.test(normMutTer$T_C_norm, normMutCock$T_C_norm,
            conf.int = TRUE)

# W = 642, p-value = 6.261e-11
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   1.811509 3.183974
# sample estimates:
#   difference in location 
# 2.490454  

# aov

aov_prep <- data.frame(T_C_norm = normMutTer[normMutTer$Family == 'Termitidae',]$T_C_norm, Group = "Higher Ter")
aov_prep <- rbind(aov_prep, data.frame(T_C_norm = normMutTer[normMutTer$Family != 'Termitidae',]$T_C_norm, Group = "Lower Ter"))
aov_prep <- rbind(aov_prep, data.frame(T_C_norm = normMutCock$T_C_norm, Group = "Cockroach"))

aov.res <- aov(T_C_norm ~ Group, data = aov_prep) ## Google about it

summary(aov.res)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# Group        2  153.3   76.67   30.25 2.18e-09 ***
# Residuals   51  129.3    2.53                     
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov.res)
#                        diff        lwr        upr     p adj
# Higher Ter-Cockroach  3.953346  2.7220573  5.1846347 0.0000000
# Lower Ter-Cockroach   1.714788  0.3872827  3.0422942 0.0082643
# Lower Ter-Higher Ter -2.238558 -3.5513401 -0.9257751 0.0004082