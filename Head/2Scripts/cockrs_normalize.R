library(dplyr)

mutCock = read.table('~/Documents/lab/TermitesMutSpectrum/Body/3Results/Cockroaches.MutSpecData.txt', header = TRUE)
mutTer = read.table('~/Documents/lab/TermitesMutSpectrum/Body/3Results/Termites.MutSpecData.txt', header = TRUE)
# External Internal 
# 662     1013

# mann-whitney test for normilized mutspec
test_prep <- function(mut) {
  normMut = mut %>%
    group_by(Species) %>% 
    mutate(
      AllSubs = length(Subs)
    ) %>% 
    mutate(
      A_T = sum(Subs == "A_T") / AllSubs,
      A_G = sum(Subs == "A_G") / AllSubs,
      A_C = sum(Subs == "A_C") / AllSubs,
      T_A = sum(Subs == "T_A") / AllSubs,
      T_G = sum(Subs == "T_G") / AllSubs,
      T_C = sum(Subs == "T_C") / AllSubs,
      G_T = sum(Subs == "G_T") / AllSubs,
      G_A = sum(Subs == "G_A") / AllSubs,
      G_C = sum(Subs == "G_C") / AllSubs,
      C_T = sum(Subs == "C_T") / AllSubs,
      C_G = sum(Subs == "C_G") / AllSubs,
      C_A = sum(Subs == "C_A") / AllSubs,
    ) %>% 
    distinct(Species, .keep_all = TRUE) %>%
    mutate(
      A_fr = A / (A + T + G + C),
      T_fr = T / (A + T + G + C),
      G_fr = G / (A + T + G + C),
      C_fr = C / (A + T + G + C),
      T_A_norm = T_A / T_fr,
      T_C_norm = T_C / T_fr,
      T_G_norm = T_G / T_fr,
      A_T_norm = A_T / A_fr,
      A_C_norm = A_C / A_fr,
      A_G_norm = A_G / A_fr,
      G_A_norm = G_A / G_fr,
      G_C_norm = G_C / G_fr,
      G_T_norm = G_T / G_fr,
      C_A_norm = C_A / C_fr,
      C_T_norm = C_T / C_fr,
      C_G_norm = C_G / C_fr,
      # Geoscapheus_robustus has 0 G_fr
      T_C_norm_fr = T_C_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm)
    )
}

normMutCock <- test_prep(mutCock)
normMutTer <- test_prep(mutTer)
normMutCock[6, "T_C_norm_fr"] = 0.123345141 #counted manually; looks like that's a seq, not whole genome!!!!!!! delete?
t.test(normMutTer$T_C_norm_fr, normMutCock$T_C_norm_fr)

# t = 8.5145, df = 54.382, p-value = 2.338e-11
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  0.1791942 0.2917100
# sample estimates:
#  mean of x  mean of y 
# 0.32780823 0.09235614 


####################
par(las=2) # make label text perpendicular to axis
par(mar = c(10,4,4,2) + 2) # increase y-axis margin.
barplot(normMutCock$T_C_norm_fr, names.arg = normMutCock$Species, )

dev.off()
#####################


wilcox.test(normMutTer$T_C_norm_fr, normMutCock$T_C_norm_fr, conf.int = TRUE, exact = FALSE)

# W = 684, p-value = 2.078e-07
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   0.1821497 0.3051049
# sample estimates:
#   difference in location 
# 0.2477359 

# aov

aov_prep <- data.frame(T_C_norm = normMutTer[normMutTer$Family == 'Termitidae',]$T_C_norm_fr, Group = "Higher Ter")
aov_prep <- rbind(aov_prep, data.frame(T_C_norm = normMutTer[normMutTer$Family != 'Termitidae',]$T_C_norm_fr, Group = "Lower Ter"))
aov_prep <- rbind(aov_prep, data.frame(T_C_norm = normMutCock$T_C_norm_fr, Group = "Cockroach"))

aov.res <- aov(T_C_norm ~ Group, data = aov_prep) ## Google about it

summary(aov.res)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# Group        2 0.7759  0.3880   26.73 8.51e-09 ***
# Residuals   54 0.7836  0.0145                    
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov.res)
#                        diff        lwr        upr     p adj
# Higher Ter-Cockroach  0.26763715  0.17794260 0.35733170 0.0000000
# Lower Ter-Cockroach   0.18824733  0.08908627 0.28740839 0.0000829
# Lower Ter-Higher Ter -0.07938982 -0.17659998 0.01782033 0.1299455

################ BOXPLOT ############

Cocks_Norm_T_C <- aov_prep[aov_prep$Group == "Cockroach",]$T_C_norm
HighTers_Norm_T_C <- aov_prep[aov_prep$Group == "Higher Ter",]$T_C_norm
LowTers_Norm_T_C <- aov_prep[aov_prep$Group == "Lower Ter",]$T_C_norm

boxplot(Cocks_Norm_T_C, HighTers_Norm_T_C, LowTers_Norm_T_C,
        at = c(1,3,2),
        names = c("Cockroach", "Higher Ter", "Lower Ter"),
        col = c("chocolate4", "indianred3", "maroon"),
        ylab = "T_C",
        notch = F
        )









