library(ggplot2)
library(dplyr)
library(ggpubr)

mutCockroaches = read.table('~/Documents/lab/TermitesMutSpectrum/Body/3Results/Cockroaches.MutSpecData.txt', header = TRUE)
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
      A_T_norm_fraction = A_T_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      A_C_norm_fraction = A_C_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      A_G_norm_fraction = A_G_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      T_A_norm_fraction = T_A_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      T_C_norm_fraction = T_C_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      T_G_norm_fraction = T_G_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      G_A_norm_fraction = G_A_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      G_C_norm_fraction = G_C_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      G_T_norm_fraction = G_T_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      C_A_norm_fraction = C_A_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      C_T_norm_fraction = C_T_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm),
      C_G_norm_fraction = C_G_norm / (T_A_norm + T_C_norm + T_G_norm + A_T_norm + A_C_norm + A_G_norm + G_A_norm + G_C_norm + G_T_norm + C_A_norm + C_T_norm + C_G_norm)
    )
}

normMutCockroaches <- test_prep(mutCockroaches)
normMutTer <- test_prep(mutTer)
normMutCockroaches[is.na(normMutCockroaches)] = 0
normMutCockroaches[6, "T_C_norm_fraction"] = 0.123345141 #calculated manually; looks like that's a seq, not a whole genome!!!!!!! delete?
t.test(normMutTer$T_C_norm_fraction, normMutCockroaches$T_C_norm_fraction)

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
barplot(normMutCockroaches$T_C_norm_fraction, names.arg = normMutCockroaches$Species, )

dev.off()
#####################


wilcox.test(normMutTer$T_C_norm_fraction, normMutCockroaches$T_C_norm_fraction, conf.int = TRUE, exact = FALSE)

# W = 684, p-value = 2.078e-07
# alternative hypothesis: true location shift is not equal to 0
# 95 percent confidence interval:
#   0.1821497 0.3051049
# sample estimates:
#   difference in location 
# 0.2477359 

# aov

aov_prep <- data.frame(T_C_norm = normMutTer[normMutTer$Family == 'Termitidae',]$T_C_norm_fraction, Group = "Higher Ter")
aov_prep <- rbind(aov_prep, data.frame(T_C_norm = normMutTer[normMutTer$Family != 'Termitidae',]$T_C_norm_fraction, Group = "Lower Ter"))
aov_prep <- rbind(aov_prep, data.frame(T_C_norm = normMutCockroaches$T_C_norm_fraction, Group = "Cockroaches"))

aov.res <- aov(T_C_norm ~ Group, data = aov_prep) ## Google about it

summary(aov.res)
#             Df Sum Sq Mean Sq F value   Pr(>F)    
# Group        2 0.7759  0.3880   26.73 8.51e-09 ***
# Residuals   54 0.7836  0.0145                    
# ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
TukeyHSD(aov.res)
#                        diff        lwr        upr     p adj
# Higher Ter-Cockroaches  0.26763715  0.17794260 0.35733170 0.0000000
# Lower Ter-Cockroaches   0.18824733  0.08908627 0.28740839 0.0000829
# Lower Ter-Higher Ter -0.07938982 -0.17659998 0.01782033 0.1299455

################ BOXPLOT ############

Cockroaches_Norm_T_C <- aov_prep[aov_prep$Group == "Cockroaches",]$T_C_norm
HighTers_Norm_T_C <- aov_prep[aov_prep$Group == "Higher Ter",]$T_C_norm
LowTers_Norm_T_C <- aov_prep[aov_prep$Group == "Lower Ter",]$T_C_norm

boxplot(Cockroaches_Norm_T_C, HighTers_Norm_T_C, LowTers_Norm_T_C,
        at = c(1,3,2),
        names = c("Cockroaches", "Higher Ter", "Lower Ter"),
        col = c("chocolate4", "indianred3", "maroon"),
        ylab = "T_C",
        notch = F
        )
write.csv(normMutCockroaches, '~/Documents/lab/TermitesMutSpectrum/Body/2Derived/normMutCockroaches.csv', row.names = FALSE)
write.csv(normMutTer, '~/Documents/lab/TermitesMutSpectrum/Body/2Derived/normMutTer.csv', row.names = FALSE)



#####12SUBBARPLOT######################



CockroachesPlot12Subs <- data.frame(Subs = factor(c('Tn_An', 'Tn_Cn', 'Tn_Gn', 'An_Tn', 'An_Cn', 'An_Gn', 'Cn_An', 'Cn_Tn', 'Cn_Gn', 'Gn_An', 'Gn_Cn', 'Gn_Tn'),
                                           levels=c('Tn_An', 'Tn_Cn', 'Tn_Gn', 'Gn_Tn', 'Gn_Cn', 'Gn_An', 'Cn_Tn', 'Cn_Gn', 'Cn_An','An_Tn', 'An_Gn', 'An_Cn')), 
                             Values = c(sum(normMutCockroaches$A_T_norm_fraction), sum(normMutCockroaches$A_G_norm_fraction), sum(normMutCockroaches$A_C_norm_fraction), 
                                         sum(normMutCockroaches$T_A_norm_fraction), sum(normMutCockroaches$T_G_norm_fraction), sum(normMutCockroaches$T_C_norm_fraction), 
                                         sum(normMutCockroaches$G_T_norm_fraction), sum(normMutCockroaches$G_A_norm_fraction), sum(normMutCockroaches$G_C_norm_fraction), 
                                         sum(normMutCockroaches$C_T_norm_fraction), sum(normMutCockroaches$C_G_norm_fraction), sum(normMutCockroaches$C_A_norm_fraction)))

TerPlot12Subs <- data.frame(Subs = factor(c('Tn_An', 'Tn_Cn', 'Tn_Gn', 'An_Tn', 'An_Cn', 'An_Gn', 'Cn_An', 'Cn_Tn', 'Cn_Gn', 'Gn_An', 'Gn_Cn', 'Gn_Tn'),
                                          levels=c('Tn_An', 'Tn_Cn', 'Tn_Gn', 'Gn_Tn', 'Gn_Cn', 'Gn_An', 'Cn_Tn', 'Cn_Gn', 'Cn_An','An_Tn', 'An_Gn', 'An_Cn')), 
                            Values = c(sum(normMutTer$A_T_norm_fraction), sum(normMutTer$A_G_norm_fraction), sum(normMutTer$A_C_norm_fraction), 
                                       sum(normMutTer$T_A_norm_fraction), sum(normMutTer$T_G_norm_fraction), sum(normMutTer$T_C_norm_fraction), 
                                       sum(normMutTer$G_T_norm_fraction), sum(normMutTer$G_A_norm_fraction), sum(normMutTer$G_C_norm_fraction), 
                                       sum(normMutTer$C_T_norm_fraction), sum(normMutTer$C_G_norm_fraction), sum(normMutTer$C_A_norm_fraction)))
#'Tn_An', 'Tn_Cn', 'Tn_Gn', 'An_Tn', 'An_Cn', 'An_Gn', 'Cn_An', 'Cn_Tn', 'Cn_Gn', 'Gn_An', 'Gn_Cn', 'Gn_Tn'
#'A_T', 'A_G', 'A_C', 'T_A', 'T_G', 'T_C', 'G_T', 'G_A', 'G_C', 'C_T', 'C_G', 'C_A'
Cockroaches12 <- ggplot(CockroachesPlot12Subs, aes(x=Subs, y=Values)) + 
  geom_bar(stat = "identity", aes(fill=Subs)) + theme(legend.position="none", axis.text.x = element_text(size=15)) + ggtitle('Cockroaches') + 
  xlab('Substitution types') + ylab('Substitution rate') + ylim(0,11.5)

Ter12 <- ggplot(TerPlot12Subs, aes(x=Subs, y=Values)) + ggtitle('Termites') + 
  geom_bar(stat = "identity", aes(fill=Subs)) + theme(legend.position="none", axis.text.x = element_text(size=15)) + 
  xlab('Substitution types') + ylab('Substitution rate') 

ggarrange(Cockroaches12, Ter12,
         labels=c('A','B'),
         ncol=1,
         nrow=2)

#####STACKEDPLOT##############

sps <- c()
for (species in normMutCockroaches$Species){
  sps <- append(sps, rep(species,4))
}
#'A', 'T', 'G', 'C'
#'T', 'A', 'C', 'G'
Mutations <- rep(c('T', 'A', 'C', 'G'), 4)
val <- c()
counter <- 1
for (species in normMutCockroaches$Species){
  val <- append(val, normMutCockroaches$A[counter])
  val <- append(val, normMutCockroaches$T[counter])
  val <- append(val, normMutCockroaches$G[counter])
  val <- append(val, normMutCockroaches$C[counter])
  counter = counter + 1
}

stackedplotCockroaches <- data.frame(sps, Mutations, val)

StackedCockroaches <- ggplot(stackedplotCockroaches, aes(fill=Mutations, y=val, x=sps)) + ggtitle('Cockroaches') + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=15)) + 
  xlab('Species') + ylab('Substitution rate')

sps <- c()
for (fams in normMutTer$Family){
  sps <- append(sps, rep(fams,4))
}
Mutations <- c()
Mutations <- rep(c('T', 'A', 'C', 'G'), 37)
val <- c()
counter <- 1
for (species in normMutTer$Species){
  val <- append(val, normMutTer$A[counter])
  val <- append(val, normMutTer$T[counter])
  val <- append(val, normMutTer$G[counter])
  val <- append(val, normMutTer$C[counter])
  counter = counter + 1
}

stackedplotTer <- data.frame(sps, Mutations, val)

StackedTer <- ggplot(stackedplotTer, aes(fill=Mutations, y=val, x=sps)) + ggtitle('Termites') + 
  geom_bar(position="fill", stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=15)) +
  xlab('Families') + ylab('Substitution rate')

ggarrange(StackedCockroaches, StackedTer,
          labels=c('A', 'B'),
          ncol=2,
          nrow=1,
          align='h')




