rm(list = ls(all=TRUE))

if (!require(caper)) install.packages("caper")
if (!require(geiger)) install.packages("geiger")
if (!require(fastDummies)) install.packages("fastDummies")
if (!require(ggplot2)) install.packages("ggplot2")

library(caper)
library(geiger)
library(fastDummies)
library(ggplot2)

mut = read.table('../results/nd6_22_01/mutSpectrumFractions.txt', header = TRUE, sep='\t')

mut$TsTv = (mut$C_T + mut$A_G + mut$T_C + mut$G_A) / (mut$G_T + mut$A_T + mut$G_C +
                                                        mut$A_C + mut$C_A + mut$C_G +
                                                        mut$T_A + mut$T_G)

mut$Ts = mut$C_T + mut$A_G + mut$T_C + mut$G_A
mut$Tv = mut$G_T + mut$A_T + mut$G_C +  mut$A_C + mut$C_A + mut$C_G + mut$T_A + mut$T_G

mut = mut[mut$Tv != 0,] # for log()

tree = read.tree('../results/nd6_22_01/phylogeny_rooted/iqtree_mito.treefile')

###################################################################################

unique(mut$Taxonomy)

for(i in 1:nrow(mut)){
  if(as.character(mut[i,]$Taxonomy) %in% c('Kalotermitidae', 'Stolotermitidae', 'Serritermitidae',
                             'Archotermopsidae', 'Termitogetoninae', 'Prorhinotermitinae')){
    mut[i, 'MoreSocial'] = 0
  }
  if(as.character(mut[i,]$Taxonomy) %in% c('Coptotermitinae', 'Nasutitermitinae', 'Amitermes-group',
                             'Macrotermitinae', 'Apicotermitinae', 'Rhinotermitinae',
                             'Rhinotermitidae', 'Microcerotermes', 'Termes-group', 'Promiro',
                             'Cubitermitinae', 'Foraminitermitinae', 'Syntermitinae',
                             'Pericapritermes-group', 'Sphaerotermitinae', 'Neocapri-group',
                             'Cephalo-group', 'Mastotermitidae', 'Prorhinotermitinae',
                             'pericapritermes-group')){
    mut[i, 'MoreSocial'] = 1
  }
}

summary(as.factor(mut$MoreSocial))

data = mut[!is.na(mut$MoreSocial),]
row.names(data) = data$Species

pdf('../results/nd6_22_01/pglsTermiteFamilies.R.pdf')

ggplot(data, aes(Ts, log(Tv), col=as.factor(MoreSocial))) +
  geom_point() + scale_color_manual(values = c('red', 'gray'))

dev.off()
#############################################################################

tree$tip.label = sub('__', '_', tree$tip.label)
tree$tip.label = sub('__', '_', tree$tip.label)
tree_w = treedata(tree, data, sort=T, warnings=T)$phy

data<-as.data.frame(treedata(tree_w, data, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

data$Ts = as.numeric(as.character(data$Ts))
data$Tv = as.numeric(as.character(data$Tv))

data_comp <- comparative.data(tree_w, data, Species, vcv=TRUE)

model = pgls(log(Tv) ~ Ts + MoreSocial, data_comp, lambda="ML")
summary(model)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -4.4658 -0.8262  0.1738  1.0910  4.3674 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 0.985
# lower bound : 0.000, p = < 2.22e-16
# upper bound : 1.000, p = 0.40485
# 95.0% CI   : (0.928, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate  Std. Error t value  Pr(>|t|)    
# (Intercept)  3.31323573  0.39980020  8.2872 1.776e-15 ***
#   Ts           0.01541458  0.00063737 24.1846 < 2.2e-16 ***
#   MoreSocial1 -0.49550603  0.40676676 -1.2182    0.2239    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.505 on 395 degrees of freedom
# Multiple R-squared: 0.5992,	Adjusted R-squared: 0.5972 
# F-statistic: 295.3 on 2 and 395 DF,  p-value: < 2.2e-16 
