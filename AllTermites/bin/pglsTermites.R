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

pdf('../results/nd6_22_01/pglsTermites.R.pdf')

##################################################################################
### Workers

filter_workers = mut[!is.na(mut$Worker) & mut$Worker != '?',]

summary(filter_workers$TsTv)

row.names(filter_workers) = filter_workers$Species

tree$tip.label = sub('__', '_', tree$tip.label)
tree$tip.label = sub('__', '_', tree$tip.label)
tree_w = treedata(tree, filter_workers, sort=T, warnings=T)$phy

data<-as.data.frame(treedata(tree_w, filter_workers, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

data$Ts = as.numeric(as.character(data$Ts))
data$Tv = as.numeric(as.character(data$Tv))

workers <- comparative.data(tree_w, data, Species, vcv=TRUE)

model1 = pgls(log(Tv) ~ Ts + Worker, workers, lambda="ML")
summary(model1)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5.2238 -0.9384  0.0519  1.0507  4.2926 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 0.986
# lower bound : 0.000, p = < 2.22e-16
# upper bound : 1.000, p = 0.44477
# 95.0% CI   : (0.931, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate  Std. Error t value  Pr(>|t|)    
# (Intercept)  3.17240295  0.40698289  7.7949 5.773e-14 ***
#   Ts           0.01544451  0.00063939 24.1552 < 2.2e-16 ***
#   Worker1     -0.04704836  0.35168506 -0.1338    0.8936    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.512 on 395 degrees of freedom
# Multiple R-squared: 0.5976,	Adjusted R-squared: 0.5956 
# F-statistic: 293.3 on 2 and 395 DF,  p-value: < 2.2e-16 

ggplot(data, aes(Ts, log(Tv), col = Worker)) +
  geom_point() + scale_color_manual(values = c('red', 'gray'))


#########################################################################
### Soldiers

filter_soldiers = mut[!is.na(mut$Soldier),]
filter_soldiers$Soldier = as.factor(filter_soldiers$Soldier)

row.names(filter_soldiers) = filter_soldiers$Species

tree$tip.label = sub('__', '_', tree$tip.label)
tree$tip.label = sub('__', '_', tree$tip.label)
tree_w = treedata(tree, filter_soldiers, sort=T, warnings=T)$phy

data<-as.data.frame(treedata(tree_w, filter_soldiers, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

data$Ts = as.numeric(as.character(data$Ts))
data$Tv = as.numeric(as.character(data$Tv))

soldiers <- comparative.data(tree_w, data, Species, vcv=TRUE)

model2 = pgls(log(Tv) ~ Ts + Soldier, soldiers, lambda="ML")
summary(model2)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.9322 -0.8864  0.0897  1.1262  5.0307 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 0.987
# lower bound : 0.000, p = < 2.22e-16
# upper bound : 1.000, p = 0.47118
# 95.0% CI   : (0.931, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)  3.2472791  0.4764783  6.8152  3.51e-11 ***
#   Ts           0.0155974  0.0006431 24.2535 < 2.2e-16 ***
#   Soldier1    -0.1266079  0.2824297 -0.4483    0.6542    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.533 on 397 degrees of freedom
# Multiple R-squared: 0.5974,	Adjusted R-squared: 0.5954 
# F-statistic: 294.5 on 2 and 397 DF,  p-value: < 2.2e-16 

ggplot(data, aes(Ts, log(Tv), col = Soldier)) +
  geom_point() + scale_color_manual(values = c('red', 'gray'))

#########################################################################
### Diet

filter_diet = mut[!is.na(mut$diet.Wood.Soil),]

# 0 - Soil, 1 - Wood
data = dummy_cols(filter_diet, 'diet.Wood.Soil', remove_first_dummy = TRUE)
row.names(data) = data$Species

tree$tip.label = sub('__', '_', tree$tip.label)
tree$tip.label = sub('__', '_', tree$tip.label)
tree_w = treedata(tree, data, sort=T, warnings=T)$phy

data<-as.data.frame(treedata(tree_w, data, sort=T, warnings=T)$data)

data$Ts = as.numeric(as.character(data$Ts))
data$Tv = as.numeric(as.character(data$Tv))

diet <- comparative.data(tree_w, data, Species, vcv=TRUE)

model3 = pgls(log(Tv) ~ Ts + diet.Wood.Soil_Wood, diet, lambda="ML")
summary(model3)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -5.2018 -0.8581  0.0970  1.1165  4.9943 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 0.986
# lower bound : 0.000, p = < 2.22e-16
# upper bound : 1.000, p = 0.43931
# 95.0% CI   : (0.927, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate  Std. Error t value  Pr(>|t|)    
# (Intercept)          3.23279294  0.40575133  7.9674 1.732e-14 ***
#   Ts                   0.01557973  0.00064322 24.2214 < 2.2e-16 ***
#   diet.Wood.Soil_Wood -0.11129159  0.13544719 -0.8217    0.4118    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.528 on 397 degrees of freedom
# Multiple R-squared: 0.598,	Adjusted R-squared: 0.5959 
# F-statistic: 295.2 on 2 and 397 DF,  p-value: < 2.2e-16

ggplot(data, aes(Ts, log(Tv), col = diet.Wood.Soil_Wood)) +
  geom_point() + 
  scale_colour_manual(values = c('red', 'gray'), labels = c('Soil', 'Wood'))

dev.off()
