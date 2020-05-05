rm(list = ls(all=TRUE))

# library(devtools)
# install_github("JeroenSmaers/evomap")
library(evomap)
library(phytools)
library(geiger)
library(nlme)


mut = read.table('../results/nd6_22_01/mutSpectrumFractions.txt', header = TRUE, sep='\t')

mut$TsTv = (mut$C_T + mut$A_G + mut$T_C + mut$G_A) / (mut$G_T + mut$A_T + mut$G_C +
                                                        mut$A_C + mut$C_A + mut$C_G +
                                                        mut$T_A + mut$T_G)

mut$Ts = mut$C_T + mut$A_G + mut$T_C + mut$G_A
mut$Tv = mut$G_T + mut$A_T + mut$G_C +  mut$A_C + mut$C_A + mut$C_G + mut$T_A + mut$T_G

tree = read.tree('../results/nd6_22_01/phylogeny_rooted/iqtree_mito.treefile')

##################################################################################
### Workers

filter_workers = mut[!is.na(mut$Worker),]

summary(filter_workers$TsTv)

row.names(filter_workers) = filter_workers$Species

tree$tip.label = sub('__', '_', tree$tip.label)
tree$tip.label = sub('__', '_', tree$tip.label)
tree_w = treedata(tree, filter_workers, sort=T, warnings=T)$phy


data<-as.data.frame(treedata(tree_w, filter_workers, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

nrow(data[data$Worker == 0,]) # 72

Workers <- which(data$Worker == 1)
NonWorkers <- which(data$Worker == 0)

data = data[, c('Ts', 'Tv')]
data$Ts = as.integer(data$Ts)
data$Tv = as.integer(data$Tv)

plot(data$Tv ~ data$Ts, col="white", pch=19, xlab="", ylab="", asp=1, cex.lab=2) 
pGLS.plotGrade("Tv", "Ts", data, tree_w, model="lambda", group=Workers, col="green", lwd=5, cex=1, pch=19) 
pGLS.plotGrade("Tv", "Ts", data, tree_w, model="lambda", group=NonWorkers, col="grey", lwd=5, cex=1, pch=19)


#For differences in slope: 
grpS<-rep("A",length(rownames(data))) 
grpS[Workers]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(data)

#For differences in intercept: 
grpI<-rep("A",length(rownames(data))) 
grpI[Workers]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

Model<-model.matrix(as.formula(Tv ~ Ts),data)

#(1) Differences in slopes, holding intercept constant: 
Model_S<-model.matrix(as.formula(Tv ~ grpS:Ts), data) 

#(2) Differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Tv ~ grpI + Ts), data) 

#(3) Differences in slopes and differences in intercept: 
Model_SI<-model.matrix(as.formula(Tv ~ grpI + grpS:Ts),data)


#(1) Differences in slopes, holding intercept constant:
gls.ancova(Tv ~ Ts, vcv(tree_w), Model, Model_S)

# df       Sum Sq Mean Sum Sq F value Pr(>F)
# FullModel     3 4856574.2508  12171.8653  2.4634 0.1173
# ReducedModel  2 4886558.3988   12216.396  

#(2) Differences in intercept, holding slopes constant:
gls.ancova(Tv ~ Ts, vcv(tree_w), Model, Model_I)

# df       Sum Sq Mean Sum Sq F value Pr(>F)
# FullModel     3 4883697.2535  12239.8427  0.2338  0.629
# ReducedModel  2 4886558.3988   12216.396 

#(3) Differences in slopes and differences in intercept:
gls.ancova(Tv ~ Ts, vcv(tree_w), Model, Model_SI)

# df       Sum Sq Mean Sum Sq F value Pr(>F)
# FullModel     4 4856571.0075  12202.4397  1.2287 0.2938
# ReducedModel  2 4886558.3988   12216.396 

# check phylogenetic signal
phylosig(tree_w, data$Ts, method = 'lambda') # lambda : 0.315021 
 
phylosig(tree_w, data$Tv, method = 'lambda') # lambda : 0.124083 

###################################################################################
### Soldiers 

filter_soldiers = mut[!is.na(mut$Soldier),]

summary(filter_soldiers$TsTv)

row.names(filter_soldiers) = filter_soldiers$Species

tree$tip.label = sub('__', '_', tree$tip.label)
tree$tip.label = sub('__', '_', tree$tip.label)
tree_w = treedata(tree, filter_soldiers, sort=T, warnings=T)$phy


data<-as.data.frame(treedata(tree_w, filter_soldiers, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

nrow(data[data$Soldier == 0,]) # 45

Soldier <- which(data$Soldier == 1)
NonSoldier <- which(data$Soldier == 0)

data = data[, c('Ts', 'Tv')]
data$Ts = as.integer(data$Ts)
data$Tv = as.integer(data$Tv)

plot(data$Tv ~ data$Ts, col="white", pch=19, xlab="", ylab="", asp=1, cex.lab=2) 
pGLS.plotGrade("Tv", "Ts", data, tree_w, model="lambda", group=Soldier, col="green", lwd=5, cex=1, pch=19) 
pGLS.plotGrade("Tv", "Ts", data, tree_w, model="lambda", group=NonSoldier, col="grey", lwd=5, cex=1, pch=19)


#For differences in slope: 
grpS<-rep("A",length(rownames(data))) 
grpS[Soldier]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(data)

#For differences in intercept: 
grpI<-rep("A",length(rownames(data))) 
grpI[Soldier]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

Model<-model.matrix(as.formula(Tv ~ Ts),data)

#(1) Differences in slopes, holding intercept constant: 
Model_S<-model.matrix(as.formula(Tv ~ grpS:Ts), data) 

#(2) Differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Tv ~ grpI + Ts), data) 

#(3) Differences in slopes and differences in intercept: 
Model_SI<-model.matrix(as.formula(Tv ~ grpI + grpS:Ts),data)


#(1) Differences in slopes, holding intercept constant:
gls.ancova(Tv ~ Ts, vcv(tree_w), Model, Model_S)

# df       Sum Sq Mean Sum Sq F value Pr(>F)
# FullModel     3 4880628.6548   12232.152  0.4848 0.4867
# ReducedModel  2 4886558.3988   12216.396   

#(2) Differences in intercept, holding slopes constant:
gls.ancova(Tv ~ Ts, vcv(tree_w), Model, Model_I)

# df       Sum Sq Mean Sum Sq F value Pr(>F)
# FullModel     3 4878183.9498  12226.0249   0.685 0.4084
# ReducedModel  2 4886558.3988   12216.396 

#(3) Differences in slopes and differences in intercept:
gls.ancova(Tv ~ Ts, vcv(tree_w), Model, Model_SI)

# df       Sum Sq Mean Sum Sq F value Pr(>F)
# FullModel     4 4867431.1834  12229.7266   0.782 0.4582
# ReducedModel  2 4886558.3988   12216.396   

###################################################################################
### Diet 

filter_diet = mut[!is.na(mut$diet.Wood.Soil),]

summary(filter_diet$TsTv)

row.names(filter_diet) = filter_diet$Species

tree$tip.label = sub('__', '_', tree$tip.label)
tree$tip.label = sub('__', '_', tree$tip.label)
tree_w = treedata(tree, filter_diet, sort=T, warnings=T)$phy


data<-as.data.frame(treedata(tree_w, filter_diet, sort=T, warnings=T)$data)
data$Species = as.character(data$Species)

nrow(data[data$diet.Wood.Soil == 'Wood',]) # 267

WoodFeeders <- which(data$diet.Wood.Soil == 'Wood')
SoilFeeders <- which(data$diet.Wood.Soil == 'Soil')

data = data[, c('Ts', 'Tv')]
data$Ts = as.integer(data$Ts)
data$Tv = as.integer(data$Tv)

plot(data$Tv ~ data$Ts, col="white", pch=19, xlab="", ylab="", asp=1, cex.lab=2) 
pGLS.plotGrade("Tv", "Ts", data, tree_w, model="lambda", group=WoodFeeders, col="green", lwd=5, cex=1, pch=19) 
pGLS.plotGrade("Tv", "Ts", data, tree_w, model="lambda", group=SoilFeeders, col="grey", lwd=5, cex=1, pch=19)


#For differences in slope: 
grpS<-rep("A",length(rownames(data))) 
grpS[SoilFeeders]<-"B" 
grpS<-as.factor(grpS) 
names(grpS)<-rownames(data)

#For differences in intercept: 
grpI<-rep("A",length(rownames(data))) 
grpI[SoilFeeders]<-"B" 
grpI<-as.factor(grpI) 
names(grpI)<-rownames(data)

Model<-model.matrix(as.formula(Tv ~ Ts),data)

#(1) Differences in slopes, holding intercept constant: 
Model_S<-model.matrix(as.formula(Tv ~ grpS:Ts), data) 

#(2) Differences in intercept, holding slopes constant: 
Model_I<-model.matrix(as.formula(Tv ~ grpI + Ts), data) 

#(3) Differences in slopes and differences in intercept: 
Model_SI<-model.matrix(as.formula(Tv ~ grpI + grpS:Ts),data)


#(1) Differences in slopes, holding intercept constant:
gls.ancova(Tv ~ Ts, vcv(tree_w), Model, Model_S)

# df       Sum Sq Mean Sum Sq F value Pr(>F)
# FullModel     3 4848667.9054  12152.0499   3.118 0.0782
# ReducedModel  2 4886558.3988   12216.396               

#(2) Differences in intercept, holding slopes constant:
gls.ancova(Tv ~ Ts, vcv(tree_w), Model, Model_I)

# df       Sum Sq Mean Sum Sq F value Pr(>F)
# FullModel     3 4883155.9875  12238.4862   0.278 0.5983
# ReducedModel  2 4886558.3988   12216.396    

#(3) Differences in slopes and differences in intercept:
gls.ancova(Tv ~ Ts, vcv(tree_w), Model, Model_SI)

# df       Sum Sq Mean Sum Sq F value Pr(>F)
# FullModel     4 4825503.5546  12124.3808  2.5179 0.0819
# ReducedModel  2 4886558.3988   12216.396 
