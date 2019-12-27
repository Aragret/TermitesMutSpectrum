rm(list=ls(all=TRUE))

# install.packages("scales")
# install.packages('ggsci')
library(ggsci)
library(scales)
library(ggplot2)

mut = read.table('../results/cockroaches11_19/MutSpecCockroachesTsTv.txt', sep='\t', header=TRUE)

plot(mut$sumOfSubs, mut$BranchLength)
plot(mut$sumOfSubs, mut$TsTv)
cor.test(mut$sumOfSubs, mut$TsTv, method = 'spearman')


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

pdf('../results/cockroaches11_19/TsTvNormalization.R.pdf', width = 11)

plot(mut$sumOfSubs, mut$TsTv)
plot(mut$sumOfSubs, mut$Ts)
plot(mut$sumOfSubs, mut$Tv)
plot(mut$Ts, mut$Tv)

mut$Cockroaches = as.factor(mut$Cockroaches)

ggplot(mut, aes(Ts, Tv, col = Cockroaches)) + 
  geom_point()

ggplot(mut, aes(Ts, Tv, col = Cockroaches)) + 
  geom_point(aes(size = HigherTermites))

ggplot(mut[mut$HigherTermites == 1,], aes(Ts, Tv, col = Taxonomy)) +
  geom_point(aes(size = 1.5)) + # scale_colour_manual(values=pal_simpsons(palette = c("springfield"), alpha = 1))
  scale_color_simpsons(palette = c("springfield"))

unique(mut[mut$Cockroaches == 0, 'Taxonomy'])

manualcolors<-c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30',
                '#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695',
                '#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')

ggplot(mut[mut$Cockroaches == 0,], aes(Ts, Tv, col = Taxonomy)) +
  geom_point(aes(size = 1.5)) +
  scale_color_manual(values = manualcolors)

ggplot(mut[mut$Cockroaches == 0,], aes(Ts, Tv, col = Taxonomy)) +
  geom_point(aes(size = HigherTermites)) +
  scale_color_manual(values = manualcolors)

dev.off()

cor.test(mut$Ts, mut$Tv, method = 'spearman')

a = lm(Ts ~ Tv, mut)
summary(a)

res = a$residuals

cor.test(mut$BranchLength, res, method = 'spearman')

summary(glm(Cockroaches ~ TsTv, data=mut))
summary(glm(Cockroaches ~ sumOfSubs, data=mut))
summary(glm(Cockroaches ~ TsTv + sumOfSubs, data=mut))
# summary(glm(Cockroaches ~ TsTv*sumOfSubs, data=mut))

summary(glm(TsTv ~ Cockroaches + sumOfSubs, data=mut))

