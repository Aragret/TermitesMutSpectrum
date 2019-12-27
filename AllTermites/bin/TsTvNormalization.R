rm(list=ls(all=TRUE))

# install.packages("scales")
install.packages('ggsci')
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
  geom_point() + # scale_colour_manual(values=pal_simpsons(palette = c("springfield"), alpha = 1))
  scale_color_simpsons(palette = c("springfield"))

unique(mut[mut$Cockroaches == 0, 'Taxonomy'])

manualcolors<-c('black','forestgreen', 'red2', 'orange', 'cornflowerblue', 
                'magenta', 'darkolivegreen4',  
                'indianred1', 'tan4', 'darkblue', 
                'mediumorchid1','firebrick4',  'yellowgreen', 'lightsalmon', 'tan3',
                "tan1",'darkgray', 'wheat4', '#DDAD4B', 'chartreuse', 'seagreen1',
                'moccasin', 'mediumvioletred', 'seagreen','cadetblue1', "darkolivegreen1")

ggplot(mut[mut$Cockroaches == 0,], aes(Ts, Tv, col = Taxonomy)) +
  geom_point() +
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

