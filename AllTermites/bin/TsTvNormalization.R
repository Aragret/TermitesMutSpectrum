rm(list=ls(all=TRUE))

library(ggplot2)

mut = read.table('../results/cockroaches11_19/MutSpecCockroachesTsTv.txt', sep='\t', header=TRUE)

plot(mut$sumOfSubs, mut$BranchLength)
plot(mut$sumOfSubs, mut$TsTv)
cor.test(mut$sumOfSubs, mut$TsTv, method = 'spearman')


higher_termites = c('Apicotermitinae', 'Cubitermitinae', 'Foraminitermitinae',
                    'Macrotermitinae', 'Nasutitermitinae', 'Sphaerotermitinae',
                    'Syntermitinae', 'Termes-group')

for(i in 1:nrow(mut)){
  # i = 1
  if(mut[i, 'Taxonomy'] %in% higher_termites){
    mut[i, 'Higher_termites'] = 1
  }
  else{mut[i, 'Higher_termites'] = 0}
}

mut$Higher_termites = as.factor(mut$Higher_termites)

pdf('../results/cockroaches11_19/TsTvNormalization.R.pdf')

plot(mut$sumOfSubs, mut$TsTv)
plot(mut$sumOfSubs, mut$Ts)
plot(mut$sumOfSubs, mut$Tv)
plot(mut$Ts, mut$Tv)

mut$Cockroaches = as.factor(mut$Cockroaches)

ggplot(mut, aes(Ts, Tv, col = Cockroaches)) + 
  geom_point()

ggplot(mut, aes(Ts, Tv, col = Cockroaches)) + 
  geom_point(aes(size = Higher_termites))


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

