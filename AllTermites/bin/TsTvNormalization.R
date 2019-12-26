rm(list=ls(all=TRUE))

mut = read.table('../results/cockroaches11_19/MutSpecCockroachesTsTv.txt', sep='\t', header=TRUE)

plot(mut$sumOfSubs, mut$BranchLength)
plot(mut$sumOfSubs, mut$TsTv)
cor.test(mut$sumOfSubs, mut$TsTv, method = 'spearman')

pdf('../results/cockroaches11_19/TsTvNormalization.R.pdf')

plot(mut$sumOfSubs, mut$Ts)
plot(mut$sumOfSubs, mut$Tv)
plot(mut$Ts, mut$Tv)
plot(mut$sumOfSubs, mut$TsTv)

dev.off()

cor.test(mut$Ts, mut$Tv, method = 'spearman')

a = lm(Ts ~ Tv, mut)
summary(a)

res = a$residuals

cor.test(mut$BranchLength, res, method = 'spearman')

summary(glm(Cockroaches ~ TsTv + sumOfSubs, data=mut))
summary(glm(Cockroaches ~ TsTv*sumOfSubs, data=mut))
