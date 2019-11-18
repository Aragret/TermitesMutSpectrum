

mut = read.table('results/4foldSubsRatesBrLen02.txt', header=TRUE, sep='\t')

plot(mut$C_T_rate, mut$BranchLength)

lin_model = lm(C_T_rate ~ BranchLength, mut)
summary(lin_model)

lin_model$residuals

res = as.data.frame(matrix(nrow = nrow(mut)))
for(i in 15:26){
  lin_model = lm(mut[, i] ~ mut$BranchLength)
  res = cbind(res, lin_model$residuals)
}

res = res[, -1]

names(res) = sub(' ', '', paste(names(mut[15:26]), '.res'))

mut = cbind(mut, res)

plot(mut$C_T_rate.res, mut$BranchLength)
plot(mut$A_G_rate, mut$BranchLength)
plot(mut$A_G_rate.res, mut$BranchLength)

cor.test(mut$A_G_rate, mut$BranchLength)
cor.test(mut$A_G_rate.res, mut$BranchLength)

write.table(mut, 'results/4foldSubsRatesBrLen02Residuals.txt')
