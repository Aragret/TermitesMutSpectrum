

mut = read.table('results/4foldSubsRatesBrLen02.txt', header=TRUE, sep='\t')

lin_model = lm(C_T_rate ~ BranchLength, mut)
summary(lin_model)

lin_model$residuals

plot(mut$C_T_rate, mut$BranchLength)
abline(lm(mut$C_T_rate ~ mut$BranchLength))

res = as.data.frame(matrix(nrow = nrow(mut)))
for(i in 15:26){
  lin_model = lm(mut[, i] ~ mut$BranchLength)
  res = cbind(res, lin_model$residuals)
}

res = res[, -1]

names(res) = sub(' ', '', paste(names(mut[15:26]), '.res'))

mut = cbind(mut, res)

pdf('results/residuals.pdf')
plot(mut$C_T_rate, mut$BranchLength)
abline(lm(mut$C_T_rate ~ mut$BranchLength))

plot(mut$C_T_rate.res, mut$BranchLength)
abline(lm(mut$C_T_rate.res ~ mut$BranchLength))
dev.off()

# rate != 0
notZero = mut[mut$C_T_rate != 0,]
plot(notZero$C_T_rate, notZero$BranchLength)
abline(lm(notZero$C_T_rate ~ notZero$BranchLength))
lm(notZero$C_T_rate ~ notZero$BranchLength)

# plot(mut$C_T_rate.res, mut$BranchLength)
plot(mut$A_G_rate, mut$BranchLength)
plot(mut$A_G_rate.res, mut$BranchLength)

cor.test(mut$A_G_rate, mut$BranchLength)
cor.test(mut$A_G_rate.res, mut$BranchLength)

write.table(mut, 'results/4foldSubsRatesBrLen02Residuals.txt', row.names = FALSE, quote = FALSE, sep='\t')
