b_aminoacids <- read.csv("~/Documents/lab/TermitesMutSpectrum/Body/3Results/Blattodea_aminoacids.txt")
sp <- c(b_aminoacids$Species)
# GAINERS ################
pro <- c(b_aminoacids$Pro)
his <- c(b_aminoacids$His)
gln <- c(b_aminoacids$Gln)
asn <- c(b_aminoacids$Asn)
lys <- c(b_aminoacids$Lys)
thr <- c(b_aminoacids$Thr)
# LOSERS ################
leu <- c(b_aminoacids$Leu)
phe <- c(b_aminoacids$Phe)
cys <- c(b_aminoacids$Cys)
trp <- c(b_aminoacids$Trp)
gly <- c(b_aminoacids$Gly)
val <- c(b_aminoacids$Val)
##########################
gainers <- c()
species <- c()
losers <- c()
for (i in 1:length(sp)){
  gainers <- c(gainers, pro[i] + his[i] + gln[i] + asn[i] + lys[i] + thr[i])
  losers <- c(losers, leu[i] + phe[i] + cys[i] + trp[i] + gly[i] + val[i])
  species <- c(species, sp[i])
}

table <- data.frame(Species = species, Gainers = gainers, Losers = losers)
write.csv2(table, "~/Documents/lab/TermitesMutSpectrum/Body/3Results/Gains_and_Loses.txt", row.names = FALSE)
# CLEAN ###################
rm(list = ls())