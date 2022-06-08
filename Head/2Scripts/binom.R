library(dplyr)
normMut <- read.csv('~/Documents/lab/TermitesMutSpectrum/Body/2Derived/normMutTer.csv')

directional_pairs = list(
  c("A>C", "C>A"),
  c("A>G", "G>A"),
  c("A>U", "U>A"),
  c("C>G", "G>C"),
  c("C>U", "U>C"),
  c("G>U", "U>G")
)
reciprocal_pairs = list(
  c("A>C", "U>G"),
  c("A>G", "U>C"),
  c("A>U", "U>A"),
  c("C>G", "G>C"),
  c("C>U", "G>A"),
  c("G>U", "C>A")
)

binom <- function(normMut) {
  
}
binom.test(round(58.73343605547 / 11.2317224287485), round(58.73343605547 / 11.2317224287485) + round(170.872687704026 / 11.2317224287485)
)

NucSubst = c(
  "T>A",
  "T>C",
  "T>G",
  "A>T",
  "A>C",
  "A>G",
  "G>A",
  "G>C",
  "G>T",
  "C>A",
  "C>T",
  "C>G"
)

for (i in c(1:37)) { #swap 37 to the num of species

ObsToExp <- as.numeric(normMut[i, 38:49])

binomPrep <- data.frame(NucSubst)
binomPrep$ObsToExp = ObsToExp

sp_name <- normMut [i,1]
path <- sprintf('/home/glebo/Documents/lab/TermitesMutSpectrum/Body/2Derived/Cockroaches/binom_subs_per_sp/%s.csv', sp_name)
write.csv(binomPrep, path, row.names = FALSE)
}








