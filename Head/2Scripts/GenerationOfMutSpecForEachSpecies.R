rm(list=ls(all=TRUE))

### neutral ATGC
NeutralATGC = read.table('../../Body/3Results/Cockroaches.Normalization.NeutralATGC.txt', header = TRUE, sep='\t')

####### READ
MUT = read.table("../../Body/3Results/Mutational_spectra_in_Cockroaches.txt", header = TRUE)
length(unique(MUT$Species)) #37

##### FILTER 1: to take only normal substitutions and filter out species with too high fraction (> 5%) of unnormal substitutions
VecOfNormalSubstitutions <- c('A_C','C_A',
                              'A_G','G_A',
                              'C_G','G_C',
                              'C_T','T_C',
                              'G_T','T_G',
                              'T_A','A_T')
nrow(MUT)
table(MUT$Subs)   # MANY CRAPPY SUBSTITUTIONS!!!!!!!!!!!!!!!!!! WHY?????????????????????
SP = data.frame(table(MUT$Species)); names(SP) = c('Species','NumberOfAllSubst')
SPN = data.frame(table(MUT[MUT$Subs %in% VecOfNormalSubstitutions,]$Species)); names(SPN) = c('Species','NumberOfNormalSubst')
SP = merge(SP,SPN); SP$FractionOfNormal = SP$NumberOfNormalSubst/SP$NumberOfAllSubst
hist(SP$FractionOfNormal)
summary(SP$FractionOfNormal) # how many to delete? ask to have more than 95% of substitutions as normal
# SpeciesToDelete = SP[SP$FractionOfNormal <=0.95,]$Species; length(SpeciesToDelete)
# MUT = MUT[!MUT$Species %in% SpeciesToDelete,]
# MUT = MUT[MUT$Subs %in% VecOfNormalSubstitutions,]

##### FILTER 2: Synonymous Substitutions
nrow(MUT) # 1884
MUT = MUT[as.character(MUT$AncestralAA) == as.character(MUT$DescendantAA),]; nrow(MUT) # 1675
table(MUT$AncestralAA)
table(MUT$DescendantAA)

MUT = merge(MUT,NeutralATGC[,c(1:6)], by = c("Species", "Gene"))
nrow(MUT) # 1675

taxons = read.csv('../../Body/2Derived/CockroachesFamilies.csv')
MUT = merge(MUT, taxons[, 1:2])

write.table(MUT, '../../Body/3Results/Cockroaches.MutSpecData.txt', quote = FALSE, row.names = FALSE)