rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

library("Biostrings")
library("seqinr")
library(dplyr)

align_dir = '../../AllTermites/data/forAlign/'

setwd(align_dir)

filelist <- list.files()

VecOfSynFourFoldDegenerateSites <- c('CTT', 'CTC', 'CTA', 'CTG', 
                                     'GTT', 'GTC', 'GTA', 'GTG', 
                                     'TCT', 'TCC', 'TCA', 'TCG', 
                                     'CCT', 'CCC', 'CCA', 'CCG', 
                                     'ACT', 'ACC', 'ACA', 'ACG', 
                                     'GCT', 'GCC', 'GCA', 'GCG', 
                                     'CGT', 'CGC', 'CGA', 'CGG', 
                                     'GGT', 'GGC', 'GGA', 'GGG',
                                     'AGT', 'AGC', 'AGA', 'AGG')

majorGenes = c('ATP6', 'ATP8', 'COI', 'COII', 'COIII', 'CyB', 'ND3', 'ND6', 'ND2')
minorGenes = c('ND1', 'ND4', 'ND4L', 'ND5')

final = c()

for(f in filelist){
  one_line = c()
  # f = "ATP6.fasta"
  geneName = sub('.fasta', '', f)
  genes = read.fasta(f, as.string = TRUE)
  genes_df = data.frame(Species=names(genes), Seqs=unlist(getSequence(genes, as.string=T)))
  genes_df$Codons = lapply(as.character(genes_df$Seqs), s2c)
  genes_df$Codons = lapply(genes_df$Codons, splitseq)
  for(i in 1:nrow(genes_df)){
    # i = 1
    speciesName = as.character(genes_df$Species[i])
    codons = unlist(genes_df$Codons[i])
    numbers4f = which(toupper(codons) %in% VecOfSynFourFoldDegenerateSites)
    codons4f = toupper(codons[numbers4f])
    lastLetter = substring(codons4f, first = 3, last = 3)
    a = sum(nchar(gsub("[^A]", "", lastLetter)))
    t = sum(nchar(gsub("[^T]", "", lastLetter)))
    g = sum(nchar(gsub("[^G]", "", lastLetter)))
    c = sum(nchar(gsub("[^C]", "", lastLetter)))
    one_line = rbind(one_line, c(speciesName, geneName, a, t, g, c))
  }
  one_line = as.data.frame(one_line)
  names(one_line) = c('Species', 'Gene', 'A', 'T', 'G', 'C')
  if(f == "ATP6.fasta"){
    final = one_line
  }
  else{
    final = full_join(final, one_line)
  }
}

data = final %>% mutate(MajorStrand = case_when(.$Gene %in% majorGenes ~ 1, 
                                         .$Gene %in% minorGenes ~ 0))

write.table(data, '../../results/nucleotide_content06_20/ATGCforEachGene4fold.txt',
            sep = '\t', row.names = FALSE, quote = FALSE)

##################################################################################
# check Reticulitermes genes to be sure

atp6 = read.fasta('ATP6.fasta', as.string = TRUE)

seq = atp6$Reticulitermes_flavipes_IS58[1]

translate(s2c(seq), numcode = 5)

# [1] "M" "M" "S" "N" "L" "F" "S" "I" "F" "D" "P" "T" "T" "E" "I" "N" "S" "L" "P" "M" "N"
# [211] "S" "Y" "V" "F" "A" "I" "L" "S" "T" "L" "Y" "S" "S" "E" "V" "N" "*"

# /translation="MSNLFSIFDPTTEINSLPMNWTSTMVGLLLIPTSIWLTPSRNSM
#                      ALNLLMNKLHAEMKTILSKGNQNKGNSFIFTSLFLMILMNNFLGLFPYIFTSTSHLTL
#                      TLTLALPLWMTFMLFGWIKNTNHMFEHLVPQGTPTMLMPFMVIIETISNLIRPGTLAV
#                      RLTANMIAGHLLLTLLGNNGPNMSHTLLTVLIIAQILLLILESAVAIIQSYVFAILST
#                      LYSSEVN"

nd5 = read.fasta('ND5.fasta', as.string = TRUE)

seq2 = nd5$Reticulitermes_flavipes_IS58[1]

translate(s2c(seq2), numcode = 5)

# [1] "M" "M" "P" "I" "S" "I" "C" "F" "A" "S" "F" "I" "F" "L" "F" "G" "L" "G" "W" "V" "C"
# [568] "V" "L" "L" "F" "V" "L" "I" "Y"

# /translation="MMPISICFASFIFLFGLGWVCCFLGVYLVLSDLVYFVDWGIVSL
#                      NGSSVIMTFLFDWMSLLFLGFVFIISSLVILYSDDYMSGDFNIFRFIMLVLMFVVSMM
#                      FLIISPNMISILLGWDGLGLVSYCLVIYYQNVSSYGAGMLTVLSNRIGDVALLMVIAW
#                      MINFGSWNFIYYLEFMAGSVEMELISFLVVLAAMTSSAQIPFSSWLPAAMAAPTPVSA
#                      LVHSSTLVTAGVYLLIRFSPSFSCLLNTILLLVSALTMFMAGLGANFEYDLKKIIALS
#                      TLSQLGLMIMTVSVGLSSLAFFHLLTHALFKALLFMCAGGVIHSMGDSQDIRFMGGLS
#                      VYMPFTSSCLMVSSFALCGMPFLAGFYSSDFILEMISMSYVNVFGFLLLFVSTGLTVC
#                      YSFRLFYFVLCGDFNFVSLYSMVDTNYNMVMGMVGLLVLSVLGGGALMWLICPTPSVI
#                      CLPYYLSFLTFFFVSLGGFIGYEMAGFNFGDYLLSMYYYKVSSFSGSMWFMPFFSTYG
#                      VSFGPLGFGYSSMRVFDSGWMEYFGGQGLYWVLFNLGKFNQWVQHSSLKLFLGFFVMW
#                      VVLLFVLIY"

nd4 = read.fasta('ND4.fasta', as.string = TRUE)

seq3 = nd4$Reticulitermes_flavipes_IS58[1]

translate(s2c(seq3), numcode = 5)

# [1] "M" "L" "S" "F" "L" "C" "F" "L" "F" "F" "L" "T" "P" "L" "C" "I" "F" "P" "G" "S" "W"
# [442] "A" "V" "F" "W" "M" "*"

seq3

