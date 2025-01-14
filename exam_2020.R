



apa <- ape::read.FASTA("~/Downloads/DNAseq_20200824.fasta")


dna <- as.character(apa)
dna$DNAsequence[dna$DNAsequence == "t"] <- "u"
dna

dna_collapse <- paste(dna$DNAsequence, collapse = "")
dna_collapse
for(i in 1:3){
  print(ape::trans(ape::as.DNAbin(dna_collapse),codonstart = i))
}

ape::as.DNAbin(dna)


poop <- ape::trans(ape::as.DNAbin(dna_collapse),code = 1,codonstart = 1)
poop$DNAsequence
dna_collapse <- paste(dna$DNAsequence, collapse = "")

BiocManager::install("Biostrings")
  
library(Biostrings)
res <-Biostrings::RNAString(dna_collapse)


translate(res)


