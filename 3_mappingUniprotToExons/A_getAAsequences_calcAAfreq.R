####This file maps AA sequences onto exons, as well as calculating AA frequencies
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(Biostrings)
library(readr)

setwd("C:/tissueAbSplicing")

####Read file with UniProt IDs, transcript types, other protein coding related parameters
#Download from UCSC table browser
#All GENCODE v47
#-> hg38.wgEncodeGencodeCompV47 (allow selection from other tables) (also check all columns)
#-> .wgEncodeGencodeAttrsV47 fields (get 'transcriptType'), hg38.wgEncodeGencodeUniProtV47 (get 'acc')
GENCODE_uniprot <- read.table("C:/dataprocessing/database/uniprot/GENCODEv47_comprehensive_uniprot", sep = "\t", header = TRUE, comment = "$")
colnames(GENCODE_uniprot) <- gsub("^.*\\.", "", colnames(GENCODE_uniprot))

#split UCSC table browser tx-level information into exon-level information
#also obtain CDS sequences
GENCODE_split <- GENCODE_uniprot %>%
  filter(chrom %in% paste0("chr", c(1:22,"X","Y"))) %>%
  filter(transcriptType %in% paste0("protein_coding")) %>%
  select(geneSym = name2, txID = name, chr = chrom, strand, txStart, txEnd, cdsStart, cdsEnd, cdsStartStat, cdsEndStat, exonStart = exonStarts, exonEnd = exonEnds, exonFrame = exonFrames, transcriptType, GENCODE.level = level) %>%
  mutate(txID = gsub("\\..*?$", "", txID)) %>%
  mutate_at(.vars = c("exonStart", "exonEnd", "exonFrame"), .funs = function(x){gsub(",$", "", x)}) %>%
  separate_rows(., exonStart, exonEnd, exonFrame, sep = ",") %>%
  filter(exonFrame != -1) %>%
  group_by(txID) %>%
  mutate(totalExNo = n(), exNo = ifelse(strand == "+", seq.int(1, n(), 1), seq.int(n(), 1, -1))) %>%
  arrange(geneSym, txID, exNo) %>%
  ungroup() %>%
  mutate(
    exonCDSstart = pmax(as.numeric(cdsStart), as.numeric(exonStart))+1,
    exonCDSend = pmin(as.numeric(cdsEnd), as.numeric(exonEnd))
  ) %>%
  mutate(
    coords = paste0(chr, ":", exonCDSstart, "-", exonCDSend, ":", strand)
  ) %>%
  mutate(
    seq = as.character(getSeq(
      BSgenome.Hsapiens.UCSC.hg38,
      GRanges(coords)
    )) #exon CDS sequences
    ) %>%
  mutate(seq = case_when(
    exonFrame == 0 ~ seq,
    exonFrame == 1 ~ substr(seq, 3, nchar(seq)),
    exonFrame == 2 ~ substr(seq, 2, nchar(seq)))
  ) %>% #remove leading nucleotides that are part of a codon separated by exon junctions
  mutate(
    seq = substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
  ) %>% #remove trailing nucleotides that are part of a separated codon
  mutate(
    AAseq = as.character(translate(DNAStringSet(seq)))
  ) #get AA seq for each exon (minus separated codons)

#summary of AA sequence for exons of all protein coding transcript
GENCODE.out <- GENCODE_split %>%
  select(geneSym, chr, exonStart, exonEnd, strand, cdsStartStat, cdsEndStat,seq, AAseq) %>%
  mutate(AAlen = nchar(gsub("\\*", "", AAseq)), exonGC = round(stringr::str_count(seq, "G|C")/nchar(seq)*100, 2))

#count AA for each translated exon CDS
AAs <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

GENCODE.AAcomp <- lapply(
  1:length(AAs),
  function(x){
    stringr::str_count(GENCODE.out$AAseq, AAs[x])
  }
)
#combine results into a matrix
GENCODE.AAcomp <- do.call(cbind, GENCODE.AAcomp)
colnames(GENCODE.AAcomp) <- AAs

#add AA counts to summary dataframe
GENCODE.out <- GENCODE.out %>%
  cbind(., GENCODE.AAcomp)

#export file
write.table(GENCODE_out, "./GENCODEv47_comprehensive_exons_aaSeqs.txt", sep = "\t", row.names = FALSE, quote = FALSE)
