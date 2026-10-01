library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(tidyr)
library(stringr)

setwd("C:/DataProcessing/Database/Gencode")

####read GTF file####
#downloaded from GENCODE website (GENCODE v47, primary, comprehensive)
gtf <- read.table("gencode.v47.primary_assembly.annotation.gtf", sep = "\t")

####retrieve exon/transcript/gene attributes####
#this function retrieves information from the attribute (9th) column of GTF files
retrieveAttr <- function(x, regex){
  matches <- regexpr(regex, x, perl = TRUE)
  strings <- rep(NA, length(x))
  strings[matches != -1] <- regmatches(x, matches)
  return(strings)
}

exons <- gtf %>%
  filter(V3 == "exon") %>% #keep only 'exon' rows
  mutate(
    V4 = V4-1,
    "geneSym" = retrieveAttr(V9, "(?<=gene_name ).*?(?=;)"),
    "geneID" = gsub("\\..*?$", "", retrieveAttr(V9, "(?<=gene_id ).*?(?=;)")),
    "HGNC.ID" = retrieveAttr(V9, "(?<=hgnc_id ).*?(?=;)"),
    "txID" = gsub("\\..*?$", "", retrieveAttr(V9, "(?<=transcript_id ).*?(?=;)")),
    "txName" = retrieveAttr(V9, "(?<=transcript_name ).*?(?=;)"),
    "exID" = retrieveAttr(V9, "(?<=exon_id ).*?(?=;)"),
    "exNo" = retrieveAttr(V9, "(?<=exon_number ).*?(?=;)"),
    "geneType" = retrieveAttr(V9, "(?<=gene_type ).*?(?=;)"),
    "txType" = retrieveAttr(V9, "(?<=transcript_type ).*?(?=;)"),
    "gencodeLevel" = retrieveAttr(V9, "(?<= level ).*?(?=;)"),
    "txSupportLevel" = retrieveAttr(V9, "(?<=transcript_support_level ).*?(?=;)"),
    "tags" = gsub("(tag )|(; havana.*?$)", "", retrieveAttr(V9, "(?<=tag ).*(?=;)"))
  ) %>%
  group_by(geneID) %>%
  mutate(earliestCoord = min(V4)) %>% #find the earliest genomic starting point for each gene (5' end for + strand, 3' end for - strand)
  mutate(exNo = as.numeric(exNo)) %>%
  arrange(V1, earliestCoord, txID, exNo) %>% #sort by transcripts, in ascending order based on gene-wise exon numbering
  group_by(txID) %>%
  mutate(totalExNo = max(exNo)) %>% #get total number of exons in transcript
  ungroup() %>%
  select(
    geneSym, geneID, HGNC.ID, txName, txID, exID, chr = V1, exStart = V4, exEnd = V5, strand = V7, exNo, totalExNo, geneType, txType, gencodeLevel, txSupportLevel, tags
  ) %>% #rename columns accordingly
  mutate(
    exon5pProblem = grepl("mRNA_start_NF", tags) & exNo == 1,
    exon3pProblem = grepl("mRNA_end_NF", tags) & exNo == totalExNo
  ) %>% #to check for incomplete transcripts (starts or ends not known)
  filter(
    chr %in% paste0("chr",c(1:23,"X","Y"))
  ) %>% #keep only the main chromosomes
  group_by(txID) %>%
  mutate(
    coords = paste0(chr, ":", exStart, "-", exEnd, ":", strand),
    intronCoords.upst = paste0(chr, ":", ifelse(strand == "+", lag(exEnd), exEnd), "-",ifelse(strand == "+", exStart, lag(exStart)), ":", strand),
    intronCoords.downst = paste0(chr, ":", ifelse(strand == "+", exEnd, lead(exEnd)), "-",ifelse(strand == "+", lead(exStart),exStart), ":", strand)
  ) %>% #coordinates of exons and flanking introns (upstream and downstream)
  mutate(
    intronCoords.upst = ifelse(exNo == 1, NA, intronCoords.upst),
    intronCoords.downst = ifelse(exNo == totalExNo, NA, intronCoords.downst)
  ) %>% #set false introns to NA
  ungroup()

####calculating exon physical properties####

hash4 <- function(X){
  hash4_total <- 0
  for (i in 1:nchar(X[is.na(X) == FALSE])[1]){
    hash4_total <- hash4_total + as.numeric(substr(chartr("ACGT", "0123", X), i, i)) * 4^(nchar(X[is.na(X) == FALSE])[1]-i)
  }
  return(hash4_total)
} #part of MaxEnt function

calcMAXENT_3SS <- function(X, me2x3){
  any.na <- any(is.na(X))
  if(any.na == TRUE){
    num_rows <- length(X)
    filled_rows <- which(!is.na(X))
    
    X <- X[!is.na(X)]
  }
  
  wt_3SS <- data.frame(Base = c("A", "C", "G", "T"),
                       bgd = c(0.27, 0.23, 0.23, 0.27),
                       N1 = c(0.9903, 0.0032, 0.0034, 0.0030),
                       N2 = c(0.0027, 0.0037, 0.9905, 0.0030)
  )
  
  scores <- log2(
    ((wt_3SS$N1[match(substr(X, 19, 19), wt_3SS$Base)] * wt_3SS$N2[match(substr(X, 20, 20), wt_3SS$Base)]) /
       (wt_3SS$bgd[match(substr(X, 19, 19), wt_3SS$Base)] * wt_3SS$bgd[match(substr(X, 20, 20), wt_3SS$Base)])) *
      ((me2x3[match(hash4(substr(X, 1, 7))+1, me2x3$Row), "acc1"] *
          me2x3[match(hash4(substr(X, 8, 14))+1, me2x3$Row), "acc2"] *
          me2x3[match(hash4(paste0(substr(X, 15, 18), substr(X, 21, 23)))+1, me2x3$Row), "acc3"] * 
          me2x3[match(hash4(substr(X, 5, 11))+1, me2x3$Row), "acc4"] *
          me2x3[match(hash4(substr(X, 12, 18))+1, me2x3$Row), "acc5"]) / 
         (me2x3[match(hash4(substr(X, 5, 7))+1, me2x3$Row), "acc6"] *
            me2x3[match(hash4(substr(X, 8, 11))+1, me2x3$Row), "acc7"] *
            me2x3[match(hash4(substr(X, 12, 14))+1, me2x3$Row), "acc8"] *
            me2x3[match(hash4(substr(X, 15, 18))+1, me2x3$Row), "acc9"]))
  )
  
  if(any.na == TRUE){
    X <- as.vector(matrix(NA, nrow = num_rows))
    X[filled_rows] <- scores
    return(X)
  } else {
    return(scores)
  }
}  #MaxEntScan implemented in R

calcMAXENT_5SS <- function(X, me2x5){
  any.na <- any(is.na(X))
  if(any.na == TRUE){
    num_rows <- length(X)
    filled_rows <- which(!is.na(X))
    
    X <- X[!is.na(X)]
  }
  
  wt_5SS <- data.frame(Base = c("A", "C", "G", "T"),
                       bgd = c(0.27, 0.23, 0.23, 0.27),
                       N1 = c(0.004, 0.0032, 0.9896, 0.0032),
                       N2 = c(0.0034, 0.0039, 0.0042, 0.9884)
  )
  
  scores <- log2(
    (wt_5SS$N1[match(substr(X, 4, 4), wt_5SS$Base)] * wt_5SS$N2[match(substr(X, 5, 5), wt_5SS$Base)]) /
      (wt_5SS$bgd[match(substr(X, 4, 4), wt_5SS$Base)] * wt_5SS$bgd[match(substr(X, 5, 5), wt_5SS$Base)]) *
      (me2x5$Value[match(paste0(substr(X, 1, 3), substr(X, 6, 9)), me2x5$Seq)])
  )
  
  if(any.na == TRUE){
    X <- as.vector(matrix(NA, nrow = num_rows))
    X[filled_rows] <- scores
    return(X)
  } else {
    return(scores)
  }
}  #MaxEntScan implemented in R

#define genomic sequence & read files required for MaxEnt calculations
genome <- BSgenome.Hsapiens.UCSC.hg38
me2x5 <- read.table("./me2x5.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE) #for calcMAXENT_3SS
me2x3 <- read.table("./me2x3.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE) #for calcMAXENT_5SS

#calculate splice site coordinates
SSranges <- exons %>%
  ungroup() %>%
  mutate(
    fiveSS_start = ifelse(strand == "+", exEnd - 2, exStart + 1 - 6),
    fiveSS_end = ifelse(strand == "+", exEnd + 6, exStart + 1 + 2),
    threeSS_start = ifelse(strand == "+", exStart + 1 - 20, exEnd - 2),
    threeSS_end = ifelse(strand == "+", exStart + 1 + 2, exEnd + 20),
    index = seq.int(1,nrow(.),1)
  )

#retrieve splice site sequences
GR_3SS <- makeGRangesFromDataFrame(SSranges %>% select(chr, start = threeSS_start, end = threeSS_end, strand, index), keep.extra.columns = TRUE)
seq3SS <- as.character(getSeq(genome, GR_3SS))

GR_5SS <- makeGRangesFromDataFrame(SSranges %>% select(chr, start = fiveSS_start, end = fiveSS_end, strand, index), keep.extra.columns = TRUE)
seq5SS <- as.character(getSeq(genome, GR_5SS))

#calculate intronic GC
intronCoords <- data.frame(intronCoords = unique(c(exons$intronCoords.upst, exons$intronCoords.downst))[-1])
gr.intronCoords <- GRanges(intronCoords$intronCoords)
start(gr.intronCoords) <- start(gr.intronCoords)+1 #convert from 0-based to 1-based coordinates

seqIntrons <- as.character(getSeq(genome, gr.intronCoords))

intronCoords$intGC <- round(str_count(seqIntrons, "G|C") / nchar(seqIntrons), 4) #calculaate GC%

#calculate exonic GC
exonCoords <- data.frame(exonCoords = unique(exons$coords))
gr.exonCoords <- GRanges(exonCoords$exonCoords)
start(gr.exonCoords) <- start(gr.exonCoords)+1 #convert from 0-based to 1-based coordinates

seqExons <- as.character(getSeq(genome, gr.exonCoords))

exonCoords$exGC <- round(str_count(seqExons, "G|C") / nchar(seqExons), 4) #calculaate GC%

#add physical parameters to main dataframe
exons <- exons %>%
  group_by(txID) %>%
  mutate(
    exLen = exEnd - exStart,
    upIntLen = ifelse(strand == "+", exStart - lag(exEnd), lag(exStart) - exEnd),
    downIntLen = ifelse(strand == "+", lead(exStart) - exEnd, exStart - lead(exEnd))
  ) %>%
  ungroup() %>%
  mutate(
    upIntGC = intronCoords$intGC[match(.$intronCoords.upst, intronCoords$intronCoords)],
    exGC = exonCoords$exGC[match(.$coords, exonCoords$exonCoords)],
    downIntGC = intronCoords$intGC[match(.$intronCoords.downst, intronCoords$intronCoords)]
  ) %>%
  mutate(
    seq3SS = ifelse(is.na(upIntLen), NA, seq3SS),
    seq5SS = ifelse(is.na(downIntLen), NA, seq5SS)
  ) %>%
  mutate(
    MAXENT_3SS = calcMAXENT_3SS(seq3SS, me2x3),
    MAXENT_5SS = calcMAXENT_5SS(seq5SS, me2x5)
  )

####define alternative splicing status of exons####

##defining cassette exon events
#make genomicRange object of all exons
GR_exons <- makeGRangesFromDataFrame(
  exons %>%
    mutate(
      exPos = case_when(exNo == 1 ~ "first", exNo == totalExNo ~ "last", TRUE ~ "int") #mark exon position (first, last, internal)
    ) %>%
    select(geneID, exStart, exEnd, coords, txID, exPos),
  keep.extra.columns = TRUE, seqnames.field = "geneID"
)
#find number of range overlaps between all exons
overlapExons <- as.data.frame(findOverlaps(GR_exons, GR_exons, type=c("any"), select=c("all")))
#make genomicRange object of all transcripts
GR_tx <- makeGRangesFromDataFrame(
  exons %>%
    group_by(geneID, txID) %>%
    summarise(start = min(exStart), end = max(exEnd)) %>%
    select(geneID, start, end, txID),
  keep.extra.columns = TRUE, seqnames.field = "geneID"
)
#find number of range overlaps between exons and full transripts
overlapTx <- as.data.frame(findOverlaps(GR_exons, GR_tx, type=c("any"), select=c("all")))
#count number of transcripts that overlap each exon (i.e. exonic genomic sequence is physically within transcript start/stop genomic positions)
txCount <- overlapTx %>%
  group_by(queryHits) %>%
  arrange(queryHits) %>%
  summarise(n.tx = n()) %>%
  cbind(., as.data.frame(GR_exons)[.$queryHits,] %>% select(geneID = seqnames, coords)) %>%
  select(-queryHits) %>%
  distinct()
#retrieve exon/transcript information for all exon-exon overlap-pairs
overlapCount <- cbind(
  as.data.frame(GR_exons) %>% slice(overlapExons$queryHits) %>% select(geneID = seqnames, coords.query = coords, exPos.query = exPos),
  as.data.frame(GR_exons) %>% slice(overlapExons$subjectHits) %>% select(coords.subject = coords, exPos.subject = exPos , txID.subject = txID)
) %>%
  distinct() %>%
  mutate(terminal.subject = exPos.subject %in% c("first","last")) %>%
  left_join(., txCount, by = c("coords.query" = "coords", "geneID" = "geneID"))

cassetteExons <- overlapCount %>%
  arrange(geneID, coords.query, txID.subject, terminal.subject) %>%
  filter(!duplicated(paste0(coords.query, "_", txID.subject))) %>%
  group_by(coords.query, n.tx) %>%
  arrange(coords.query) %>%
  mutate(
    n.ov = n(),
    n.terminalMatches = sum(terminal.subject)
    ) %>%
  select(geneID, exPos.query, coords = coords.query, n.ov, n.terminalMatches, n.tx) %>%
  distinct() %>% filter(exPos.query == "int")
#n.tx: number of transcripts that this exon is found in
#n.ov: number of exons from different transcripts that overlap with a particular exon
#cassette exons: internal exons with a number of overlaps lower than number of transcripts in which the exon can be found in (minus any matches to terminal exons)
cassetteExoncoords <- cassetteExons %>% filter((n.ov < (n.tx-n.terminalMatches)) & n.tx != 1) %>% pull(coords)

exons <- exons %>% mutate(SE = coords %in% cassetteExoncoords) #annotate cassette exons (SE / skippable exons)

####Define exon type (w.r.t. to functional properties)####
#read additional transcript annotations (CDS start/end coordinates)
#download from UCSC table browser (GENCODEv47 > knownGeneV47 - all columns)
addtlAnnotations <- read.table(gzfile("gencodeV47.comprehensive.txAnnotations.gz"), sep = "\t", comment.char = '&', header = TRUE)
colnames(addtlAnnotations) <- gsub(".*\\..*\\.", "", colnames(addtlAnnotations))
addtlAnnotations <- addtlAnnotations %>%
  select(txID = name, txStart, txEnd, cdsStart, cdsEnd) %>%
  mutate(txID = gsub("\\..*", "", txID))

#add transcript CDS details to main dataframe
exons <- exons %>%
  left_join(., addtlAnnotations, by = "txID")
#define exon types based on protein coding annotations
exons <- exons  %>%
  mutate(
    exType = ifelse(
      strand == "+",
      case_when(
        cdsStart == cdsEnd ~ "NCE", #non-coding exon
        totalExNo == 1 ~ "SET", #single-exon / intronless transcript
        exEnd < cdsStart ~ "5UTR",  #5'UTR
        exStart > cdsEnd ~ "3UTR",  #3'UTR
        exStart <= cdsStart & exEnd >= cdsEnd ~ "SSE", #translation start+stop exon
        exStart <= cdsStart & exEnd >= cdsStart ~ "FCE", #first coding exon
        exStart <= cdsEnd & exEnd >= cdsEnd ~ "LCE", #last coding exon
        exStart > cdsStart & exEnd < cdsEnd ~ "ICE", #internal coding exon
        TRUE ~ "NCE"
      ),
      case_when(
        cdsStart == cdsEnd ~ "NCE",
        totalExNo == 1 ~ "SET",
        exStart > cdsEnd ~ "5UTR",
        exEnd < cdsStart ~ "3UTR",
        exStart <= cdsStart & exEnd >= cdsEnd ~ "SSE",
        exStart <= cdsEnd & exEnd >= cdsEnd ~ "FCE",
        exStart <= cdsStart & exEnd >= cdsStart ~ "LCE",
        exStart > cdsStart & exEnd < cdsEnd ~ "ICE",
        TRUE ~ "NCE"
      )
    )
  )

####add UCSC-based transcript ranks (used when non-redundant exons by coordinates are required; properties calculated in this script will be based on the top-ranked transcript)
#download from UCSC table browser (All GENCODEv47 > comprehensive (WgEncodeGencodeCompV47) - all columns)
txRanks <-  read.table(gzfile("gencodeV47.comprehensive.txRanks.gz"), sep = "\t", header = TRUE, comment.char = '&') %>%
  mutate(transcriptId = gsub("\\..*?$", "", transcriptId)) %>%
  select(txID = transcriptId, txRank = transcriptRank)
#sort transcripts by rank and exons by exon number
exons <- exons %>%
  left_join(., txRanks, by = "txID") %>%
  arrange(geneSym, txRank, exNo)

####export file to use as reference file in other processes
write.table(exons, "./GENCODEv47_comprehensive_exons_v2.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#this version only keeps one copy of each non-rendundant exon (smaller file that is useful when full transcript-wise information is not necessary)
write.table(exons %>% filter(!duplicated(coords)), "./GENCODEv47_comprehensive_exons_concise_v2.txt", sep = "\t", row.names = FALSE, quote = FALSE)
