####this script helps to collate the individual HPA immune cell type expression specificity files
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

setwd("C:/data/TissueAbSplicing")

####collate enriched/enhanced/nonspecific genes from Human Protein Atlas for immune cell lineages (Human Blood Atlas to be specific)
#generate list of files + metadata
files <- list.files("./HPAblood", pattern = "(enriched\\.tsv)|(enhanced\\.tsv)|(nonspecific\\.tsv)|(nonspecificOtherEnhanced\\.tsv)", full.names = TRUE)

metadata <- data.frame("file" = files) %>%
  mutate(sample = gsub("(\\.\\/HPAblood\\/)|(\\.tsv)", "", file)) %>%
  separate(col = sample, sep = "_", into = c("tissue","expression"), remove=FALSE)

#read & process files (keeping required columns)
tempList <- vector("list", length = nrow(metadata))

for(i in 1:nrow(metadata)){
  tempList[[i]] <- read.table(metadata$file[i], sep = "\t", header = TRUE) %>%
    select(Gene, Ensembl, Uniprot, Protein.class, RNA.tissue.specificity, RNA.tissue.distribution)
}
names(tempList) <- metadata$sample

#count number of genes in each file and append to metadata object
metadata$n.genes <- unlist(lapply(X = tempList, FUN = nrow))

####checking overlaps between group-enriched genes of the immune cell types

allTissues <- unique(metadata$tissue) #all immune cell types in data

#pull out all group-enriched genes for each immune cell type
mergedTissues <- vector("list", length = length(allTissues))
for(i in 1:length(allTissues)){
  mergedTissues[[i]] <- tempList[which(grepl(paste0("^",allTissues[i], "_groupenriched"), names(tempList)))] %>% do.call(rbind, .)
}
names(mergedTissues) <- allTissues

list.genes <- lapply(
  mergedTissues, function(x){
    x$Ensembl}
  )

#calculate group-enriched gene pair-wise overlaps between cell types
overlap.list <- lapply(
  list.genes,
  function(refSet){
    lapply(list.genes, function(comparedSet){sum(comparedSet %in% refSet)})
  }
)

overlap.df <- data.frame(
  "refSet" = stringr::str_split_i(names(unlist(overlap.list)), pattern = "\\.", 1),
  "comparedSet" = stringr::str_split_i(names(unlist(overlap.list)), pattern = "\\.", 2),
  "n.overlap" = unlist(overlap.list)
) %>%
  group_by(refSet) %>%
  mutate(pct.overlap = round((n.overlap/max(n.overlap)*100),2)) %>% # %overlap of refSet in comparedSet
  ungroup() %>%
  filter(pct.overlap != 100) #drop self-comparisons

#overlap% matrix
overlap.df.wide <- overlap.df %>%
  select(-3) %>%
  pivot_wider(names_from = comparedSet, values_from = pct.overlap)

####output cell type + expression specificity matrix

tempList.comb <- tempList

#combine '_nonspecific' with '_nonspecificOtherEnhanced' (essentially the same category)
#assumes that '_nonspecific' is always the element before '_nonspecificOtherEnhanced' (sorted earlier)
for(i in 1:length(tempList.comb)){
  if(grepl("_nonspecific$", names(tempList.comb)[i])){
    tempList.comb[[i]] <- rbind(tempList.comb[[i]], tempList.comb[[i+1]])
  } else if (grepl("_nonspecificOtherEnhanced$", names(tempList.comb)[i])){
    tempList.comb[[i]] <- NA
  }
}
tempList.comb <- tempList.comb[which(lengths(tempList.comb) != 1)] #drop now-empty groups

#keep only geneIDs & mark presence in a cell type expression subset with 1
tempListGenes <- lapply(tempList.comb, FUN = function(x){x %>% mutate(present = 1) %>% select(geneID = Ensembl, present)})
#all genes in this dataset
allGenes <- do.call(rbind, tempListGenes)$geneID
allGenes <- allGenes[order(allGenes)]
allGenes <- unique(allGenes[allGenes != ""])

#collate the speificity matrix
res.out <- c(list(data.frame(geneID = allGenes)), tempListGenes) %>%
  purrr::reduce(., .f = left_join, by = "geneID") %>%
  mutate(across(.cols = -1, .fns = function(x){ifelse(is.na(x), FALSE, TRUE)}))
colnames(res.out) <- c("geneID", names(tempListGenes))

#save matrix
write.table(res.out, "./HPA_bloodSpecificityMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)
