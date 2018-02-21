20171219

<!---
  This tests whether we can remove the rows that returned "LRG..." in the 'ensemble_gene_stable_id' column without losing any information.
```{r}
test <- SNP_genes_full[grep("LRG", SNP_genes_full$ensembl_gene_stable_id), ] # creates data frame with just "LRG" rows (161 rows total)
test2 <- SNP_genes_full[-grep("LRG", SNP_genes_full$ensembl_gene_stable_id), ] # creates data frame with no "LRG" rows
setdiff(test$refsnp_id, test2$refsnp_id) # tests whether there are any differences between the refsnp_ids of 'test' and 'test2' (that is, if there are differences, the "LRG" rows would not be duplicates). Returns 'character(0)'
which((test$refsnp_id %in% test2$refsnp_id) == "FALSE") # just to make sure, we test again whether there were any non-matches. Returns 'integer(0), i.e. all refsnp_ids are in both 'test' and 'test2'.
rm(test)
rm(test2)
```
--->
  
  <!-- Again, testing whether we can remove the rows that returned "LRG..." in the 'ensemble_gene_stable_id' column without losing any information.
```{r}
test <- proxy_genes[grep("LRG", proxy_genes$ensembl_gene_stable_id), ] # creates data frame with just "LRG" rows (161 rows total)
test2 <- proxy_genes[-grep("LRG", proxy_genes$ensembl_gene_stable_id), ] # creates data frame with no "LRG" rows
setdiff(test$refsnp_id, test2$refsnp_id) # tests whether there are any differences between the refsnp_ids of 'test' and 'test2' (that is, if there are differences, the "LRG" rows would not be duplicates). Returns 'character(0)'
which((test$refsnp_id %in% test2$refsnp_id) == "FALSE") # just to make sure, we test again whether there were any non-matches. Returns 'integer(0), i.e. all refsnp_ids are in both 'test' and 'test2'.
rm(test)
rm(test2)
```
--->
  
<!-- Trying to iterate for each res

```{r}

for (i in ls(pattern = "res")){
  x <- paste("signif_", i, sep = "")
  assign(x, dplyr::filter(get(i, envir = as.environment(-1)), P.Value < 0.05))
}

#####
pull_signif_genes <- function(x, contrast){
  assign(paste("signif_", contrast, sep = ""), dplyr::filter(x, P.Value < 0.05))
}
pull_signif_genes(x = res1, contrast = "apoe4_for_PET_yes") # doesn't work?


```
-->
  
20171220

<!--
  ```{r, eval = FALSE}
```{r, results = "hide"}
library(VennDiagram)
length(which(!is.na(sig_all2$P.Value.c1)))
length(which(!is.na(sig_all2$P.Value.c2)))
length(which(!is.na(sig_all2$P.Value.c3)))
length(which(!is.na(sig_all2$P.Value.c4)))

length(which(!is.na(sig_all2$P.Value.c1) & !is.na(sig_all2$P.Value.c2)))
length(which(!is.na(sig_all2$P.Value.c1) & !is.na(sig_all2$P.Value.c3)))
length(which(!is.na(sig_all2$P.Value.c1) & !is.na(sig_all2$P.Value.c4)))
length(which(!is.na(sig_all2$P.Value.c2) & !is.na(sig_all2$P.Value.c3)))
length(which(!is.na(sig_all2$P.Value.c2) & !is.na(sig_all2$P.Value.c4)))
length(which(!is.na(sig_all2$P.Value.c3) & !is.na(sig_all2$P.Value.c4)))

#123
length(which(!is.na(sig_all2$P.Value.c1) & !is.na(sig_all2$P.Value.c2) & !is.na(sig_all2$P.Value.c3)))
#234
length(which(!is.na(sig_all2$P.Value.c2) & !is.na(sig_all2$P.Value.c3) & !is.na(sig_all2$P.Value.c4)))
#134
length(which(!is.na(sig_all2$P.Value.c1) & !is.na(sig_all2$P.Value.c3) & !is.na(sig_all2$P.Value.c4)))
#124
length(which(!is.na(sig_all2$P.Value.c1) & !is.na(sig_all2$P.Value.c2) & !is.na(sig_all2$P.Value.c4)))

#1234

length(which(!is.na(sig_all2$P.Value.c1) & !is.na(sig_all2$P.Value.c2) & !is.na(sig_all2$P.Value.c3) & !is.na(sig_all2$P.Value.c4)))
```
-->

20170112

<!--
  Import data and create map file
```{r}
ped <- read.table("C:/Users/bre227/Documents/R/AD_project/Working/genotype_data.txt", header = T)
map <- data.frame(ped[, 1:4], stringsAsFactors = F)
map <- dplyr::select(map, chromosome, snp, position)
map$chromosome <- as.character(map$chromosome)
map$snp <- as.character(map$snp)
# check that there are no unexpected values
unique(map$chromosome)
# noticed that there is an empty row - remove from map and ped file
map <- map[-which(is.na(map$chromosome)), ]
ped <- ped[-which(is.na(ped$chromosome)), ]
# replace chr "X" with "23" as required by format
map$chromosome[map$chromosome == "X"] <- "23"
```

Find duplicated SNP entries in the PED and determine whether all the genotype data is also duplicated (in which case we can delete).
```{r}
length(unique(map$snp)) # only 2,088 out of 2,150 SNPs are unique?
# get duplicated rs IDs from MAP file
map$snp[which(duplicated(map$snp))]
# create new data frame with genotype data from PED file matching those rs IDs
dupe <- ped[ped$snp %in% map$snp[which(duplicated(map$snp))], ]

# By visualising them, we see that an individual's entries for the duplicated SNP can be different. JD suggested to remove them all.

map <- map[!map$snp %in% dupe$snp, ] # remove from map
ped <- ped[!ped$snp %in% dupe$snp, ] # remove from ped
```

Write map file
```{r}
map <- data.frame(lapply(map, function(x){
  gsub(" ", "", x)
}), stringsAsFactors = F)
write.table(map, "C:/Users/bre227/Documents/R/AD_project/Working/snp_data.map", row.names = F, col.names = F, quote = F, sep = "\t")
```

Reformat ped into PED format
```{r}
ped2 <- ped[, -c(1:4)] # remove extraneous columns
meta <- read.table("C:/Users/bre227/Documents/R/AD_project/Working/key_metadata.txt", header = T) # read in metadata
colnames(ped2) <- gsub("X", "", colnames(ped2)) # remove "X" from colnames
ped2 <- data.frame(t(ped2), stringsAsFactors = F) # transpose ped2
ped2$AIBL.Id <- rownames(ped2) # make a new column with the row names (AIBL Ids) to use to bind with metadata
meta$AIBL.Id <- as.character(meta$AIBL.Id) # convert to characters to allow binding
library(dplyr)
ped2 <- left_join(ped2, meta, by = "AIBL.Id")
ped3 <- select(ped2, AIBL.Id, Demographic.Sex, PET, everything())
ped3[, grep("Age|apoe4", colnames(ped3), ignore.case = T)] <- NULL # remove Age and apoe4 columns


# convert male/female to 1/2
ped3$Demographic.Sex <- as.character(ped3$Demographic.Sex)
ped3$Demographic.Sex[ped3$Demographic.Sex == "Male"] <- "1"
ped3$Demographic.Sex[ped3$Demographic.Sex == "Female"] <- "2"

# convert PET status from POS/NEG to 1/0
ped3$PET <- as.character(ped3$PET)
ped3$PET[ped3$PET == "POS"] <- "2"
ped3$PET[ped3$PET == "NEG"] <- "1"
ped3$PET[is.na(ped3$PET)] <- "0"

# clean data
ped4 <- data.frame(lapply(ped3, function(x){
  gsub(" ", "", x)
  gsub("_", " ", x)
}), stringsAsFactors = F)
ped4[is.na(ped4)] <- "0 0"

# write file
write.table(ped4, "C:/Users/bre227/Documents/R/AD_project/Working/snp_data.ped", row.names = F, col.names = F, quote = F, sep = "\t")
```
-->
  
20180116

<!--
  See if we can get more annotation information by using the locations of the SNPs to get the gene names, rather than the ensembl IDs.

```{r}
library(GenomicRanges)
# include "chr" prefix before each chromosome name
cnv <- SNP_genes
cnv$chr_name <- paste("chr", cnv$chr_name, sep = "")

# make GRanges object for SNPs
cnv1 <- makeGRangesFromDataFrame(cnv,
                                 keep.extra.columns = F,
                                 seqnames.field = "chr_name",
                                 start.field = "chrom_start",
                                 end.field = "chrom_end")

# create function to match overlaps
#splitColumnByOverlap <- function(query, subject, column = "ENTREZID", ...){
#    olaps <- findOverlaps(query, subject, ...)
#    f1 <- factor(subjectHits(olaps),
#                 levels=seq_len(subjectLength(olaps)))
#    splitAsList(mcols(query)[[column]][queryHits(olaps)], f1)
#}

# load library
library(Homo.sapiens)
gns <- genes(Homo.sapiens, columns = c("ENTREZID", "SYMBOL", "ENSEMBL"))
head(gns)

# run function
ols <- findOverlaps(gns, cnv1, select = "all")
ols
unique(queryHits(ols)) # only returns hits for 227 SNPs??

ols1 <- subsetByOverlaps(gns, cnv1)
ols1 # confirmed?

ols2 <- data.frame(ols1)

```
-->
  
  