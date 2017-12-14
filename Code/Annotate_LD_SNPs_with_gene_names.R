##################################################
# Author: Ian Brettell
# Date: 20171208
# Purpose: To determine whether the SNPs that are in LD with those 
#           targeted fall within (or without) the same genes
##################################################

##################################################
# IMPORT FILES
##################################################

setwd("~/R/AD Project/Data/CSIROData/SNP/Expanded SNP set")
AIBLSNPlist <- read.delim("AIBLgene_SNP_LIST_04032015.tsv")
SNAP_edit <- read.delim("SNAPResults_edit.txt")
SNAP <- read.delim("SNAPResults.txt")

##################################################
# INSTALL PACKAGES
##################################################

source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("biomaRt")
library(biomaRt)
library(tidyverse)
biocLite("GenomicRanges")
library(GenomicRanges)
biocLite("Homo.sapiens")
library(Homo.sapiens)
biocLite("stephenturner/annotables")
library(annotables)
biocLite("ensemblVEP")
library(ensemblVEP)

##################################################
# COMBINE DATA SETS
##################################################

# note: 'SNAP_edit' has removed all SNPs from 'SNAP' that did not have proxies  

length(unique(SNAP_edit$SNP)) # number of SNPs = from 2,084 reduced to 1527 (557 removed)
length(unique(AIBLSNPlist$GENE)) # check number of genes = 298
which(AIBLSNPlist$GENE == "") # shows that there are 16 SNPs with no associated gene
which(AIBLSNPlist$GENE == "-") # shows that there is 1 additional SNP with no associated gene
length(setdiff(AIBLSNPlist$SNP, SNAP_edit$SNP)) # 557: correct

AIBLdummy <- dplyr::select(AIBLSNPlist, SNP, GENE) # creates dummy data frame with just two variables
SNAP2 <- left_join(SNAP_edit, AIBLdummy, by = "SNP") # combine data frames
rm(AIBLdummy) # remove dummy
SNAP2 <- select(SNAP2, SNP, GENE, everything()) # reorder columns

##################################################
# GET ANNOTATION FOR PROXY SNPS
##################################################

# create vector with unique SNPs
uq_proxies <- as.vector(unique(SNAP2$Proxy)) 

# access SNP mart
listMarts() # lists marts available
mart <- useMart('ENSEMBL_MART_SNP')
listDatasets(mart) # lists datasets available
mart <- useMart('ENSEMBL_MART_SNP', 'hsapiens_snp')

# get loci for SNPs
proxy_loci <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"),
                    filters = c("snp_filter"),
                    values = c(uq_proxies),
                    mart = mart)

##################################################

#create vectors for start and end loci
#proxy_start <- as.vector(proxy_loci$chrom_start)
#proxy_end <- as.vector(proxy_loci$chrom_end)

# access gene mart
#listMarts()
#gene_mart <- useMart('ENSEMBL_MART_ENSEMBL')
#listDatasets(gene_mart)
#gene_mart <- useMart('ENSEMBL_MART_ENSEMBL', 'hsapiens_gene_ensembl')

# get genes associated with proxy SNPs
#proxy_gns <- getBM(attributes = c("start_position", "hgnc_symbol"),
                   #filters = c("start", "end"),
                   #values = list(proxy_start, proxy_end),
                   #mart = gene_mart) # query: command takes four hours to run?

## output was 49 million observations of 2 variables??

##################################################
# Alternative method to the one above

genelist <- grch38 %>% 
  dplyr::select(symbol, chr, start, end)
gns <- makeGRangesFromDataFrame(genelist,
                                keep.extra.columns = T,
                                ignore.strand = T,
                                seqnames.field = "chr",
                                start.field = "start",
                                end.field = "end")

proxygr <- makeGRangesFromDataFrame(proxy_loci,
                                    ignore.strand = T,
                                    seqnames.field = "chr_name",
                                    start.field = "chrom_start",
                                    end.field = "chrom_end")

overlapGenes <- data.frame(findOverlaps(proxygr, gns))
regionsWIthHits <- data.frame(proxygr[overlapGenes$queryHits])
regionsWIthHits$genes <- gns$symbol[overlapGenes$subjectHits]
length(unique(regionsWIthHits$genes)) # returns 734
##################################################

##################################################

## test using online BioMart

#mart <- useMart('ENSEMBL_MART_SNP', 'hsapiens_snp')
#proxy_loci <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "chrom_strand"),
                    #filters = c("snp_filter"),
                    #values = c(uq_proxies),
                    #mart = mart) # showed that all ~18,000 SNPs lie on the forward strand??

#proxies_500 <- head(SNAP2, n = 500) %>% 
#  dplyr::select(Proxy, Chromosome, Coordinate_HG18)
#colnames(proxies_500)[colnames(proxies_500) == "Coordinate_HG18"] <- "Start"
#proxies_500$Strand <- 1
#chr <- as.character(proxies_500$Chromosome)
#chr <- gsub("chr", "", chr)
#proxies_500$Chromosome <- chr
#rm(chr)
#proxies_500$loci <- paste(proxies_500$Chromosome, proxies_500$Start, proxies_500$Start, "1", sep = ":")
#write.table(proxies_500$loci, file = "biomart.test.txt", quote = F, row.names = F, col.names = F)

## only returned a handfull of entries, instead of the expectd 500

##################################################

# test using ensemblVEP package

library(ensemblVEP)
write.table(uq_proxies, "~/R/AD Project/Working/uq_proxies.txt", quote = F, row.names = F, col.names = F)
## submit file to "http://www.ensembl.org/Tools/VEP" and use output as follows

uq_proxies_output <- read.delim("~/R/AD Project/Working/uq_proxies_VEP_output.txt", header = T)

# extract unique SNPid and gene names
prxy_snpsandgns <- dplyr::select(uq_proxies_output, X.Uploaded_variation, SYMBOL)
uq_prxy_snpsandgns <- prxy_snpsandgns[!duplicated(prxy_snpsandgns), ]

# produce same file but for AIBL SNPs

AIBL_SNPs <- as.vector(unique(AIBLSNPlist$SNP))
write.table(AIBL_SNPs, "~/R/AD Project/Working/AIBL_SNPs.txt", quote = F, row.names = F, col.names = F)

# extract unique SNPid and gene names

AIBLSNPs_output <- read.delim("~/R/AD Project/Working/AIBL_SNPs_VEP_output.txt", header = T)
AIBL_snpsandgns <- dplyr::select(AIBLSNPs_output, X.Uploaded_variation, SYMBOL)
uq_AIBL_snpsandgns <- AIBL_snpsandgns[!duplicated(AIBL_snpsandgns), ]

########################

# through bioconductor
library(ensemblVEP)
myparam <- VEPFlags(version = max(unlist(currentVEP())),
                    scriptPath = character(),
                    flags <- list(vcf = F, 
                                 everything,
                                 output_file = "~/R/AD Project/Working/AIBL_SNPs.txt"))
