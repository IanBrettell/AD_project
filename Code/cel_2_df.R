##########################################
By James Doecke, CSIRO
##########################################

### Bring in cel files and create eset ###
rm(list=ls())
setwd("~/Documents/Work/AIBLData/AIBL_Expression_Arrays/All_Arrays_B1_B6/")
source("http://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("oligo")
biocLite("limma")
biocLite("pd.huex.1.0.st.v2")
biocLite("huex10sttranscriptcluster.db")

# load the affy library
library(oligo)
library(pd.huex.1.0.st.v2)
library(huex10sttranscriptcluster.db)

# Read in the CEL files in the directory
celFiles <- list.celfiles()[1:5]
affyRaw <- read.celfiles(celFiles)

# You might need to install and load a package for the specific array you are using
# It may try to load it automatically, but may fail.  Install & load the library manually if this happens.
eset <- rma(affyRaw)

dim(eset)
# Features  Samples 
#    22011      218 
   
# Finally, save the data to an output file to be used by other programs, etc (Data will be log2 transformed and normalized)
write.exprs(eset,file="data.txt")

data <- exprs(eset)

# http://www.affymetrix.com/support/help/exon_glossary/index.affx#rnaclass
# get the biological annotation via the getNetAffx function
featureData(eset) <- getNetAffx(eset, "probeset")

# get overview
eset

# The biological annotation made available through NetAffx is now stored in the feature- Data slot \
featureData(eset)
varLabels(featureData(eset))

# the gene assignment for the 2 first metaprobesets (core) can be obtained as follows
pData(featureData(eset))[1:2, "geneassignment"]
	
# get sample names
sns <- sampleNames(eset)


### pull colnames and convert to AIBL.IDs ###
colnames(data)[20] <- "BW_CSIRO_125_(HuEx-1_0-st-v2)2.CEL"
colnames(data)[142] <- "BW_CSIRO_71_(HuEx-1_0-st-v2)2.CEL"
s1 <- colnames(data)[1:147]
s2 <- colnames(data)[148:length(colnames(data))]
samps1 <- unlist(strsplit(s1,split="_"))[seq(3,length(unlist(strsplit(s1,split="_"))),5)]
samps2 <- gsub("CSIRO","",unlist(strsplit(s2,split="_"))[seq(2,length(unlist(strsplit(s2,split="_"))),4)])
samples <- c(samps1,samps2)
colnames(data) <- samples

# plot the log ratio's vs intensities
xl <- c(2.8, 4)
yl <- c(-1, 1)
MAplot(eset[, 1:3], pairs=TRUE, ylim=yl, xlim=xl)

### pull gene symbols and match to probe ids to replace rownames ###
x <- huex10sttranscriptclusterSYMBOL
# Get the probe identifiers that are mapped to a gene symbol
mapped_probes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_probes])
if(length(xx) > 0) {
# Get the SYMBOL for the first five probes
xx[1:5]
# Get the first one
xx[[1]]
}

xx1 <- as.matrix(unlist(xx))

gns <- rownames(data)
gnsl <- xx1[rownames(xx1) %in% gns]
length(gnsl)
#[1] 15232


### Find ENSEMBL IDs ###
x3 <- huex10sttranscriptclusterENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID
mapped_genes <- mappedkeys(x3)
# Convert to a list
xx3 <- as.list(x3[mapped_genes])

xx3.1 <- as.matrix(unlist(xx3))
ens <- xx3.1[rownames(xx3.1) %in% gns]
length(ens)
#[1] 14497

### Find Entrez IDs ###
x4 <- huex10sttranscriptclusterENTREZID
# Get the probe identifiers that are mapped to an ENTREZ Gene ID
mapped_probes <- mappedkeys(x4)
# Convert to a list
xx4 <- as.list(x4[mapped_probes])

xx4.1 <- as.matrix(unlist(xx4))
entz <- xx4.1[rownames(xx4.1) %in% gns]
length(entz)
#[1] 15232

### Read the ID mapping file ###
idmap <- read.csv("~/Dropbox/transfer/Dec2017/All.Prep.recoveriesFor Bill..csv")

idmap1 <- idmap[idmap$AIBL.ID!=0 & idmap$CSIRO..!=0,c(1,2)]
idmap2 <- idmap1[idmap1$CSIRO..%in%colnames(data),]

data2 <- data[,as.character(idmap2$CSIRO..)]
colnames(data2) <- idmap2$AIBL.ID

write.table(data2,row.names = T,col.names = T,"~/Dropbox/eQTL/Data/AIBL_expression_set/AIBL_Gene_Expression_UpdtdDec2017.txt",sep = "")
write.table(colnames(data2),row.names = T,col.names = T,"~/Dropbox/eQTL/Data/AIBL_expression_set/AIBL_Gene_Expression_IDs_UpdtdDec2017.txt")

