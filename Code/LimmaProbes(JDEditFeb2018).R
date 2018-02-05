
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)

ll <- lapply(cluster_genexp, function (x) {

xm <- as.matrix(log(x[,-1]))
svd.decomp2 <- svd(xm)
svd.v <- svd.decomp2$v
X.pca <- xm %*% svd.v
xout <-  X.pca[,1]
return(xout)
})

lld <- do.call(rbind.data.frame,ll)
lld2 <- t(lld)
rownames(lld2) <- cluster_genexp[[1]][,1]
colnames(lld2) <- names(cluster_genexp)

genemap <- read.csv("~/Dropbox/eQTL/Results/cluster_data/gene_clusters_hsa_map.csv")
keymeta <- read.delim("~/Dropbox/transfer/Jan2018/AAIC Abstracts/key_metadata.txt",stringsAsFactors = F)

lld3 <- lld2[as.character(keymeta$AIBL.Id),]
kp <- ifelse(is.na(keymeta$PET)==T,FALSE,TRUE)
keymeta1 <- keymeta[kp,]

 Y1 <- factor(as.character(keymeta1$PET))
 lld4 <- lld3[kp,]
 treat1<-factor(c(Y1), labels=c("Normal","Cancer"))
 design_treat1<-model.matrix(~0+treat1)
 colnames(design_treat1) <- c("Normal","Cancer")
 contrast.matrix_treat1<-makeContrasts(Normal-Cancer, levels=design_treat1)
## run LIMMA
# Transpose data into matrix for limma fit
 x1 <- t(lld4)
 fit.lm1<-lmFit(x1,design_treat1)   #using the matrix x1 from above
 fit1.c <-contrasts.fit(fit.lm1,contrast.matrix_treat1)
 fit.eb1<- eBayes(fit1.c)
## comparison treated - untreated
topALL1<-toptable(fit=fit1.c, eb=fit.eb1, adjust="fdr", n=244, genelist=colnames(lld4))
topALL.1 <- topALL1[,c("ID","P.Value")]

## Merge in meta data ##
names(genemap)[1]="ID"
res.NetW <- merge(topALL.1,genemap,by.id="ID")

### Set up the LIMMA for the genes alone and combine results ###

cluster_genexp2 <- lapply(cluster_genexp, function(x) x[,-1])
fulldata <- do.call(cbind.data.frame,cluster_genexp2)

cn <- unlist(strsplit(x=colnames(fulldata),split = ".X"))[seq(2,length(unlist(strsplit(x=colnames(fulldata),split = ".X"))),2)]
colnames(fulldata)=cn

### This line does not work
fulldata1 <- fulldata[,as.character(genemap$root_gene_probe_id)]

length(intersect(as.character(genemap$root_gene_probe_id),colnames(fulldata)))