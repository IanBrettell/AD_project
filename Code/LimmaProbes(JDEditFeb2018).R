
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)

load("C:/Users/bre227/Dropbox/eQTL/Results/cluster_data/gene_clusters_expression.RData")
ll <- lapply(cluster_genexp, function (x) {

xm <- as.matrix(log(x[,-1]))
svd.decomp2 <- svd(xm)
svd.v <- svd.decomp2$v
X.pca <- xm %*% svd.v
xout <-  X.pca[,1]
return(xout)
})

# Each of the items in the 'll' list is a vector of 218 eigenvalues (one for each individual)

lld <- do.call(rbind.data.frame,ll)
lld2 <- t(lld)
rownames(lld2) <- cluster_genexp[[1]][,1]
colnames(lld2) <- names(cluster_genexp)

genemap <- read.csv("C:/Users/bre227/Dropbox/eQTL/Results/cluster_data/gene_clusters_hsa_map.csv")
keymeta <- read.table("C:/Users/bre227/Documents/R/AD_project/Working/key_metadata.txt", header = T)

lld3 <- lld2[as.character(keymeta$AIBL.Id),]
kp <- ifelse(is.na(keymeta$PET)==T,FALSE,TRUE) # create vector of TRUE/FALSE for recorded PET status
keymeta1 <- keymeta[kp,] 

 Y1 <- factor(as.character(keymeta1$PET))
 lld4 <- lld3[kp,] # remove individuals that do not have PET status recorded
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

fulldata1 <- fulldata[, which(unique(colnames(fulldata)) %in% genemap$root_gene_probe_id)] # take only those genes that are in genemap
length(intersect(as.character(genemap$root_gene_probe_id),colnames(fulldata)))

#fulldata2 <- fulldata1[kp, ] # remove individuals with no PET status recorded
#x2 <- t(fulldata2) # transpose data into matrix for limma fit
#fit.lm2 <- lmFit(x2, design_treat1)
#fit2.c <- contrasts.fit(fit.lm2, contrast.matrix_treat1)
#fit.eb2 <- eBayes(fit2.c)
### comparison treated - untreated
#topALL2 <- toptable(fit = fit2.c, eb = fit.eb2, adjust = "fdr", n = 244, genelist = colnames(fulldata1))

pbs <- genemap$root_gene_probe_id[genemap$ID %in% topALL1$ID]
length(unique(pbs)) # shows there are 102 unique gene probes here
topALL.2 <- results_PETposvneg[results_PETposvneg$probe_id %in% pbs, ]
topALL.2 <- dplyr::select(topALL.2, P.Value, probe_id)

# combine data frames
res.NetW$root_gene_probe_id <- as.character(res.NetW$root_gene_probe_id)
res.all <- dplyr::left_join(res.NetW, topALL.2, by = c("root_gene_probe_id" = "probe_id"))
res.all$svd_better <- res.all$P.Value.x < res.all$P.Value.y # create new column telling us whether the SVD p-value is lower than the initial limma p-value
length(which(res.all$svd_better == T))
# only for 21/244 networks does the SVD predict PET status better than the root gene alone

res_svd_better <- dplyr::filter(res.all, svd_better == T) %>% 
  dplyr::select(ID, P.Value.x, P.Value.y, name, root_gene_symbol, root_gene_probe_id)
res_svd_better$P.Value.diff <- res_svd_better$P.Value.x - res_svd_better$P.Value.y




