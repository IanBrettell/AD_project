### Bring in Expression and IDS data ###

#################################################
### Create modified IDS for Roche import ###
#################################################
ids <- read.delim("~/Documents/Work/2017/Nov2017/Misc/aibl-ids-6-7/aibl-ids-6.0.0-201712010300.txt")
ids1 <- data.frame(ids[,c(1,2,grep("Amyloid",names(ids),fixed=F),
																				grep("Simple",names(ids),fixed=F), 
																				grep("Progress.Summary.Date.of.NP.assessment",names(ids),fixed=F), 
																				grep("ApoE.genotype", names(ids),fixed=F),
																				grep("Site", names(ids),fixed=F),
																				grep("Sex",names(ids),fixed=F),
																				grep("Amyloid",names(ids),fixed=F))])

names(ids1)[1] <- "AIBL.ID"

#################################################
### create age and apoe4 ###
#################################################
ids1$Age <- as.numeric((as.Date(ids1$Progress.Summary.Date.of.NP.assessment,format=c("%d/%m/%Y")))-(as.Date(paste("15",substr(ids1$Demographic.YearMonthOfBirth,start=5,stop=6),substr(ids1$Demographic.YearMonthOfBirth,start=1,stop=4),sep="/"),format = c("%d/%m/%Y"))))/365.25

ids1$apoe4 <- ifelse(ids1$Demographic.ApoE.genotype=="E2/E2"|ids1$Demographic.ApoE.genotype=="E3/E2"|ids1$Demographic.ApoE.genotype=="E3/E3",0,ifelse(ids1$Demographic.ApoE.genotype=="E4/E2"|ids1$Demographic.ApoE.genotype=="E4/E3"|ids1$Demographic.ApoE.genotype=="E4/E4",1,NA))

### Create binary PET status ###
ids1$PET <- ifelse(ids1$Image.PET.Amyloid.PIB_NAV.Status == "Positive" | ids1$Image.PET.Amyloid.Florbetapir.Status== "Positive" | ids1$Image.PET.Amyloid.Flutemetamol.Status == "Positive", "POS", ifelse(ids1$Image.PET.Amyloid.PIB_NAV.Status == "Negative" | ids1$Image.PET.Amyloid.Florbetapir.Status== "Negative" | ids1$Image.PET.Amyloid.Flutemetamol.Status == "Negative", "NEG",NA))

#################################################
### Create a combined SUVR Variable (BeCKeT) ###
#################################################
ids1$suvr.com <- ifelse(ids1$Image.PET.Amyloid.PIB_NAV.Status != "",ids1$Image.PET.Amyloid.PIB_NAV.SUVR,
									 ifelse(ids1$Image.PET.Amyloid.Florbetapir.Status!= "",ids1$Image.PET.Amyloid.Florbetapir.BeCKeT,
									 			ifelse(ids1$Image.PET.Amyloid.Flutemetamol.Status != "",ids1$Image.PET.Amyloid.Flutemetamol.BeCKeT ,NA)))
ids1$Tracer <- ifelse(ids1$Image.PET.Amyloid.PIB_NAV.Status!="","PiB",ifelse(ids1$Image.PET.Amyloid.Florbetapir.Status!="","Flor",ifelse(ids1$Image.PET.Amyloid.Flutemetamol.Status!="","Flute",NA)))


### Bring in expression data ###
exp <- read.delim("C:/Users/bre227/Dropbox/eQTL/Data/AIBL_expression_set/AIBL_Gene_Expression_UpdtdDec2017.txt",sep="")
e1 <- t(exp[1:100,])
id.numbers <- gsub("X","",rownames(e1))
		
ids2 <- ids1[ids1$AIBL.ID %in% id.numbers,]



