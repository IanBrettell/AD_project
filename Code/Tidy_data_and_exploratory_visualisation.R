####################################################
# Author: Ian Brettell
# Date created: 171207
# Purpose: exploratory analysis of mRNA expression data
####################################################


####################################################
# IMPORT FILES
####################################################

exdata <- read.delim("AIBL_Gene_Expression.txt", sep = " ", header = T)
metadata <- read.delim("aibl-ids-6.0.0-201712010300.txt", sep = "\t", header = T)
ids <- read.delim("AIBL_Gene_Expression_IDs.txt", header = F)

####################################################
# LOAD PACKAGES
####################################################

library(tidyverse)
library(rafalib)
library(dummies)

####################################################
# TIDY AND COMBINE METADATA
####################################################

# remove 'X' from columns in 'exdata'

excols <- colnames(exdata)
excols <- gsub("X", "", excols)
colnames(exdata) <- excols

# sort IDs and rename for joining

ids <- data.frame(sort(ids))
colnames(ids) <- "AIBL.Id"

# missing_samples <- setdiff(1:287, ids)

# extract metadata for columns of interest and join to IDs

meta1 <- filter(metadata, Collection == "1") %>% 
  select(AIBL.Id,
         Demographic.Sex, 
         Image.PET.Amyloid.PIB_NAV.Status, 
         Demographic.ApoE.genotype, 
         Demographic.YearMonthOfBirth,
         Progress.Summary.Date.of.NP.assessment) %>% 
  right_join(ids, by = "AIBL.Id")

## create subset of data on which to test plots etc.

test.dat <- exdata[1:100, ] # extract first 100 rows
test.dat <- data.frame(t(test.dat)) # transpose
excols <- colnames(test.dat) # assign column names to vector
excols <- gsub("X", "", excols) # remove the X prefixes
colnames(test.dat) <- excols # replace the column names 

test.dat$AIBL.Id <- as.integer(rownames(test.dat))
test.dat <- left_join(meta1, test.dat, by = "AIBL.Id")

####################################################
# CREATE INITIAL EXPLORATORY PLOTS
####################################################

ggplot(data = test.dat) +
  geom_point(mapping = aes(x = AIBL.Id, y = test.dat$`2315554`, colour = Image.PET.Amyloid.PIB_NAV.Status))


####################################################
# CREATE PCA PLOTS
####################################################

# using the 'dummies' package, create a new data frame converting all variables into integers
test.dat1 <- dummy.data.frame(test.dat, names = c("Demographic.Sex", 
                                                  "Image.PET.Amyloid.PIB_NAV.Status",
                                                  "Demographic.ApoE.genotype",
                                                  "Progress.Summary.Date.of.NP.assessment"))

# divide the new data into test and train

pca.train <- test.dat1[1:nrow(test.dat1),]
pca.test <- test.dat1[-(1:nrow(test.dat1)),]

# principal component analysis

prin_comp <- prcomp(pca.train, scale. = T)


####################################################

test.dat2 <- exdata[1:100, ]
test.dat2 <- data.frame(t(test.dat2))

prin_comp <- prcomp(test.dat2, scale. = T)
names(prin_comp)

#outputs the mean of variables
prin_comp$center

#outputs the standard deviation of variables
prin_comp$scale

#outputs the principal component loading
prin_comp$rotation
prin_comp$rotation[1:5,1:4]

dim(prin_comp$x)

biplot(prin_comp, scale = 0)
std_dev <- prin_comp$sdev
pr_var <- std_dev^2
pr_var[1:10]

prop_varex <- pr_var/sum(pr_var)
prop_varex[1:20]

plot(prop_varex, 
     xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

plot(cumsum(prop_varex), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")
