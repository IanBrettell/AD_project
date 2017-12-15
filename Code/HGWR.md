Happy git with R test
================
Ian Brettell
14 December 2017

``` r
exdata <- read.delim("~/R/AD_project/Data/Expression/AIBL_expression_set/AIBL_Gene_Expression.txt", sep = " ", header = T)
metadata <- read.delim("~/R/AD_project/Data/Expression/aibl-ids-6.0.0-201712010300.txt", sep = "\t", header = T)
ids <- read.delim("~/R/AD_project/Data/Expression/AIBL_Gene_Expression_IDs.txt", header = F)
```
