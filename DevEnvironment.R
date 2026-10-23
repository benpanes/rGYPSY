pack.list <- c("devtools",
               "roxygen2",
               "testthat")

lapply(pack.list,
       library,
       character.only=T)

#Rcpp::compileAttributes() # Run if C++ files are changed / before pushing
setwd("H:/Shared drives/FGrOW/Projects/FRA/GYPSY R Version/R/rGYPSY")
getwd()


## Need to add zero filling
## Need to convert from species group level rows to columns for each species group's parameters
## Order will always be the same, zero filling with missing groups/values then skipping 0's
## Need to add Y2BH and Y2SH as output options from SI conversion function

library(data.table)
library(openxlsx)
df <- fread("H:/Shared drives/FGrOW/Projects/FRA/GYPSY R Version/R/rGYPSY/scratch/data/for_gypsy_test_bp.csv")
df2 <- read.xlsx("H:/Shared drives/FGrOW/Projects/FRA/GYPSY R Version/R/rGYPSY/scratch/data/test_gypsy_bp.xlsx")
