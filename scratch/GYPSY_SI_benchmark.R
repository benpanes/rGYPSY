# =============================================================================
# GYPSY_SI — Rcpp wrapper and benchmark
# =============================================================================
# Requirements:
#   install.packages(c("Rcpp", "data.table", "microbenchmark"))
#   A working C++ compiler:
#     Windows -> Rtools (https://cran.r-project.org/bin/windows/Rtools/)
#     macOS   -> Xcode Command Line Tools  (`xcode-select --install`)
#     Linux   -> gcc/g++ (usually pre-installed)
# =============================================================================

#setwd("scratch")
library(Rcpp)
library(data.table)
library(microbenchmark)

# --- 1. Compile the C++ source ---
#   Place GYPSY_SI_cpp.cpp in your working directory, then:
sourceCpp("src/GYPSY_SI_cpp.cpp")
sourceCpp("src/GYPSY_Age_cpp.cpp")
#   This compiles once and makes GYPSY_SI_cpp() available in your R session.
#   Recompilation is only needed when the .cpp file changes.

# --- 2. R wrapper (preserves the original interface exactly) ----------------
#   Drop-in replacement for GYPSY_SI_vec(). Passes all arguments through to
#   the compiled C++ function; validation is handled in C++.

GYPSY_SI <- function(spp,
                     ht,
                     age,
                     agetype = "T",
                     ndx = 2, # 1 for total, 2  for B, 3 for Y2BH, 4 for Y2SH
                     iter = 1000,
                     tol = 1e-8,
                     bh = 1.3,
                     sh = 0.3){
  GYPSY_SI_cpp(
    spp     = as.character(spp),
    ht      = as.numeric(ht),
    age     = as.numeric(age),
    agetype = as.character(agetype),
    ndx     = as.integer(ndx),
    iter    = as.integer(iter),
    tol     = as.numeric(tol),
    bh      = as.numeric(bh),
    sh      = as.numeric(sh)
  )
}

GYPSY_Age <- function(spp,
                      si,
                      ht,
                      sitype = "B",
                      ndx = 1, # 1 for T, 2 for B, 3 for S
                      iter = 1000,
                      tol = 1e-8,
                      bh = 1.3,
                      sh = 0.3){
  GYPSY_Age_cpp(
    spp     = as.character(spp),
    si      = as.numeric(si),
    h      = as.numeric(ht),
    sitype  = as.character(sitype),
    ndx     = as.integer(ndx),
    iter    = as.integer(iter),
    tol     = as.numeric(tol),
    bh      = as.numeric(bh),
    sh      = as.numeric(sh)
  )
}

GYPSY_SI_convert <- function(spp,
                             si,
                             sitype = "B", # Switch, use "T" for total SI to BH
                             ndx = 4,
                             iter = 1000,
                             tol = 1e-8,
                             bh = 1.3,
                             sh = 0.3){
  GYPSY_Age_cpp(
    spp     = as.character(spp),
    si      = as.numeric(si),
    sitype  = as.character(sitype),
    ndx     = as.integer(ndx),
    iter    = as.integer(iter),
    tol     = as.numeric(tol),
    bh      = as.numeric(bh),
    sh      = as.numeric(sh)
  )
}

# --- 3. Verify outputs match the R version ----------------------------------
#   Run this block once after compiling to confirm numerical equivalence.
#   Tolerance of 1e-6 is well within the convergence tolerance of both versions.

#source("scratch/GYPSY_SI_vec.R")   # load original R function

# --- Test data
trees <- fread("scratch/data/trees.csv")
trees_meas <- fread("scratch/data/trees_measurement.csv")
trees_merge <- trees_meas[trees,
                          on=c("company",
                               "company_plot_number",
                               "tree_number")]
# Uppercase species
trees_merge[,species:=toupper(species)]
# Assign species groups
trees_merge[,spp_grp:=fcase(
  species%in%c("AW","BW","PB","AX"),"AW",
  species%in%c("PL","PJ","PW","PF","PX","PA","LT","LA","LW"),"PL",
  species%in%c("SW","SE","SX","FB","FA","FD"),"SW",
  species=="SB","SB"
)]
# Convert breast height and stump age to total


#trees_merge[dbh_age>0 & total_age%in%c(0,NA),
#            total_age:=fcase(
#              spp_grp=="AW",dbh_age+4,
#              spp_grp=="PL",dbh_age+8,
#              spp_grp=="SB",dbh_age+15,
#              spp_grp=="SW",dbh_age+12)]
#trees_merge[stump_age>0 & total_age%in%c(0,NA),
#            total_age:=fcase(
#              spp_grp=="AW",round(stump_age+0.5),
#              spp_grp=="PL",stump_age+3,
#              spp_grp=="SB",stump_age+5,
#              spp_grp=="SW",stump_age+4
#            )]
# Filter to eligible trees
age_trees <- trees_merge[height>0.1 &
                           ((spp_grp=="AW" &
                             (total_age>=7.5 |
                              stump_age>=8   |
                              dbh_age  >=11.5))|
                            (spp_grp=="PL" &
                             (total_age>=7.5 |
                              stump_age>=10.5|
                              dbh_age  >=15.5))|
                            (spp_grp=="SB" &
                             (total_age>=9.5 |
                              stump_age>=14.5|
                              dbh_age  >=24.5))|
                            (spp_grp=="SW" &
                             (total_age>=9.5 |
                              stump_age>=13.5|
                              dbh_age  >=21.5)))]

test_spp   <- c("AW", "PL", "SW", "SB")
test_topht <- c(15.2, 12.8, 18.4, 10.1)
test_age   <- c(10, 25.0, 40.0, 35.0)
test_sibh  <- c(33.97015, 23.69902, 23.43431, 15.91066)
test_sit   <- c(33.79709, 22.42038, 21.78578, 14.11113)

GYPSY_SI_convert(test_spp, test_sibh) ### Test output ages next

#cpp_result <- age_trees[,GYPSY_SI(spp_grp,height,total_age,3)]


#test_spp   <- "SB"
#test_topht <- 1
#test_age   <- 60

#r_result <- GYPSY_SI_vec(test_spp, test_topht, test_age, n=2)
GYPSY_SI(test_spp, test_topht, test_age, "T", 1)


age_trees[!is.na(total_age),
          Y2BH:=GYPSY_SI(spp_grp,height,total_age,ndx=3)]
age_trees[!is.na(total_age),
          Y2SH:=GYPSY_SI(spp_grp,height,total_age,ndx=4)]

age_trees[is.na(total_age) & !is.na(stump_age),
          Y2BH:=GYPSY_SI(spp_grp,height,stump_age,"S",3)]
age_trees[is.na(total_age) & !is.na(stump_age),
          Y2SH:=GYPSY_SI(spp_grp,height,stump_age,"S",4)]

age_trees[is.na(total_age) & is.na(stump_age) & !is.na(dbh_age),
          Y2BH:=GYPSY_SI(spp_grp,height,dbh_age,"B",3)]
age_trees[is.na(total_age) & is.na(stump_age) & !is.na(dbh_age),
          Y2SH:=GYPSY_SI(spp_grp,height,dbh_age,"B",4)]

age_trees[,range(Y2BH),by="spp_grp"]

age_trees[,.(spp_grp,Y2SH=GYPSY_SI(spp_grp,height,total_age,"S",3))][,range(Y2SH,na.rm=T),by="spp_grp"]

system.time(
  r_result   <- age_trees[,GYPSY_SI_vec(spp_grp,height,total_age,3)]
)
system.time(
  cpp_result <- age_trees[,GYPSY_SI(spp_grp,height,total_age,3)]
)

cat("R result:   ", round(r_result,   8), "\n")
cat("C++ result: ", round(cpp_result, 8), "\n")
cat("Max diff:   ", max(abs(r_result - cpp_result)), "\n")
stopifnot(all(abs(r_result - cpp_result) < 1e-6))
cat("Outputs match.\n\n")

# --- 4. Benchmark -----------------------------------------------------------
#   Adjust n to reflect a realistic dataset size for your use case.

set.seed(42)
n <- 10000   # <-- change this to your typical row count

spp_pool <- c("AW", "PL", "SW", "SB")
age_min  <- c(AW = 8, PL = 8, SW = 10, SB = 10)

bm_spp   <- sample(spp_pool, n, replace = TRUE)
bm_topht <- runif(n, 5, 30)
bm_age   <- mapply(function(s) runif(1, age_min[s] + 1, 100), bm_spp)

results <- microbenchmark(
  R_version   = GYPSY_SI_vec(bm_spp, bm_topht, bm_age, ndx = 2),
  Cpp_version = GYPSY_SI(bm_spp, bm_topht, bm_age, ndx = 2),
  times = 20L
)

print(results, unit = "ms")
# Reported in milliseconds. Divide by 1000 for seconds.
# The "median" column is the most reliable comparison point.
