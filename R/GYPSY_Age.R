# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

GYPSY_Age <- function(spp,
                      si,
                      ht,
                      sitype = "B", # T for Total or B for Breast Height
                      ndx = 1, # 1 for T, 2 for S, 3 for B
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
                             sitype = "B", # Input SI switch, use "T" for total
                             ndx = 4, # 4 for alternate SI, 5 for Y2BH, 6 for Y2SH
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
