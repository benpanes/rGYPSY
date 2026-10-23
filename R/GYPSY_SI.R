# GYPSY_SI (Calculate Site Index from GYPSY)
# Calculate site index in m from GYPSY model
# ndx(total, bh, years to BH, years to SH)
# 1,2,3,4

GYPSY_SI <- function(spp,
                     ht,
                     age,
                     agetype = "T", # T = Total, B = BH, S = SH
                     ndx = 2, # Output: 1 = SI_t, 2 = SI_BH, 3 = Y2BH, 4 = Y2SH
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

