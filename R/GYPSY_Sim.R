# Project stand with GYPSY
# Variables:
#   SI, Age, TopHt, BA, N




GYPSY_PSI <- function(stand,
                      aw_ps, # First function that requires data widening
                      aw_age_t,
                      aw_si_t,
                      sb_ps,
                      sb_age_t,
                      sb_si_t,
                      sw_ps,
                      sw_age_t,
                      sw_si_t,
                      pl_ps,
                      pl_age_t,
                      pl_si_t){
  GYPSY_PSI_cpp(stand,
                aw_ps,
                aw_age_t,
                aw_si_t,
                sb_ps,
                sb_age_t,
                sb_si_t,
                sw_ps,
                sw_age_t,
                sw_si_t,
                pl_ps,
                pl_age_t,
                pl_si_t)
}

GYPSY_SDF <- function(stand,
                      aw_age_b,
                      aw_n,
                      aw_si_b,
                      aw_psi,
                      sb_age_t,
                      sb_n,
                      sb_si_b,
                      sb_psi,
                      sw_age_t,
                      sw_n,
                      sw_si_b,
                      sw_psi,
                      pl_age_t,
                      pl_n,
                      pl_si_b,
                      pl_psi){
  GYPSY_SDF_cpp(stand,
                aw_age_b,
                aw_n,
                aw_si_b,
                aw_psi,
                sb_age_t,
                sb_n,
                sb_si_b,
                sb_psi,
                sw_age_t,
                sw_n,
                sw_si_b,
                sw_psi,
                pl_age_t,
                pl_n,
                pl_si_b,
                pl_psi)
}

GYPSY_Ht <- function(spp,
                     si_t,
                     age_t){
  GYPSY_Ht_cpp(spp,
               si_t,
               age_t)
}

GYPSY_N <- function(){

}

GYPSY_BA <- function(){

}

GYPSY_Vol <- function(){

}

GYPSY_Sim <- function(stand,  # Outputs: stand id, stand age, total age, ba, ps, den, mden, topht, qmd, tvol, mvol
                      stand_age,
                      spp,
                      si      = NA_real_,
                      age     = NA_real_,
                      ht      = NA_real_,
                      ps      = NA_real_,
                      ba      = NA_real_,
                      sitype  = "B",   # T for Total or B for Breast Height
                      agetype = "T",   # "T", # T = Total, B = BH, S = SH
                      simlen  = 200,   # Simulation length
                      psflag  = 1,     # Percent stocking override, set to 0 turn off PS model
                      baflag  = 1,     # BA localization override, set to 0 turn off BA localization
                      bh      = 1.3,   # Optional breast height specification
                      sh      = 0.3,   # Optional stump height specification
                      mspp    = "ALL", # Species for merchantability criteria, either specify by species or use same standards for all species
                      sdob    = 15,    # Optional minimum stump diameter specification
                      tdib    = 10,    # Optional minimum top diameter specification
                      iter    = 1000,  # Maximum iterations for hill-climbing algorithms
                      tol     = 1e-8){ # Minimum difference between iterations before extracting results from hill-climbing algorithms

  input_length <- length(stand)
  inputs <- list(si  = si,
                 age = age,
                 ht  = ht,
                 ps  = ps,
                 ba  = ba)


  ### add stand age checks
  ### Most error checks should be part of the R code

  for(i in seq_along(inputs)){ # Replicate blank inputs (allows for partial blank columns)
    if(any(is.na(inputs[[i]])) & length(inputs[[i]])==1){
      inputs[[i]] <- rep(NA_real_, input_length)
    }
  }
  si  <- si
  age <- age
  ht  <- ht
  ps  <- ps
  ba  <- ba

  miss_input <- c() # Check for missing input data
  for(i in 1:input_length){
    if(sum(is.na(si[i]) + is.na(age[i]) + is.na(ht[i]))>1){
      miss_input <- append(miss_input, i)
    }
  }
  if(length(miss_input>0)){
    input_message <- paste(miss_input, collapse=", ")
    stop(paste0("Simulation requires 2 of SI, Age, and Top Height. Check row(s): ",input_message))
  }
  if(input_length!=length(stand_age) |
     input_length!=length(spp) |
     input_length!=length(si) |
     input_length!=length(age)|
     input_length!=length(ht) |
     input_length!=length(ps) |
     input_length!=length(ba)){
    stop("Inputs must be the same length")
  }
  dt1 <- data.table("stand"     = stand,     # Combine inputs into single table
                    "stand_age" = stand_age,
                    "spp"       = toupper(spp),
                    "si"        = si,
                    "age"       = age,
                    "ht"        = ht,
                    "ps"        = ps,
                    "ba"        = ba)
  i <- dt1[, .I[.N > 1], by = .(stand, spp)] # Find duplicate stand species groups
  if(length(i > 0)){
    duplicate_message <- paste(i, collapse=", ")
    stop(paste0("Duplicate stand species groups in row(s): ", duplicate_message))
  }
  i <- dt1[, .I[spp%notin%c("AW", )]]

  dt1[!is.na(ht) & !is.na(age), # Following established logic, recalculate SI if age and top height are entered.
      si_t:=GYPSY_SI(spp, ht, age, agetype = "T", ndx = 1, iter, tol, bh)]
  dt1[!si_t,
      si_bh:=GYPSY_SI_convert(spp, si_t, sitype = "T", ndx = 4, iter, tol, bh)] # Convert to SI_BH
  if(sitype==1){
    dt1[is.na(ht) | is.na(age), # SI must be supplied if age or height are blank
        si_t:=si]
    dt1[is.na(ht) | is.na(age),
        si_bh:=GYPSY_SI_convert(spp, si_t, sitype = "T", ndx = 4, iter, tol, bh)]

  } else if(sitype==2){
    dt1[is.na(ht) | is.na(age),
        si_t:=GYPSY_SI_convert(spp, si, sitype = "B", ndx = 4, iter, tol, bh)]
    dt1[is.na(ht) | is.na(age),
        si_bh:=si]
  }
  # Age calculations
  if(agetype=="T"){
    dt1[!is.na(age),
        age_t:=age] # Fill from input if provided
  } else if(agetype=="B"){
    dt1[!is.na(age),
        age_b:=age]
  }
  dt1[is.na(age_t), # Else calculate
      age_t:=GYPSY_Age(spp, si_t, ht, sitype, ndx = 1, iter, tol, bh)] # Starting from SI_total doesn't require convergence
  dt1[is.na(age_b),
      age_b:=GYPSY_Age(spp, si_t, ht, sitype, ndx = 2, iter, tol, bh)]
  dt1[,y2bh:=GYPSY_SI_convert(spp, si, sitype = "T", ndx = 5)]
  dt1[is.na(ht),
      ht:=GYPSY_Ht(spp, age_t, si_t)] ## This might not be used for anything, possibly leave out

  # Pivot species wider
  dt2 <- dt1[,.(stand_age=max(stand_age(na.rm=T)))]
  src_col <- c("age_t",
               "age_b",
               "y2bh",
               "si_t",
               "si_b",
               "ht",
               "n",
               "ps",
               "ba")
  add_col <- paste("aw_", src_col)
  dt2[dt1[spp=="AW"],
               on="stand",
               (add_col) := mget(src_col)]
  add_col <- paste("sb_", src_col)
  dt2[dt1[spp=="SB"],
               on="stand",
               (add_col) := mget(src_col)]
  add_col <- paste("sw_", src_col)
  dt2[dt1[spp=="SW"],
               on="stand",
               (add_col) := mget(src_col)]
  add_col <- paste("pl_", src_col)
  dt2[dt1[spp=="PL"],
      on="stand",
      (add_col) := mget(src_col)]

  # Zero-fill
  fil_cols <- grep("(_ps|_n|_ba)$", names(dt2), value = T)

  dt2[, (grep("(_ps|_n|_ba)$", names(dt2), value = T)) := lapply(.SD, function(x) { # Grab all ps, n, and ba columns ### Should this grab all columns? Might want to zero fill things like SI or age
    x[is.na(x)] <- 0 # Replace NA with 0
    x # Return x
  }), .SDcols = fil_cols]

  # Calculate other parameters (PSI, SDF, n0, sc0, bainc0)
  dt2[,c("aw_psi",
         "sb_psi",
         "sw_psi",
         "pl_psi") := GYPSY_PSI()]

  dt2[,c("aw_sdf",
         "sb_sdf",
         "sw_sdf",
         "pl_sdf") := GYPSY_SDF()]

  # Grow forward
  #--- SAS age calculation equivalents

  # Calculate values for each year





}
