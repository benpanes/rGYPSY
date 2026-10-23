# GYPSY_SI (Calculate Site Index from GYPSY)
# Calculate site index in m from GYPSY model
# ndx(total, bh, years to BH, years to SH)
# 1,2,3,4

GYPSY_SI_vec <- function(spp,topht,age,ndx=2,iter=1000,tol=0.00000001){
  # Input error checks
  errors <- character()
  if(length(spp)>1){
    # Vector inputs
    bad_topht <- which(topht<0.1 | is.na(topht))
    if(length(bad_topht)>0) {
      errors <- c(errors,
                  paste("Top height must be greater than 0.1; check row(s):",
                        paste(bad_topht, collapse = ", ")))
    }
    bad_spp <- which(!spp%in%c("AW","PL","SW","SB"))
    if (length(bad_spp)>0) {
      errors <- c(errors,
                  paste("Species must be AW, PL, SW, or SB; check row(s):",
                        paste(bad_spp, collapse = ", ")))
    }
    # Age
    bad_age_awpl <- which(spp%in%c("AW","PL") & age<=7.5)
    if (length(bad_age_awpl)>0) {
      errors <- c(errors,
                  paste("Age must be >= 7.5 for AW or PL; check row(s):",
                        paste(bad_age_awpl, collapse = ", ")))
    }
    bad_age_swsb <- which(spp%in%c("SW","SB") & age<=9.5)
    if (length(bad_age_swsb)>0) {
      errors <- c(errors,
                  paste("Age must be >= 9.5 for SW or SB; check row(s):",
                        paste(bad_age_swsb, collapse = ", ")))
    }

  }
  if(length(spp)==1){
    # Single value inputs
    if(topht<0.1 |  is.na(topht)){
      errors <- c(errors,"Top height must be greater than 0.1")
    }
    if(!spp%in%c("AW","PL","SW","SB")){
      errors <- c(errors,"Species must be AW, PL, SW, or SB")
    }
    if(spp%in%c("AW","PL") & age<=7.5){
      errors <- c(errors,"Age must be >= 7.5 for AW or PL")
    }
    if(spp%in%c("SW","SB") & age<=9.5){
      errors <- c(errors,"Age must be >=9.5 for SW or SB")
    }
  }
  if(!ndx%in%c(1,2,3)){
    errors <- c(errors,"Invalid ndx. Must be 1, 2, or 3")
  }
  # Final
  if(length(errors)>0) {
    stop(paste(errors, collapse = "\n"))
  }

  SI_model <- function(spp,age,topht,iter,tol){
    AW_SI <- function(age,topht,si0){
      b1 = 9.908888
      b2 = -3.92451
      b3 = -0.32778
      b4 = 0.134376

      k1 = exp(b1+b2*sqrt(log(50+1))+b3*(log(si0))^2+b4*sqrt(50))
      k2 = si0^(b3*log(si0))
      k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
      Y2BH = exp((log(k3)/b2)^2)-1
      x10 = (1+exp(b1+b2*sqrt(log(0+age+1.0))+b3*log(si0)^2+b4*sqrt(50)))
      x20 = (1+exp(b1+b2*sqrt(log(50+1.0))+b3*log(si0)^2+b4*sqrt(50)))
      si1 = topht*x10/x20
      si0 = (si0+si1)/2
      si_bh = si1*(1+exp(b1+b2*sqrt(log(50+1))+b3*(log(si1))^2+b4*sqrt(50)))/
        (1+exp(b1+b2*sqrt(log(50+Y2BH+1))+b3*(log(si1))^2+b4*sqrt(50)))
      return(list(si0,si1,si_bh,Y2BH))
    }
    SB_SI <- function(age,topht,si0){
      b1 = 14.56236
      b2 = -6.04705
      b3 = -1.53715
      b4 = 0.240174

      k1 = exp(b1+b2*sqrt(log(50+1))+b3*(log(si0))^1+b4*sqrt(50))
      k2 = si0^b3
      k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
      Y2BH = exp((log(k3)/b2)^2)-1
      x10 = (1+exp(b1+b2*sqrt(log(0+age+1))+b3*log(si0)^1+b4*sqrt(50)))
      x20 = (1+exp(b1+b2*sqrt(log(50+1))+b3*log(si0)^1+b4*sqrt(50)))
      si1 = topht*x10/x20
      si0 = (si0+si1)/2
      si_bh = si1*(1+exp(b1+b2*sqrt(log(50+1))+b3*(log(si1))^1+b4*sqrt(50)))/
        (1+exp(b1+b2*sqrt(log(50+Y2BH+1))+b3*(log(si1))^1+b4*sqrt(50)))
      return(list(si0,si1,si_bh,Y2BH))
    }
    SW_SI <- function(age,topht,si0){
      b1 = 12.14943
      b2 = -3.77051
      b3 = -0.28534
      b4 = 0.165483

      k1 = exp(b1+b2*sqrt(log(50^2+1))+b3*log(si0)^2+b4*sqrt(50))
      k2 = si0^(b3*log(si0))
      k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
      Y2BH = sqrt(exp((log(k3)/b2)^2)-1)
      x10 = (1+exp(b1+b2*sqrt(log(0+age^2+1))+b3*log(si0)^2+b4*sqrt(50)))
      x20 = (1+exp(b1+b2*sqrt(log(50^2+1))+b3*log(si0)^2+b4*sqrt(50)))
      si1 = topht*x10/x20
      si0 = (si0+si1)/2
      si_bh = si1*(1+exp(b1+b2*sqrt(log(50^2+1))+b3*(log(si1))^2+b4*sqrt(50)))/
        (1+exp(b1+b2*sqrt(log((50+Y2BH)^2+1))+b3*(log(si1))^2+b4*sqrt(50)))
      return(list(si0,si1,si_bh,Y2BH))
    }
    PL_SI <- function(age,topht,si0){
      b1 = 12.84571
      b2 = -5.73936
      b3 = -0.91312
      b4 = 0.150668

      k1 = exp(b1+b2*sqrt(log(50+1))+b3*(log(si0))^1+b4*sqrt(50))
      k2 = si0^b3
      k3 = (si0*(1+k1)/1.3-1)/(exp(b1)*exp(b4*sqrt(50))*k2)
      Y2BH = exp((log(k3)/b2)^2)-1
      x10 = (1+exp(b1+b2*sqrt(log(0+age+1))+b3*log(si0)^1+b4*sqrt(50)))
      x20 = (1+exp(b1+b2*sqrt(log(50+1))+b3*log(si0)^1+b4*sqrt(50)))
      si1 = topht*x10/x20
      si0 = (si0+si1)/2
      si_bh = si1*(1+exp(b1+b2*sqrt(log(50+1))+b3*(log(si1))^1+b4*sqrt(50)))/
        (1+exp(b1+b2*sqrt(log((50+Y2BH)+1))+b3*(log(si1))^1+b4*sqrt(50)))
      return(list(si0,si1,si_bh,Y2BH))
    }

    df <- data.table::data.table("spp"=spp,
                                 "age"=age,
                                 "topht"=topht,
                                 "si0"=10,
                                 "si1"=100,
                                 "diff"=90)
    for(i in 1:iter){ # AW
      idx <- df[, spp=="AW" & diff > tol]
      if(!any(idx)) break
      df[idx,
         c("si0","si1","si_bh","Y2BH"):=
           AW_SI(age,topht,si0)]
      df[idx,
         diff:=abs(si0-si1)]
    }
    for(i in 1:iter){ # SB
      idx <- df[, spp=="SB" & diff > tol]
      if(!any(idx)) break
      df[idx,
         c("si0","si1","si_bh","Y2BH"):=
           SB_SI(age,topht,si0)]
      df[idx,
         diff:=abs(si0-si1)]
    }
    for(i in 1:iter){ # SW
      idx <- df[, spp=="SW" & diff > tol]
      if(!any(idx)) break
      df[idx,
         c("si0","si1","si_bh","Y2BH"):=
           SW_SI(age,topht,si0)]
      df[idx,
         diff:=abs(si0-si1)]
    }
    for(i in 1:iter){ # PL
      idx <- df[, spp=="PL" & diff > tol]
      if(!any(idx)) break
      df[idx,
         c("si0","si1","si_bh","Y2BH"):=
           PL_SI(age,topht,si0)]
      df[idx,
         diff:=abs(si0-si1)]
    }
    return(df[,.(si0,si1,si_bh,Y2BH)])
  }
  return(SI_model(spp, age, topht,iter,tol)[[ndx+1]])
}
