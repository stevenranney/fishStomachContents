
#Walleye standard weight (Ws) equation

calc_wae_wr <- function(weight, length){
  
  (weight/(10^(-5.453+3.180*(log10(length)))))*100

  }
