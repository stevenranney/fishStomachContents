
#Calculate Walleye relative weight with standard weight (Ws) equation
calc_wae_wr <- function(weight, length){
  
  (weight/(10^(-5.453+3.180*(log10(length)))))*100

}

#Calculate smallmouth relative weight with standard weight (Ws) equation
calc_smb_wr <- function(weight, length){
  
  (weight/(10^(-5.329+3.2*(log10(length)))))*100
  
}

# Assign walleye length categories
assign_wae_psd <- function(data){
  
  ifelse((data>=250)&(data<380), "S-Q",
         ifelse((data>=380)&(data<510), "Q-P",
                ifelse((data>=510)&(data<630), "P-M",
                       ifelse((data>=630)&(data<760), "M-T",
                              ifelse(data>=760, ">T", "SS")))))
}

# Assign smb length categories
assign_smb_psd <- function(data){
  
  ifelse((data>=180)&(data<280), "S-Q",
                        ifelse((data>=280)&(data<350), "Q-P",
                               ifelse((data>=350)&(data<430), "P-M",
                                      ifelse((data>=430)&(data<510), "M-T",
                                             ifelse(data>=510, ">T", "SS")))))
}


# Assign length categories
assign_length_cat <- function(data, vec){
  
  ifelse((data>=vec[1])&(data<vec[2]), "S-Q",
         ifelse((data>=vec[2])&(data<vec[3]), "Q-P",
                ifelse((data>=vec[3])&(data<vec[4]), "P-M",
                       ifelse((data>=vec[4])&(data<vec[5]), "M-T",
                              ifelse(data>=vec[5], ">T", "SS")))))
  
  
}

# Helper for assinging length classes

round_down <- function(x,to=10)
{
  to*(x %/% to)
}

#Simple fx to apply a linear model
apply_lm <- function(x, model){
  
  model$coefficients[[2]]*x+model$coefficients[[1]]
  
}

# Calculate R^2 
R2 <- function(x, y, model){

  #Defines the total sum of squares for both models (totalSS=(n-1)*var(y))
  totalSS <- (length(x)-1)*var(y)
  #residSS=deviance(nls model)
  residSS <- deviance(model)

  return(1-(residSS/totalSS))

}

# Calculate R^1 for a fitted quantile regression model; requires the full model and an 
# intercept-only model
R1 <- function(q_mod, q_mod_null){

  1 - (q_mod$rho/q_mod_null$rho)

}

# Calculate percent difference between two numbers
calc_perc_diff <- function(x, y){
  
  ((x-y)/y)*100
}