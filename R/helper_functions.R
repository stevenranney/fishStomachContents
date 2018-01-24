
#Walleye standard weight (Ws) equation

calc_wae_wr <- function(weight, length){
  
  (weight/(10^(-5.453+3.180*(log10(length)))))*100

}

# Helper for assinging length classes

round_down <- function(x,to=10)
{
  to*(x %/% to)
}

# Calculate R^2 from a non-linear model