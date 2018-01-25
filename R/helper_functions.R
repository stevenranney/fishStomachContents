
#Calculate Walleye relative weight with standard weight (Ws) equation

calc_wae_wr <- function(weight, length){
  
  (weight/(10^(-5.453+3.180*(log10(length)))))*100

}

# Helper for assinging length classes

round_down <- function(x,to=10)
{
  to*(x %/% to)
}

# Calculate R^2 
R2 <- function(x, y, model){

  #Defines the total sum of squares for both models (totalSS=(n-1)*var(y))
  totalSS <- (length(x)-1)*var(y)
  #residSS=deviance(nls model)
  residSS <- deviance(model)

  return(1-(residSS/totalSS))

}
