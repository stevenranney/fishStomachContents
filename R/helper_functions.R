
#Calculate Wr from weight, length, and vector of coefficients
calc_rel_weight <- function(weight, length, vec){

    (weight/(10^(vec[1]+vec[2]*(log10(length)))))*100

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
round_down <- function(x, to=10)
{
  to*(x %/% to)
}

#Simple fx to apply a linear model
apply_lm <- function(x, model){
  
  model$coefficients[[2]]*x+model$coefficients[[1]]
  
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

#Write two figures at once 
save_as_png_tiff <- function(filepath, plot = last_plot(), width = NA, height = NA, units = c("in", "cm", "mm")){
  
  ggsave(paste0(filepath, ".png"), plot = plot, width = width, height = height, units = "in")
  ggsave(paste0(filepath, ".tiff"), plot = plot, width = height, height = height, units = "in")

  }