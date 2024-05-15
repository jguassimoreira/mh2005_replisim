###########################################################################
## function to solve for upper and lower bounds of imbalanced L1 samples ##

#inputs: n_l1 - ideal level 1 sample size given no imbalance, or the average level 1 sample size given imbalance
#         imbalance
#         mod_obj - mlm power model object, used as part of data generation
#         imbalance - ratio of size between largest and smallest groups. Ideal/average n_l1 size falls in the middle. 

#outputs: Upper and lower bounds of level 1 sample sizes, to be used to draw L1 sample sizes from unif distribution

get_range_imbalance = function(n_l1, imbalance) { 
  
  # system of equation used to solve for the range follows
  # 2*n_l1 = lower + upper
  # imbalance = upper / lower
  # Rearranging using arithmetic...
  # upper = imbalance * lower
  # doing some substitution...
  
  lower = (2 * n_l1) / (1 + imbalance)
  
  upper = lower * imbalance
  
  return(list(lower = lower, upper = upper))
  
}
