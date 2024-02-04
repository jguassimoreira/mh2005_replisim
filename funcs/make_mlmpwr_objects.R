##################################################################
## function to make model objects that will be used to sim data ##

#inputs: mod - type of model (i.e., 1 L1 slope & cl int, 2 such slopes & int, or 3 such slopes & ints)
#        sim_icc - Y variable ICC
#        w_r2 - proportion of variance explained by the within variables
#        b_r2 - proportion of variance explained by the between variables
#        cl_r2 - proportion of variance explained by the cross level interactions (product terms)
#        raneff_cor - correlation between random intercept and slope (non-zero for only mod1)

#outputs: mlm_power object 


make_mlmpwr_objects = function(mod, sim_icc, w_r2, b_r2, cl_r2, raneff_cor) {

  if (mod == 1) {
    
    # NOTE: The 'weight' argument for some of the functions below refers to the contribution of a given variable
    # to the variance explained metric. All weights are set to one, ensuring that multi-variable models have
    # vars that contribute equally to the variance explained. 
    
    model = (
      effect_size(
        icc = sim_icc, #icc
        within = w_r2, #var explained by within vars
        between = b_r2, #var explained by between vars
        product = cl_r2, #var explained by cl ints (product)
        random_slope = 0.1 #var explained by random slopes - kept constant for the sim
      )
      + outcome('y', mean = 1, sd = .95) #mean and SD are hard coded, consistent with original M&H sim
      + within_predictor('x1', icc = 0, weight = 1) #specify one within var, icc will always be zero because we're assuming cwc 
      + between_predictor('z1', weight = 1) #label for only between var
      + product('x1','z1', weight = 1)  #product var
      + random_slope('x1', weight = 1)  #random slope
      + correlations(randeff = raneff_cor) #correlation between random effects
    )
    
  } else if (mod == 2) {
    
    model = (
      effect_size(
        icc = sim_icc, #icc
        within = w_r2, #var explained by within vars
        between = b_r2, #var explained by between vars
        product = cl_r2, #var explained by cl ints (product)
        random_slope = 0.1 #var explained by random slopes
      )
      + outcome('y', mean = 1, sd = .95) #mean and SD are hard coded, consistent with original M&H sim
      + within_predictor('x1', icc = 0, weight = 1) #specify first var, icc will always be zero because we're assuming cwc 
      + within_predictor('x2', icc = 0, weight = 1) #specify second var
      + between_predictor('z1', weight = 1) #label for only between var
      + product('x1','z1', weight = 1)  #first product var
      + product('x2','z1', weight = 1)  #second product var
      + random_slope('x1', weight = 1)  #first random slope
      + random_slope('x2', weight = 1)  #second random slope
      + correlations(randeff = 0) #correlation between random effects
    )
    
  } else if (mod == 3) {
    
    model = (
      effect_size(
        icc = sim_icc, #icc
        within = w_r2, #var explained by within vars
        between = b_r2, #var explained by between vars
        product = cl_r2, #var explained by cl ints (product)
        random_slope = 0.1 #var explained by random slopes
      )
      + outcome('y', mean = 1, sd = .95) #mean and SD are hard coded, consistent with original M&H sim
      + within_predictor('x1', icc = 0, weight = 1) #specify first var, icc will always be zero because we're assuming cwc 
      + within_predictor('x2', icc = 0, weight = 1) #specify second var
      + within_predictor('x3', icc = 0, weight = 1) #specific third var
      + between_predictor('z1', weight = 1) #label for only between var
      + product('x1','z1', weight = 1)  #first product var
      + product('x2','z1', weight = 1)  #second product var
      + product('x3','z1', weight = 1)  #third product var
      + random_slope('x1', weight = 1)  #first random slope
      + random_slope('x2', weight = 1)  #second random slope
      + random_slope('x3', weight = 1)  #third random slope
      + correlations(randeff = 0) #correlation between random effects
    )
    
  } else {
    
    return('invalid model configuration')
    
  }
  
  return(model) #return model
  
}

