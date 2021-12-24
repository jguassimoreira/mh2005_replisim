######################################################################
## function to augment simulated datasets per M&H2005's description ##

#intputs: n_l1 - level 1 sample size
#         N_l2 - level 2 sample size
#         xdat - fixed level 1 (within cluster) data
#         zdat - fixed level 2 (between cluster) data

#outputs: dat - Augmented design matrix; i.e., data with repeated observations for larger sample sizes


augment_dat_legacy = function(n_l1, N_l2, xdat, zdat){
  
  #extend the design matrix depending on the level 1 and level 2 sample sizes
  
  if (n_l1 == 30) {
    
    dat = data.frame(subID = rep(xdat$subID, each = 6),
                     x = rep(xdat$x, each = 6))  #repeat the level 1 observations 6 times to get L1 n = 30
    
  } else if (n_l1 == 50) {
    
    dat = data.frame(subID = rep(xdat$subID, each = 10),
                     x = rep(xdat$x, each = 10)) #repeat the level 1 observations to get the L1 n = 50
    
    #M&H didn't describe how they 'repeated' their data in great detail, so this is our best guess
    #one consequence of this approach is that, at this point, there are essentially 6 or 10 duplicates
    #of each original level 1 variable per level 2 unit 
    
  } else {
    
    dat = data.frame(subID = xdat$subID, x = xdat$x)
  }
  
  dat = dplyr::inner_join(dat, zdat, by = 'subID') #join the newly replicated xdata (now in the 'dat' dataframe) w level 2 zdata 
  
  
  if (N_l2 == 50) {
    
    #repeat the level 2 observations
    #do so by repeating rows as many times as there are level 1 observations
    #since we start with 30 L2 units (the smallest L2 N in M&H), we need to add 20 more
    dat = rbind(dat, data.frame(subID = rep(paste0('sub',31:50), each = n_l1), x = dat$x[1:(n_l1*20)], z = dat$z[1:(n_l1*20)])) #adding 20
    
  } else if (N_l2 == 100) {
    
    #do the same thing here to get to 100, but this time we need to add 70 L2 units
    #It's unclear what M&H did in the original; I'm assuming they repeated the original data until they got to 100
    dat = rbind(dat,
                data.frame(subID = rep(paste0('sub',31:60), each = n_l1), x = dat$x[1:(n_l1*30)], z = dat$z[1:(n_l1*30)]), #adding 30 (60)
                data.frame(subID = rep(paste0('sub',61:90), each = n_l1), x = dat$x[1:(n_l1*30)], z = dat$z[1:(n_l1*30)]), #adding 30 (90)
                data.frame(subID = rep(paste0('sub',91:100), each = n_l1), x = dat$x[1:(n_l1*10)], z = dat$z[1:(n_l1*10)])) #adding 10 (100)
    
  }
  
  #The issue with setting up the data like this is that I think the ICC for X is either 0(/close to zero)
  #OR the data have a non-zero ICC but the predictor values are implausibly similar for different cases 
  #(which maybe results in an ICC that's too high?)
  
  return(dat)
  
}
