##########################################################################
## function to compute Wald confidence intervals of variance parameters ##

#intputs: model - lmer model object

#outputs: wci - Wald confindence intervals for elements of variance parameters (in SD/corr scale)

#adapted from: https://rpubs.com/bbolker/varwald

compute_wald_ci_varcomp = function(model){ 
  
  #if (isREML(model)) {
  #  warning("refitting model with ML") #refit model with maximum likelihood
  #  model = lme4::refitML(model)
  #}
  
  dd = lme4::devfun2(model,useSc=TRUE,signames=FALSE) #extract the deviance wrt the variance parameters
  vv = as.data.frame(lme4::VarCorr(model)) #extract the variance parameters themselves
  
  vnames = c(sprintf(".sig%02d",1:3),".sigma") #create vector of names for var parameters
  pars = setNames(vv[c(1,3,2,4),"sdcor"],vnames) #set the names of the variance parameters
  
  hh1 = numDeriv::hessian(dd,pars) #take the hessian (2nd deriv) of the deviance for the var params
  vv2 = 2*solve(hh1) #take inverse matrix to find the SEs
  
  dimnames(vv2) = list(vnames,vnames) #these next few lines calculate and extract the CIs for us
  vhack = list(coefficients=pars,vcov=vv2)
  wci = confint.default(vhack)
  
  return(wci)
  
  }
