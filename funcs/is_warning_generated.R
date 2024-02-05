###################################################
## function to check for lmer convergence issues ##

#inputs: m - fitted model

#outputs: Logical value indicated whether model failed to converge (TRUE) or not (FALSE)

is_warning_generated = function(m) {
  df = summary(m) #save summary of model
  !is.null(df$optinfo$conv$lme4$messages) && #check for warning and 'failed to converge' text
    grepl('failed to converge', df$optinfo$conv$lme4$messages)
}