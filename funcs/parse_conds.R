#######################################################################
## quick helper function to parse string with simulation conditions  ##

#intputs: path - path to rds file that holds simulation reps for one cell in the study

#outputs: The ICC level (numeric) used in the given condition and the other conditions (N and n)

parse_conds = function(path) {
  
  #split string to the rds file, get file name (last element in the flattened list)
  rds_file_name = unlist(strsplit(path, "/"))[length(unlist(strsplit(path, "/")))]
  #get the icc label (first object in the flattened vector)
  icc_lab = unlist(strsplit(rds_file_name, "_"))[1]
  #convert icc label to icc value by using case_when
  icc_val = dplyr::case_when(icc_lab == "icc1" ~ 0.1,
                             icc_lab == "icc2" ~ 0.2,
                             icc_lab == "icc3" ~ 0.3)
  #extract other labels
  n_lab = unlist(strsplit(rds_file_name, "_"))[2]
  N_lab = unlist(strsplit(rds_file_name, "_"))[3]
  
  return(list(icc_val, n_lab, N_lab))
  
}