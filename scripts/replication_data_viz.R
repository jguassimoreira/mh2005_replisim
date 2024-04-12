###Let's plot the results of our literature search/review

#first load the ggplot2 package which we'll use to plot the results 
library(ggplot2)
library(ggpubr)

#########
## Import helper functions, set path to data
#########
studyPath = file.path("/Users", "jguassimoreira", "Documents", "mh2005_replisim") #declare the path to the study
scriptPath = file.path(sprintf("%s", studyPath), "funcs") #path to helper funcs
file.sources = Sys.glob(file.path(sprintf("%s", scriptPath), "*.R")) #Construct paths to helper funcs, save in vector
invisible(sapply(file.sources,FUN=source)) #use sapply to source each file path in the helper func vector
datPath = file.path(sprintf("%s", studyPath), "data") #path where we'll eventually save our data

#########
## Manually input values from Maas & Hox Table 4
#########
ngrps = rep(c(30, 50, 100), each = 9)
grpsize = rep(rep(c(5, 30, 50), each = 3), 3)
icc = rep(c(0.1, 0.2, 0.3), 9)
u0  = c(.103, .086, .080, .088, .079, .095, .084, .082, .107,
        .065, .070, .080, .083, .080, .082, .062, .069, .073,
        .060, .060, .061, .053, .059, .056, .058, .071, .058)
u1  = c(.116, .093, .093, .064, .075, .089, .082, .091, .086,
        .072, .066, .084, .060, .080, .065, .072, .076, .071,
        .065, .059, .057, .050, .058, .056, .060, .059, .053)
e   = c(.069, .067, .051, .047, .053, .064, .070, .050, .053,
        .061, .058, .068, .051, .065, .052, .048, .052, .052,
        .060, .062, .050, .050, .043, .035, .046, .043, .048)
int = c(.064, .073, .072, .060, .064, .057, .065, .075, .043,
        .053, .059, .061, .065, .061, .053, .055, .052, .054,
        .068, .048, .061, .047, .053, .047, .054, .047, .048)
x  =  c(.060, .063, .067, .065, .049, .055, .062, .060, .062,
        .054, .057, .057, .056, .043, .070, .064, .050, .059,
        .047, .041, .048, .051, .051, .048, .041, .058, .061)
z  =  c(.056, .052, .071, .052, .041, .047, .037, .048, .060,
        .057, .051, .040, .058, .052, .040, .044, .061, .055,
        .051, .053, .065, .045, .039, .060, .037, .047, .055)
xz =  c(.058, .062, .065, .046, .051, .057, .048, .063, .053,
        .072, .047, .048, .048, .059, .054, .049, .056, .039,
        .054, .062, .043, .047, .046, .053, .041, .053, .052)

mh_t4 = data.frame(ngrps, grpsize, icc, u0, u1, e, int, x, z, xz)

#########
## Read in sim results, extract to make our own version of MH table 4
#########
conds_to_load = Sys.glob(file.path(sprintf("%s/%s", datPath, "replication"), "*.RDS")) #get file paths for all M&H simulation conditions
conds_list = lapply(conds_to_load, FUN = readRDS) #read in the conditions
params_to_analyze = c("int_var", "x_var", "sigma2", "intercept", "b_x", "b_z", "b_xz") #analyze the parameters
#extract data using the interaction analysis function
dataList = lapply(1:length(params_to_analyze), function(x) replication_interaction_analysis(conds_list, params_to_analyze[x], conds_to_load))

#put into dataframe
rep_dat = data.frame(ngrps = aggregate(non_coverage ~ icc + grpsize + ngrps, FUN = "mean", dataList[[1]]$data)$ngrps, 
                     grpsize =  aggregate(non_coverage ~ icc + grpsize + ngrps, FUN = "mean", dataList[[1]]$data)$grpsize,
                     icc =  aggregate(non_coverage ~ icc + grpsize + ngrps, FUN = "mean", dataList[[1]]$data)$icc,
                     u0 = aggregate(non_coverage ~ icc + grpsize + ngrps, FUN = "mean", dataList[[1]]$data)$non_coverage,
                     u1 = aggregate(non_coverage ~ icc + grpsize + ngrps, FUN = "mean", dataList[[2]]$data)$non_coverage,
                     e = aggregate(non_coverage ~ icc + grpsize + ngrps, FUN = "mean", dataList[[3]]$data)$non_coverage,
                     int = aggregate(non_coverage ~ icc + grpsize + ngrps, FUN = "mean", dataList[[4]]$data)$non_coverage,
                     x = aggregate(non_coverage ~ icc + grpsize + ngrps, FUN = "mean", dataList[[5]]$data)$non_coverage,
                     z = aggregate(non_coverage ~ icc + grpsize + ngrps, FUN = "mean", dataList[[6]]$data)$non_coverage,
                     xz = aggregate(non_coverage ~ icc + grpsize + ngrps, FUN = "mean", dataList[[7]]$data)$non_coverage)

#re-format the sim conditionst to be numeric
rep_dat$ngrps = dplyr::case_when(rep_dat$ngrps == "ngrps30" ~ 30,
                                 rep_dat$ngrps == "ngrps50" ~ 50,
                                 rep_dat$ngrps == "ngrps100" ~ 100)

rep_dat$grpsize = dplyr::case_when(rep_dat$grpsize == "grpsize5" ~ 5,
                                   rep_dat$grpsize == "grpsize30" ~ 30,
                                   rep_dat$grpsize == "grpsize50" ~ 50)

rep_dat$icc = dplyr::case_when(rep_dat$icc == "0.1" ~ 0.1,
                               rep_dat$icc == "0.2" ~ 0.2,
                               rep_dat$icc == "0.3" ~ 0.3)


#########
## Start plotting
#########

### u0
u0_mh = ggplot(mh_t4, aes(x = ngrps, y = u0, color = factor(icc))) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize)), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")

u0_rep = ggplot(rep_dat, aes(x = ngrps, y = u0, color = factor(icc))) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize)), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")


### u1
u1_mh = ggplot(mh_t4, aes(x = ngrps, y = u1, color = factor(icc))) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize),), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")

u1_rep = ggplot(rep_dat, aes(x = ngrps, y = u1, color = factor(icc)), show.legend = T) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize)), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")


### e
e_mh = ggplot(mh_t4, aes(x = ngrps, y = e, color = factor(icc))) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize),), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")

e_rep = ggplot(rep_dat, aes(x = ngrps, y = e, color = factor(icc)), show.legend = T) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize)), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")


### int
int_mh = ggplot(mh_t4, aes(x = ngrps, y = int, color = factor(icc))) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize),), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")

int_rep = ggplot(rep_dat, aes(x = ngrps, y = int, color = factor(icc)), show.legend = T) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize)), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")


### x
x_mh = ggplot(mh_t4, aes(x = ngrps, y = x, color = factor(icc))) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize),), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")

x_rep = ggplot(rep_dat, aes(x = ngrps, y = x, color = factor(icc)), show.legend = T) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize)), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none") 


### z
z_mh = ggplot(mh_t4, aes(x = ngrps, y = z, color = factor(icc))) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize),), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")

z_rep = ggplot(rep_dat, aes(x = ngrps, y = z, color = factor(icc)), show.legend = T) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize)), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")


### xz
xz_mh = ggplot(mh_t4, aes(x = ngrps, y = xz, color = factor(icc))) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize),), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none")

xz_rep = ggplot(rep_dat, aes(x = ngrps, y = xz, color = factor(icc)), show.legend = T) +
  geom_line(aes(shape = factor(grpsize)), size = .85) +
  geom_hline(yintercept = .075, linetype = "dashed") +
  geom_point(aes(shape = factor(grpsize)), size = 2.5) +
  scale_color_grey() +  
  labs(x = "Level-2 Sample Size", y = "Non-coverage") +
  scale_y_continuous(limits = c(0, .1175)) +
  theme_test(base_size = 13.5) + theme(legend.position = "none") 


##put into a grid, save
ggsave("/Users/jguassimoreira/Documents/mh2005_replisim/plots/figure3.png",
       ggarrange(u0_mh, u0_rep,
                 u1_mh, u1_rep,
                 e_mh, e_rep,
                 nrow = 3, ncol = 2),
       width = 7, height = 10, dpi = 900)

ggsave("/Users/jguassimoreira/Documents/mh2005_replisim/plots/figure4.png",
       ggarrange(int_mh, int_rep,
                 x_mh, x_rep,
                 z_mh, z_rep,
                 xz_mh, xz_rep,
                 nrow = 4, ncol = 2),
       width = 7, height = 10, dpi = 900)


