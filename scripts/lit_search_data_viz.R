library(ggplot2)

icc_vals = c(0.62,0.48,0.52,0.55,0.36,0.33,0.47,0.17,0.25,0.44,0.12,
             0.24,0.27,0.32,0.38,0.22,0.56,0.49,0.4,0.01,0.13,0.22,0.23,
             0.39,0.74,0.8,0.27,0.5,0.58,0.66,0.62,0.34,0.11,0.29,0.65,
             0.51,0.59,0.57,0.06,0.04,0.09,0.17,0.01,0.21,0.63,0.14,0.19,
             0.16,0.57,0.2,0.16,0.14,0.05,0.07,0.08,0.13,0.07,0.34,0.06,
             0.21,0.18,0.32,0.12,0.17,0.03,0.83)

icc_df = data.frame(icc_vals)

ggplot(icc_df, aes(x=icc_vals)) + 
  geom_density(fill = 'lightgrey') + scale_color_grey() + theme_classic()

n_vals = c(3,15,7,3,66.7,19.5882352941176,2.33333333333333,49.6071428571429,
           8.15584415584416,14,2.93846153846154,96.65,62.43,16.9117647058824,
           7.91304347826087,14,6.6,31,129.714285714286,NA,3,131.393939393939,
           19.475,4,16.3333333333333,24.6944444444444,6,4.09473684210526,40.19,
           87.32,NA,8.31685393258427,3935.06349206349,NA,66.67,66.67,233.33,118,2,
           2,14.87,2,906.375,16,4.59863945578231,20,65,2.78889990089197,
           7.91891891891892,20,397.89,25.4259259259259,2)

n_df = data.frame(n_vals)

ggplot(n_df, aes(x=n_vals)) + 
  geom_density(fill = 'lightgrey') + scale_color_grey() + theme_classic()

N_vals = c(114,190,70,419,447,110,17,102,84,77,164,130,20,
           58,68,92,87,84,146,7,1150,114,33,40,102,3,72,78,
           95,213,1091,129,109,445,63,3,3,150,83,48,39,648,
           56,33,52,147,49968,101,1009,74,74,81,81,91,54,171)

N_df = data.frame(N_vals)

ggplot(N_df, aes(x=N_vals)) + 
  geom_density(fill = 'lightgrey') + scale_color_grey() + theme_classic()

l1_slopes = c(3,6,1,2,2,16,8,3,3,4,5,2,3,3,5,1,1,3,1,3,1,4,3,2,5,3,2,6,2,8,11,
              2,1,2,0,3,1,9,6,7,8,3,3,1,6,2,1)

l1_slopes_df = data.frame(l1_slopes)

ggplot(l1_slopes_df, aes(x=l1_slopes)) + 
  geom_density(fill = 'lightgrey') + scale_color_grey() + theme_classic()

l2_slopes = c(1,9,1,2,10,12,1,4,1,2,3,5,5,5,8,1,5,7,8,0,	
              2,12,4,3,4,1,3,9,6,1,2,0,1,1,2,5,4,0,2,1,3,
              4,3,1,6,4,6,1)

l2_slopes_df = data.frame(l2_slopes)

ggplot(l2_slopes_df, aes(x=l2_slopes)) +
  geom_density(fill = 'lightgrey') + scale_color_grey() + theme_classic()

cross_lev_int = c(3,0,0,0,7,1,0,0,3,4,6,0,1,0,0,0,0,2,3,0,0,
                  11,0,0,0,1,1,2,0,0,2,0,0,0,1,2,0,0,0,3,1)

cross_lev_int_df = data.frame(cross_lev_int)

ggplot(cross_lev_int_df, aes(x=cross_lev_int)) +
  geom_density(fill = "lightgrey") + scale_color_grey() + theme_classic()

n_rfx = c(1,2,2,1,4,3,4,1,1,1,3,3,1,1,2,3,1,2,1,1,1,2,2,1,1,2,1,2)
n_rfx_df = data.frame(n_rfx)

ggplot(n_rfx_df, aes(x=n_rfx)) +
  geom_density(fill = "lightgrey") + scale_color_grey() + theme_classic()