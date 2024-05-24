
require(EpiEstim)
require(dplyr)
require(ggplot2)
require(purrr)

setwd("YOUR_WORKING_DIRECTORY")
source(".\\adj_inc_EpiEstim.R")
incidence_data <- read.csv(".\\example_data.csv")


#### INPUT
start_day <- 12 ## day for start estimating
window_length <- 40 ## smoothing window for estimaiton
si_shape <- 2 ## shape parameter for serical interval assuming gamma distribution
si_rate <- 0.5 ## rate parameter for serical interval assuming gamma distribution
si_t_max <- 14 ## maximum number of days with non-zero probability for serial interval

## mobility matrix
P <- matrix(c(0.8, 0.2, 0.1,  ## mobility matrix
              0.15, 0.6, 0.3,
              0.05, 0.2, 0.6), 3)

#### OUTPUT
## the direct output from estimate_R_adj_inc keeps all outputs the same as EpiEsim, 
## but in a list format, with each element in the list correspond to one region

## obtain result from main function
result <- estimate_R_adj_inc(incidence_data, P, si_shape, si_rate, si_t_max, 
                             start_day, window_length)

## extrat R(t) estimation from the result
dat_regions <- do.call("rbind", lapply(1:length(result), function(x){
  dat <- as.data.frame(result[[x]]$R[, c(3, 5, 11)])
  dat$region <- names(result)[x]
  dat$day <- (start_day):(length(result[[x]]$dates) - window_length)
  dat
}))
## The extracted result 'dat_regions' contains the posterior mean of 
## R(t) as well as 2.5% quantile and 97.5% quantile of the posterior distribution

#### Generate plot for R(t)
dat_plot_r <- dat_regions %>% rename(r = "Mean(R)", rl = "Quantile.0.025(R)", 
                                     ru = "Quantile.0.975(R)")

g_r <- ggplot(dat_plot_r, aes(x = day, color = region)) +
  geom_line(mapping=aes(y=r,x=day, color = region),size=1.2)+lims(y = c(0, 3.5)) +
  geom_ribbon(mapping=aes(x=day,ymin=rl,ymax=ru, fill = region),alpha=0.4) +
  geom_hline(yintercept = 1,color = "black", size=1,alpha=0.5)  +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number')# 
  # + scale_y_continuous(trans = "log10")

g_r

#### Generate plot for estimated incidence

## obtain discretized gamma distribution for serial interval
w <- sapply(1:si_t_max, function(x){
  pgamma(x, si_shape, si_rate) - pgamma(x-1, si_shape, si_rate)
})  

w <- w/sum(w)

## recompute incidence with estimated R(t)
dat_inc_est_list <- lapply(1:3, function(m){
  dat_regions_R_est <- spread(dat_regions[, c(names(dat_regions)[m], "region", "day")], 
                              region, names(dat_regions)[m])[, -1]
  inc_est <- as.data.frame(t(sapply(2:nrow(dat_regions_R_est), function(t){
    rev(w[1:min(length(w), t)]) %*% 
      as.matrix(incidence_data[c(max(t-length(w) + 1, 1): t), 1:3,drop = F]) %*% 
      diag(dat_regions_R_est[t, ]) %*% P
  })))
  names(inc_est) <- names(incidence_data)
  inc_est$day <- (start_day + 1):(length(result[[1]]$dates) - window_length)
  inc_est_long <- gather(inc_est, key = "region", value = "y_est", -day)
  names(inc_est_long)[names(inc_est_long) == "y_est"] <- c("y", "yl", "yu")[m]
  inc_est_long
})

## plot for estimated incidence
dat_plot_inc <- dat_inc_est_list %>% reduce(left_join, by = c("day", "region"))

g_inc <- ggplot(dat_plot_inc, aes(x = day, color = region)) +
  geom_line(mapping=aes(y=y,x=day, color = region),size=1.2)+
  geom_ribbon(mapping=aes(x=day,ymin=yl,ymax=yu, fill = region),alpha=0.4) +
  theme_classic() +  xlab("Days") + ylab('Estimated Incidence')

g_inc