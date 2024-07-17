library(rstan)
library(parallel)
library(invgamma)
library(tidyverse)
library(matrixStats)
library(lemon)
library(shinystan)

num_cores <- parallel::detectCores()
options(mc.cores = min(4, num_cores))
rstan_options(auto_write = TRUE)

# -------------------------------------------------------------
#### INPUT
window_length <- 1 ## smoothing window for estimaiton
si_shape <- 2       ## shape parameter for serial interval assuming gamma distribution
si_rate  <- 0.5     ## rate parameter for serial interval assuming gamma distribution
si_t_max <- 14      ## maximum number of days with non-zero probability for serial interval
beta_mu  <- 0.      ## prior mean for Rj(t) ~ N
beta_sd  <- 1.      ## prior sd   for Rj(t) ~ N
sigma_shape <- 2.   ## prior for sigma_j assuming inverse gamma
sigma_scale <- 1.   ## prior for sigma_j assuming inverse gamma

# -------------------------------------------------------------
# Observed cases, matrix Y
incidence_data <- read.csv("https://raw.githubusercontent.com/zwzhou-biostat/hrt/main/data/example_data.csv")
head(incidence_data)

## what about if I add a bunch of zeros ahead of it
N_zeros <- 20
leading_zeros <- data.frame(a = rep(0, N_zeros), 
                            b = rep(0, N_zeros),
                            c = rep(0, N_zeros))

Y <- rbind(leading_zeros, incidence_data)
head(Y)

## convert to integers
Y <- as.matrix(Y)
for(i in 1:nrow(Y)) {
  for(j in 1:ncol(Y)) {
    Y[i,j] <- as.integer(Y[i,j])
  }
}

# get N and J
N <- as.integer(nrow(Y))
J <- as.integer(ncol(Y))
head(Y)

# Transfer matrix, matrix P (random example)
P <- matrix(c(0.8, 0.2, 0.1, 
              0.15, 0.6, 0.3,
              0.05, 0.2, 0.6), 3)

# make daily
P_list <- lapply(1:N, function(x) P)
P <- do.call(rbind, P_list)

# Serial interval, vector W
S <- as.integer(si_t_max)

## obtain discretized gamma distribution for serial interval
W <- sapply(1:S, function(x){
  pgamma(x, si_shape, si_rate) - pgamma(x-1, si_shape, si_rate)
})   

W <- W/sum(W)
plot(W, type = 'l')

# Data list for Stan
stan_data <- list(
  N = N,
  J = J,
  Y = Y,
  P = P,
  S = S,
  W = W,
  beta_mu = beta_mu,
  beta_sd = beta_sd,
  sigma_shape = sigma_shape, 
  sigma_scale = sigma_scale,
  Z = as.integer(window_length) 
)

# -------------------------------------------------------------
# run the STAN code, takes ~ 10min
m_hier <- stan(file="stan_sliding_v2.stan",
             data = stan_data)

# extract output
out <- rstan::extract(m_hier)

# QC
launch_shinystan(m_hier)

# -------------------------------------------------------------
# plot timeseries

# get actual data
Ra <- function(t) (20*cos(t/500) + (0.8*t - 50)^2 - (0.115 * t)^3)/1000 + 0.8 
Rb <- function(t) (30*sin(t/150) + cos(t/20) - (t/50)^2)/8 - 0.006*t 
Rc <- function(t) (30*cos(t/150) + 2*sin(t/20) + 2*(t/50)^2)/20 - 0.005*t

ggplot() + theme_bw() +
  coord_fixed(ratio = 50) +
  geom_line(aes(x = 1:200, y = Ra(1:200)), color = 'red') +
  geom_line(aes(x = 1:200, y = Rb(1:200)), color = 'green') +
  geom_line(aes(x = 1:200, y = Rc(1:200)), color = 'blue')

Rx <- data.frame(a = Ra(1:(N - N_zeros)), 
                 b = Rb(1:(N - N_zeros)), 
                 c = Rc(1:(N - N_zeros)))
Rx <- rbind(leading_zeros, Rx)

## CHANCGE THIS TO BE BY CHAIN
data_l <- lapply(1:3, function(i) {
  
  # 
  n_chains = 4
  n_iters_per_chain = 1000
  
  sub_l <- lapply(1:n_chains, function(nc_i) {
    data.frame(
      x      = 1:nrow(Y),
      y_real = Y[,i],
      R_real = Rx[, i],
      y      = colMeans(out$M[, , i]),
      yl     = colQuantiles(out$M[, , i],probs=c(0.025)),
      yh     = colQuantiles(out$M[, , i],probs=c(0.975)),
      Rt     = colMeans(out$R[, ,i]),
      Rtl    = colQuantiles(out$R[, ,i],probs=c(0.025)),
      Rth    = colQuantiles(out$R[, ,i],probs=c(0.975)),
      region = i,
      chain = nc_i
    )
  })
  
  do.call(rbind, sub_l)
  
})

data_all <- do.call(rbind, data_l)

data_grouped <- data_all %>% group_by(., x, region, chain)

data_all_summarise <- data_grouped  %>% 
  mutate(region = paste0("region: ", region)) %>%
  summarise(.groups = 'keep',
            y_mean = mean(y),
            y_real = mean(y_real),
            Rt_real = mean(R_real),
            yl = mean(yl),
            yh = mean(yh),
            Rt_mean = mean(Rt),
            Rtl = mean(Rtl),
            Rth = mean(Rth),
            Rt_mean_test = mean(Rt>1))

data_all_summarise %>% arrange(region, x)

### plot exepected cases
ggplot(data_all_summarise) +
  coord_cartesian(expand = F) +
  geom_ribbon(aes(x = x,ymin=yl,ymax=yh, fill = region), alpha=0.3) + 
  geom_line(aes(y = y_mean, x = x, color = region),linewidth = 0.5) +
  geom_line(aes(y = y_real, x = x, group = region), color = 'black', linewidth = 0.25) +
  theme_classic() +  xlab("Days") + ylab('Cases') + 
  facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
  theme(strip.background = element_blank())

# plot R(t)
ggplot(data_all_summarise %>% filter(x > N_zeros), aes(x = x, color = region)) +
  coord_cartesian(ylim = c(0, 5)) +
  geom_ribbon(aes(x = x,ymin=Rtl,ymax=Rth, fill = region), alpha=0.3) + 
  geom_line(aes(y=Rt_mean,x=x, color = region),linewidth = 0.5)+ 
  geom_line(aes(y=Rt_real,x=x, linetype = region), color = 'black', linewidth = 0.5)+
  geom_hline(yintercept = 1,color = "black", linewidth = 0.25,alpha=0.5)  +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number') +
  facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
  theme(strip.background = element_blank())


