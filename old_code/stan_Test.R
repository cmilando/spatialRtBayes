library(rstan)
library(parallel)
library(invgamma)

num_cores <- parallel::detectCores()
options(mc.cores = num_cores - 4)
rstan_options(auto_write = TRUE)

#### INPUT
start_day <- 12 ## day for start estimating
window_length <- 40 ## smoothing window for estimaiton
si_shape <- 2 ## shape parameter for serical interval assuming gamma distribution
si_rate <- 0.5 ## rate parameter for serical interval assuming gamma distribution
si_t_max <- 14 ## maximum number of days with non-zero probability for serial interval

# Observed cases, matrix Y
incidence_data <- read.csv("https://raw.githubusercontent.com/zwzhou-biostat/hrt/main/data/example_data.csv")
head(incidence_data)

Y <- as.matrix(incidence_data)
for(i in 1:nrow(Y)) {
  for(j in 1:ncol(Y)) {
    Y[i,j] <- as.integer(Y[i,j])
  }
}
N <- as.integer(nrow(Y))
J <- as.integer(ncol(Y))
head(Y)

# Transfer matrix, matrix P (random example)
P <- matrix(c(0.8, 0.2, 0.1, 
              0.15, 0.6, 0.3,
              0.05, 0.2, 0.6), 3)

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
  beta_mu = 0.0,
  beta_sd = 1.0,
  sigma_shape = 2.0, ## vary this
  sigma_scale = 1.0,
  Z = as.integer(5)
)
# 
# m_hier <- stan(file="R/chad_stan.stan", 
#              data = stan_data, 
#              control = list(max_treedepth = 20,
#                             adapt_delta = 0.99))
# 
# saveRDS(m_hier, "stan_test_v2.RDS")
# 
# ## hmm sigma is not converging
# ## 
# 
# m_hier <- stan(file="R/stan_sliding.stan", 
#                data = stan_data, 
#                control = list(max_treedepth = 20,
#                               adapt_delta = 0.99))

#######
# re-running with tests
# The log of the inverse gamma density of y given shape alpha and scale beta
# shape = 2
hist(rinvgamma(100, shape = 1, scale = 1), breaks = c(0:50, 120))
stan_data$sigma_shape = 1
stan_data$Z = 0
stan_data
m_hier1 <- stan(file="R/stan_sliding.stan", 
               data = stan_data, 
               control = list(max_treedepth = 20,
                              adapt_delta = 0.99))
saveRDS(m_hier1, "stan_test_v3_shape1_z0.RDS")

#######
hist(rinvgamma(100, shape = 10, scale = 1))
stan_data$sigma_shape = 10
stan_data$Z = 0
stan_data
m_hier2 <- stan(file="R/stan_sliding.stan", 
               data = stan_data, 
               control = list(max_treedepth = 20,
                              adapt_delta = 0.99))

saveRDS(m_hier2, "stan_test_v3_shape10_z0.RDS")



# shape = 5
hist(rinvgamma(100, shape = 5, scale = 1))
stan_data$sigma_shape = 5
m_hier <- stan(file="R/stan_sliding.stan", 
               data = stan_data, 
               control = list(max_treedepth = 20,
                              adapt_delta = 0.99))
# launch_shinystan(m_hier)


m_hier <- readRDS("stan_test_v3_shape5.RDS")
m_hier@stanmodel
out <- rstan::extract(m_hier)

chain1 <- out$xsigma[1:1000, ]
chain2 <- out$xsigma[1001:2000, ]
chain3 <- out$xsigma[2001:3000, ]
chain4 <- out$xsigma[3001:4000, ]

hist(chain1[, 1], xlim = c(0.05, 0.15))
hist(chain2[, 1], xlim = c(0.05, 0.15))
hist(chain3[, 1], xlim = c(0.05, 0.15))
hist(chain4[, 1], xlim = c(0.05, 0.15))

dim(out$logR) # iter x day x state

chain1 <- out$logR[1:1000, 1, 1]
chain2 <- out$logR[1001:2000, 1, 1]
chain3 <- out$logR[2001:3000, 1, 1]
chain4 <- out$logR[3001:4000, 1, 1]

hist(chain1, xlim = c(2.5, 8.5))
hist(chain2, xlim = c(2.5, 8.5), add = T)
hist(chain3, xlim = c(2.5, 8.5), add = T)
hist(chain4, xlim = c(2.5, 8.5), add = T)

# library(shinystan)
# launch_shinystan(m_hier)

data_l <- lapply(1:3, function(i) {
  data.frame(
    x      = 1:nrow(incidence_data),
    y_real = incidence_data[,i],
    #R_real = Rlist(1:nrow(incidence_data), i),
    #R_real = R_rev[, i],
    y      = colMeans(out$M[, , i]),
    yl     = colQuantiles(out$M[, , i],probs=c(0.025)),
    yh     = colQuantiles(out$M[, , i],probs=c(0.975)),
    Rt     = colMeans(out$R[, ,i]),
    Rtl    = colQuantiles(out$R[, ,i],probs=c(0.025)),
    Rth    = colQuantiles(out$R[, ,i],probs=c(0.975)),
    region = i
  )
})

data_all <- do.call(rbind, data_l)

data_grouped <- data_all %>% group_by(., x, region)

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

###
library(lemon)
ggplot(data_all_summarise) +
  coord_cartesian(expand = F) +
  geom_ribbon(aes(x = x,ymin=yl,ymax=yh, fill = region), alpha=0.3) + 
  geom_line(aes(y = y_mean, x = x, color = region),linewidth = 0.5) +
  geom_line(aes(y = y_real, x = x, group = region), color = 'black', linewidth = 0.25) +
  theme_classic() +  xlab("Days") + ylab('Cases') + 
  facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
  theme(strip.background = element_blank())

#######
library(matrixStats)
library(tidyverse)

library(shinystan)
launch_shinystan(m_hier)

m_hier <- readRDS("stan_test_v1.RDS")

out <- rstan::extract(m_hier)

Ra <- function(t) (20*cos(t/500) + (0.8*t - 50)^2 - (0.115 * t)^3)/1000 + 0.8 
Rb <- function(t) (30*sin(t/150) + cos(t/20) - (t/50)^2)/8 - 0.006*t 
Rc <- function(t) (30*cos(t/150) + 2*sin(t/20) + 2*(t/50)^2)/20 - 0.005*t

plot(x = 1:200, y = Ra(1:200), type = 'l', col = 'red')
lines(x = 1:200, y = Rb(1:200), type = 'l', col = 'green')
lines(x = 1:200, y = Rc(1:200), type = 'l', col = 'blue')

Rlist = function(x, i) {
  if(i == 1) return(Ra(x))
  if(i == 2) return(Rb(x))
  if(i == 3) return(Rc(x))
}

data_l <- lapply(1:3, function(i) {
  data.frame(
    x      = 1:nrow(incidence_data),
    y_real = incidence_data[,i],
    #R_real = Rlist(1:nrow(incidence_data), i),
    #R_real = R_rev[, i],
    y      = colMeans(out$M[, , i]),
    yl     = colQuantiles(out$M[, , i],probs=c(0.025)),
    yh     = colQuantiles(out$M[, , i],probs=c(0.975)),
    Rt     = colMeans(out$R[, ,i]),
    Rtl    = colQuantiles(out$R[, ,i],probs=c(0.025)),
    Rth    = colQuantiles(out$R[, ,i],probs=c(0.975)),
    region = i
  )
})

data_all <- do.call(rbind, data_l)

data_grouped <- data_all %>% group_by(., x, region)

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

###
library(lemon)
ggplot(data_all_summarise) +
  coord_cartesian(expand = F) +
  geom_ribbon(aes(x = x,ymin=yl,ymax=yh, fill = region), alpha=0.3) + 
  geom_line(aes(y = y_mean, x = x, color = region),linewidth = 0.5) +
  geom_line(aes(y = y_real, x = x, group = region), color = 'black', linewidth = 0.25) +
  theme_classic() +  xlab("Days") + ylab('Cases') + 
  facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
  theme(strip.background = element_blank())

# adjust for some burn-in of R
ggplot(data_all_summarise %>% filter(x > 12), aes(x = x, color = region)) +
  #coord_cartesian(ylim = c(0, 5)) +
  geom_ribbon(aes(x = x,ymin=Rtl,ymax=Rth, fill = region), alpha=0.3) + 
  geom_line(aes(y=Rt_mean,x=x, color = region),linewidth = 0.5)+ 
  geom_line(aes(y=Rt_real,x=x, linetype = region), color = 'black', linewidth = 0.5)+
  geom_hline(yintercept = 1,color = "black", linewidth = 0.25,alpha=0.5)  +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number') +
  facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
  theme(strip.background = element_blank())


