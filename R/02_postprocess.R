library(rstan)
library(parallel)
library(invgamma)
library(tidyverse)
library(matrixStats)
library(lemon)
library(shinystan)
source("00_generate_data.R")

m_hier <- readRDS("stan_out_weekly.RDS")

rstan::check_divergences(m_hier)
rstan::check_hmc_diagnostics(m_hier)

out <- rstan::extract(m_hier)
# QC
# launch_shinystan(m_hier)

# -------------------------------------------------------------
# get actual data
dim(out$M)
dim(out$xsigma)
N_CHAINS = 4
N_ITER = 4000

## CHANCGE THIS TO BE BY CHAIN
data_l <- lapply(1:J, function(i) {
  
  # 
  n_chains = N_CHAINS
  n_iters_per_chain = N_ITER
  
  sub_l <- lapply(1:N_CHAINS, function(nc_i) {
    
    ii = 1 + n_iters_per_chain * (nc_i - 1)
    ix = n_iters_per_chain * (nc_i)
    
    data.frame(
      ##
      x      = 1:nrow(R_this),
      y_real = N[,i],
      R_back = R_this[, i], ## R_this is the backcalcualted, Rmatrix is true
      R_real = Rmatrix[, i],
      ##
      y      = colMeans(out$M[ii:ix, , i]),
      yl     = colQuantiles(out$M[ii:ix, , i],probs=c(0.025)),
      yh     = colQuantiles(out$M[ii:ix, , i],probs=c(0.975)),
      ##
      Rt     = colMeans(out$R[ii:ix, ,i]),
      Rtl    = colQuantiles(out$R[ii:ix, ,i],probs=c(0.025)),
      Rth    = colQuantiles(out$R[ii:ix, ,i],probs=c(0.975)),
      ##
      region = i,
      chain = nc_i
    )
  })
  
  do.call(rbind, sub_l)
  
})

data_all <- do.call(rbind, data_l)
head(data_all)

#data_grouped <- data_all %>% group_by(., x, region, chain)

data_all_summarise <- data_all  %>% 
  mutate(region = paste0("region: ", region)) %>%
  group_by(x, region) %>%
  summarise(.groups = 'keep',
            y_mean = mean(y),
            y_real = mean(y_real),
            Rt_back = mean(R_back),
            Rt_real = mean(R_real),
            yl = mean(yl),
            yh = mean(yh),
            Rt_mean = mean(Rt),
            Rtl = mean(Rtl),
            Rth = mean(Rth))

data_all_summarise %>% arrange(region, x)

### plot exepected cases
p1 <- ggplot(data_all_summarise) +
  geom_ribbon(aes(x = x,ymin=yl,ymax=yh, fill = region), alpha=0.3) + 
  geom_line(aes(y = y_mean, x = x, color = region), linewidth = 0.5) +
  geom_line(aes(y = y_real, x = x, group = region), color = 'black', linewidth = 0.25) +
  theme_classic() +  xlab("Days") + ylab('Cases') + 
  facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
  theme(strip.background = element_blank())

# plot R(t)
p2 <- ggplot(data_all_summarise, aes(x = x, color = region)) +
  geom_ribbon(aes(x = x,ymin=Rtl,ymax=Rth, fill = region), alpha=0.3) + 
  geom_line(aes(y=Rt_mean,x=x, color = region),linewidth = 0.5)+ 
  geom_line(aes(y=Rt_real,x=x), color = 'black', linewidth = 0.5, linetype = '41')+
  geom_line(aes(y=Rt_back,x=x), color = 'black', linewidth = 0.5, linetype = '11')+
  geom_hline(yintercept = 1,color = "black", linewidth = 0.25,alpha=0.5)  +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number') +
  facet_rep_wrap(~region, nrow = 3) + 
  coord_cartesian(ylim = c(0, 5)) +
  theme(strip.background = element_blank())

library(patchwork)
p1 + theme(legend.position = 'none') + p2 + patchwork::plot_layout(nrow = 1)

#############################################################################
library(splines)

spline_df = vector("list", 3)
for(this_i in 1:3) {
y <- R_this[, this_i]

x <- 1:length(y)
N_KNOTS = 4 ## inversely proportional to the length of y
k <- quantile(x, probs = seq(from = 0.01, to = 0.99, length.out  = N_KNOTS))

b1 <- lm(y ~ ns(x, knots = k))
x.pred <- predict(b1, se = T)

head(x.pred$se.fit)
spline_df[[this_i]] <- data.frame(
  region = paste0('region: ', this_i),
  x = x,
  y = x.pred$fit, 
                        ylb = x.pred$fit - 1.96 * x.pred$se.fit,
                        yub = x.pred$fit + 1.96 * x.pred$se.fit)
}
spline_df <- do.call(rbind, spline_df)

ggplot(data_all_summarise, 
       aes(x = x, color = region)) +
  geom_ribbon(aes(x = x,ymin=Rtl,ymax=Rth, fill = region), alpha=0.3) + 
  geom_line(aes(y=Rt_mean,x=x, color = region),linewidth = 0.5)+ 
  geom_line(data = spline_df, aes(x = x, y = y, color= region),
            linewidth = 2.5) +
  #geom_line(aes(y=Rt_real,x=x), color = 'black', linewidth = 0.5, linetype = '41')+
  geom_line(aes(y=Rt_back,x=x), color = 'black', linewidth = 0.5, linetype = '11')+
  geom_hline(yintercept = 1,color = "black", linewidth = 0.25,alpha=0.5)  +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number') +
  facet_rep_wrap(~region, nrow = 3) + 
  coord_cartesian(ylim = c(0, 5)) +
  theme(strip.background = element_blank())

