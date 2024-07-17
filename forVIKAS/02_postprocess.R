library(rstan)
library(parallel)
library(invgamma)
library(tidyverse)
library(matrixStats)
library(lemon)
library(shinystan)

m_hier <- readRDS("stan_out.RDS")

out <- rstan::extract(m_hier)
# QC
launch_shinystan(m_hier)

# -------------------------------------------------------------
# get actual data
dim(out$M)


## CHANCGE THIS TO BE BY CHAIN
data_l <- lapply(1:J, function(i) {
  
  # 
  n_chains = N_CHAINS
  n_iters_per_chain = N_ITER
  
  sub_l <- lapply(1:n_chains, function(nc_i) {
    
    ii = 1 + n_iters_per_chain * (nc_i - 1)
    ix = n_iters_per_chain * (nc_i)
    
    data.frame(
      ##
      x      = 1:nrow(Y),
      y_real = Y[,i],
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

data_all %>%
  filter(region ==2, x == 2)

#
ggplot(data_all %>% filter(x > 5)) +
  coord_cartesian(expand = F) +
  geom_line(aes(x = x, y = Rt,
                color = paste0(chain))) +
  geom_line(aes(x = x, y = R_real,
                color = paste0(chain))) +
  facet_grid(region~., scales = 'free_y')

# ggplot(data_all) +
#   coord_cartesian(expand = F) +
#   geom_line(aes(x = x, y = Rt,
#                 color = paste0(chain))) +
#   facet_wrap(region~., scales = 'free_y')
# 
# 
# ggplot(data_all) +
#   coord_cartesian(expand = F) +
#   geom_ribbon(aes(x = x, ymin = yl, ymax = yh, 
#                   fill = factor(region),
#                   group = paste0(region, ".", chain)), 
#               color = 'black', alpha=0.3) +
#   facet_grid(region~chain)
# 
# ggplot(data_all %>% filter(x > 5)) +
#   coord_cartesian(expand = F) +
#   geom_ribbon(aes(x = x, ymin = Rtl, ymax = Rth, 
#                   fill = paste0(region, ".", chain),
#                   group = paste0(region, ".", chain)), 
#               color = 'black', alpha=0.13) +
#   facet_grid(region~.)



# geom_line(aes(y = y_mean, x = x, color = region),linewidth = 0.5) +
# geom_line(aes(y = y_real, x = x, group = region), color = 'black', linewidth = 0.25) +
# theme_classic() +  xlab("Days") + ylab('Cases') + 
# facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
# theme(strip.background = element_blank())

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
            Rth = mean(Rth))

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
N_ZEROS = 5
ggplot(data_all_summarise %>% filter(x > N_ZEROS), aes(x = x, color = region)) +
  coord_cartesian(ylim = c(0, 5)) +
  geom_ribbon(aes(x = x,ymin=Rtl,ymax=Rth, fill = region), alpha=0.3) + 
  geom_line(aes(y=Rt_mean,x=x, color = region),linewidth = 0.5)+ 
  geom_line(aes(y=Rt_real,x=x, linetype = region), color = 'black', linewidth = 0.5)+
  geom_hline(yintercept = 1,color = "black", linewidth = 0.25,alpha=0.5)  +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number') +
  facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
  theme(strip.background = element_blank())