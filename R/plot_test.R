#######
m_hier <- readRDS("stan_test_v3_shape10.RDS"); out <- rstan::extract(m_hier)
m_hier <- readRDS("stan_test_v3_shape10_z0.RDS"); out <- rstan::extract(m_hier)
m_hier <- readRDS("stan_test_v3_shape1.RDS"); out <- rstan::extract(m_hier)
m_hier <- readRDS("stan_test_v3_shape1_z0.RDS"); out <- rstan::extract(m_hier)



chain1 <- out$xsigma[1:1000, ]
chain2 <- out$xsigma[1001:2000, ]
chain3 <- out$xsigma[2001:3000, ]
chain4 <- out$xsigma[3001:4000, ]

hist(chain1[, 1], xlim = c(0, 5))
hist(chain2[, 1], add = T)
hist(chain3[, 1], add = T)
hist(chain4[, 1], add = T)

hist(chain1[, 1], xlim = c(0.05, 0.15))
hist(chain2[, 1], xlim = c(0.05, 0.15), add = T)
hist(chain3[, 1], xlim = c(0.05, 0.15), add = T)
hist(chain4[, 1], xlim = c(0.05, 0.15), add = T)

dim(out$logR) # iter x day x state

chain1 <- out$logR[1:1000, 1, 1]
chain2 <- out$logR[1001:2000, 1, 1]
chain3 <- out$logR[2001:3000, 1, 1]
chain4 <- out$logR[3001:4000, 1, 1]

hist(chain1, xlim = c(-5, 10))
hist(chain2, xlim = c(0, 10), add = T)
hist(chain3, xlim = c(0, 10), add = T)
hist(chain4, xlim = c(0, 10), add = T, col ='red')

# library(shinystan)
# launch_shinystan(m_hier)

library(matrixStats)
library(tidyverse)

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
            #Rt_real = mean(R_real),
            yl = mean(yl),
            yh = mean(yh),
            Rt_mean = mean(Rt),
            Rtl = mean(Rtl),
            Rth = mean(Rth),
            Rt_mean_test = mean(Rt>1))

data_all_summarise %>% arrange(region, x)

###
# library(lemon)
# ggplot(data_all_summarise) +
#   coord_cartesian(expand = F) +
#   geom_ribbon(aes(x = x,ymin=yl,ymax=yh, fill = region), alpha=0.3) + 
#   geom_line(aes(y = y_mean, x = x, color = region),linewidth = 0.5) +
#   geom_line(aes(y = y_real, x = x, group = region), color = 'black', linewidth = 0.25) +
#   theme_classic() +  xlab("Days") + ylab('Cases') + 
#   facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
#   theme(strip.background = element_blank())

ggplot(data_all_summarise %>% filter(x >= 5), aes(x = x, color = region)) +
  coord_cartesian(ylim = c(0, 15)) +
  geom_ribbon(aes(x = x,ymin=Rtl,ymax=Rth, fill = region), alpha=0.3) + 
  geom_line(aes(y=Rt_mean,x=x, color = region),linewidth = 0.5)+ 
  #geom_line(aes(y=Rt_real,x=x, linetype = region), color = 'black', linewidth = 0.5)+
  geom_hline(yintercept = 1,color = "black", linewidth = 0.25,alpha=0.5)  +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number') +
  facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
  theme(strip.background = element_blank())
