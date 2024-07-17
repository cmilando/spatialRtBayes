library(rstan)
library(parallel)

num_cores <- parallel::detectCores()
options(mc.cores = num_cores - 4)

stan_data <- list()

incidence_data <- read.csv("https://raw.githubusercontent.com/zwzhou-biostat/hrt/main/data/example_data.csv")
head(incidence_data)

w <- sapply(1:14, function(x){
  pgamma(x, 2, 0.5) - pgamma(x-1, 2, 0.5)
})  

w <- w/sum(w)
w <- c(w, rep(0, nrow(incidence_data) - length(w)))

stan_data$SI <- w

stan_data$casesa <- as.integer(incidence_data$a)
stan_data$casesb <- as.integer(incidence_data$b)
stan_data$casesc <- as.integer(incidence_data$c)

stan_data$N <- as.integer(nrow(incidence_data))

m_hier <- stan(file="test.stan", 
               data = stan_data)



saveRDS(m_hier, "stan_ZW_test_v1.RDS")

####
m_hier <- readRDS("stan_ZW_test_v1.RDS")
out <- rstan::extract(m_hier)

library(matrixStats)
library(tidyverse)

data_a<-data.frame(
  x=1:nrow(incidence_data),
  y_real=stan_data$casesa,
  y = colMeans(out$predictiona),
  yl = colQuantiles(out$predictiona,probs=c(0.025)),
  yh = colQuantiles(out$predictiona,probs=c(0.975)),
  Rt_real = R_rev[, 1],
  Rt = colMeans(out$Rta),
  Rtl = colQuantiles(out$Rta,probs=c(0.025)),
  Rth = colQuantiles(out$Rta,probs=c(0.975)),
  region = "a"
)

data_b<-data.frame(
  x=1:nrow(incidence_data),
  y_real=stan_data$casesb,
  y = colMeans(out$predictionb),
  yl = colQuantiles(out$predictionb,probs=c(0.025)),
  yh = colQuantiles(out$predictionb,probs=c(0.975)),
  Rt_real = R_rev[, 2],
  Rt = colMeans(out$Rtb),
  Rtl = colQuantiles(out$Rtb,probs=c(0.025)),
  Rth = colQuantiles(out$Rtb,probs=c(0.975)),
  region = "b"
)


data_c<-data.frame(
  x=1:nrow(incidence_data),
  y_real=stan_data$casesc,
  y = colMeans(out$predictionc),
  yl = colQuantiles(out$predictionc,probs=c(0.025)),
  yh = colQuantiles(out$predictionc,probs=c(0.975)),
  Rt_real = R_rev[, 3],
  Rt = colMeans(out$Rtc),
  Rtl = colQuantiles(out$Rtc,probs=c(0.025)),
  Rth = colQuantiles(out$Rtc,probs=c(0.975)),
  region = "c"
)


data_all <- rbind(data_a, data_b, data_c)

data_grouped <- data_all %>% group_by(., x, region)

data_all_summarise <- data_grouped  %>% summarise(., y_mean = mean(y),
                                                  y_real = mean(y_real),
                                                  yl = mean(yl),
                                                  yh = mean(yh),
                                                  Rt_real = mean(Rt_real),
                                                  Rt_mean = mean(Rt),
                                                  Rtl = mean(Rtl),
                                                  Rth = mean(Rth),
                                                  Rt_mean_test = mean(Rt>1))

data_all_summarise %>% arrange(region, x)

library(lemon)
ggplot(data_all_summarise) +
  # geom_smooth(data,mapping=aes(y=y,x=x, color = region),size=1.2) +
  geom_line(aes(y = y_mean, x = x, color = region),size=1.2) +
  geom_line(aes(y = y_real, x = x, group = region), color = 'black') +
  # geom_bar(data = data, aes(y = y_real, fill = region), stat='identity', alpha=0.5) +
  geom_ribbon(aes(x=x,ymin=yl,ymax=yh, fill = region),alpha=0.3) + 
  theme_classic() +  xlab("Days") + ylab('Infections') + 
  facet_wrap(~region)

# adjust for some burn-in of R
ggplot(data_all_summarise %>% filter(x > 12), aes(x = x, color = region)) +
  #coord_cartesian(ylim = c(0, 5)) +
  geom_line(aes(y=Rt_mean,x=x, color = region),linewidth = 0.5)+ 
  geom_line(aes(y=Rt_real,x=x, linetype = region), color = 'black', linewidth = 0.5)+
  geom_hline(yintercept = 1,color = "black", linewidth = 0.25,alpha=0.5)  +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number') +
  facet_rep_wrap(~region, nrow = 3, scales = 'free_y') + 
  theme(strip.background = element_blank())

