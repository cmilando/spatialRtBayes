

args<-commandArgs(TRUE)
currind <-as.integer(args[1])


## simulate base on individual data
set.seed(2020 + currind)
library(rstan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(gridExtra)
library(parallel)
m <- 3 # number of regions
# pop <- rnorm(m, 50000, 10000) %>% round()
# initial_inf_rate <- rbeta(m, 1, 500)
# initial_inf_num <- round(pop*initial_inf_rate)

# initial_inf_num <- c(10, 10, 10)
# N_day <- 200


w <- sapply(1:14, function(x){
  pgamma(x, 2, 0.5) - pgamma(x-1, 2, 0.5)
})  

w <- w/sum(w)


P <- matrix(c(0.8, 0.2, 0.1, 
              0.15, 0.6, 0.3,
              0.05, 0.2, 0.6), 3)

N_simulated_list <- readRDS("/projectnb2/hrtgrp/zhenwei/simulation/v1/N_simulated.rds")

N_simulated_data <- N_simulated_list[[currind]]


stan_data<-list()
stan_data$N  <- nrow(N_simulated_data)
stan_data$casesa <- N_simulated_data$a
# stan_data$casesa <- rep(2, 182)
stan_data$casesb <- N_simulated_data$b
# stan_data$casesb <- rep(2, 182)
stan_data$casesc <- N_simulated_data$c
# stan_data$N2 <- 60
# 
# stan_data$logra <- log(R_est_sim_data[, 1])
# stan_data$logrb <- log(R_est_sim_data[, 2])
# stan_data$logrc <- log(R_est_sim_data[, 3])



w <- sapply(1:14, function(x){
  pgamma(x, 2, 0.5) - pgamma(x-1, 2, 0.5)
})  

w <- w/sum(w)

## ?? this is different in the NIMBLE CODE?
w <- c(w, rep(0, nrow(N_simulated_data) - length(w)))
###

stan_data$SI <- w

# saveRDS(N_simulated, "/rprojectnb2/zwzhou-thesis/rstan/simulated_data_exo.rds")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## no imported infections
stan_model1 <- stan_model(model_code = 
                            "// Stan model for simple linear regression

data {
  int<lower=1> N; // days of observed data

  int casesa[N]; // reported cases
  int casesb[N]; // reported cases
  int casesc[N]; // reported cases

  real SI[N]; // fixed SI using empirical data
  
}

transformed data {
  vector[N] SI_rev; // SI in reverse order
  for(i in 1:N)
  SI_rev[i] = SI[N-i+1];
  
}

parameters {

  vector[N] weekly_effecta;
  vector[N] weekly_effectb;
  vector[N] weekly_effectc;
  
  vector<lower=0, upper=1>[N] weekly_rhoa;
  vector<lower=0, upper=1>[N] weekly_rhob;
  vector<lower=0, upper=1>[N] weekly_rhoc;

  real<lower=0> weekly_sda;
  real<lower=0> weekly_sdb;
  real<lower=0> weekly_sdc;

}

transformed parameters {

  vector[N] predictiona=rep_vector(1e-5,N);
  vector[N] predictionb=rep_vector(1e-5,N);
  vector[N] predictionc=rep_vector(1e-5,N);
  
  vector<lower=0>[N] Rta;
  vector<lower=0>[N] Rtb;
  vector<lower=0>[N] Rtc;
  
  {
  
    Rta[1:N] = exp(weekly_effecta[1:N]);
    Rtb[1:N] = exp(weekly_effectb[1:N]);
    Rtc[1:N] = exp(weekly_effectc[1:N]);

    for (i in 2:N) {
    
      real convolutiona = dot_product(predictiona[1:(i-1)], tail(SI_rev, i-1));
      real convolutionb = dot_product(predictionb[1:(i-1)], tail(SI_rev, i-1));
      real convolutionc = dot_product(predictionc[1:(i-1)], tail(SI_rev, i-1));
      
      predictiona[i] = predictiona[i] + Rta[i] * convolutiona * 0.8 + Rtb[i] * convolutionb * 0.2 + Rtc[i] * convolutionc * 0.1;
      predictionb[i] = predictionb[i] + Rta[i] * convolutiona * 0.15 + Rtb[i] * convolutionb * 0.6 + Rtc[i] * convolutionc * 0.3;
      predictionc[i] = predictionc[i] + Rta[i] * convolutiona * 0.05 + Rtb[i] * convolutionb * 0.2 + Rtc[i] * convolutionc * 0.6;
      
    }
    
  }
}

model {

  weekly_sda ~ normal(0,1);
  weekly_sdb ~ normal(0,1);
  weekly_sdc ~ normal(0,1);
  
  weekly_rhoa ~ normal(0.8, 0.5);
  weekly_rhob ~ normal(0.8, 0.5);
  weekly_rhoc ~ normal(0.8, 0.5);
  
  weekly_effecta[1:(N)] ~ normal(weekly_rhoa[1:(N)] , weekly_sda);
  weekly_effectb[1:(N)] ~ normal(weekly_rhob[1:(N)] , weekly_sdb);
  weekly_effectc[1:(N)] ~ normal(weekly_rhoc[1:(N)] , weekly_sdc);
  
  casesa ~ poisson(predictiona);
  casesb ~ poisson(predictionb);
  casesc ~ poisson(predictionc);
}"
)




fit <- sampling(stan_model1,data=stan_data,iter=2000,warmup=1000,
                chains=5,thin=5,control = list(adapt_delta = 0.99, 
                                               max_treedepth = 30))


out <- rstan::extract(fit)

if(!dir.exists(paste0("/projectnb2/hrtgrp/zhenwei/simulation/v1/replicates/sim_", currind, "/"))){
  system(paste0("mkdir /projectnb2/hrtgrp/zhenwei/simulation/v1/replicates/sim_", currind, "/"))
}

saveRDS(out, paste0("/projectnb2/hrtgrp/zhenwei/simulation/v1/replicates/sim_", currind,"/out.rds"))

q()

library(matrixStats)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(dplyr)

data_list <- list()

i_read <- 1

for(i in 1:100){
  
  if(!file.exists(paste0("/projectnb2/hrtgrp/zhenwei/simulation/v1/replicates/sim_", i, "/out.rds"))){next()}
  
  print(i)
  out <- readRDS(paste0("/projectnb2/hrtgrp/zhenwei/simulation/v1/replicates/sim_", i, "/out.rds"))
  
  
  N_simulated_list <- readRDS("/projectnb2/hrtgrp/zhenwei/simulation/v1/N_simulated.rds")
  
  N_simulated_data <- N_simulated_list[[i]]
  
  stan_data<-list()
  stan_data$N  <- nrow(N_simulated_data)
  stan_data$casesa <- N_simulated_data$a
  # stan_data$casesa <- rep(2, 182)
  stan_data$casesb <- N_simulated_data$b
  # stan_data$casesb <- rep(2, 182)
  stan_data$casesc <- N_simulated_data$c
  # stan_data$N2 <- 60
  # 
  # stan_data$logra <- log(R_est_sim_data[, 1])
  # stan_data$logrb <- log(R_est_sim_data[, 2])
  # stan_data$logrc <- log(R_est_sim_data[, 3])
  
  w <- sapply(1:14, function(x){
    pgamma(x, 2, 0.5) - pgamma(x-1, 2, 0.5)
  })  
  
  w <- w/sum(w)
  w <- c(w, rep(0, nrow(N_simulated_data) - length(w)))
  
  stan_data$SI <- w
  
  
  data_a<-data.frame(
    x=N_simulated_data$t,
    y_real=stan_data$casesa,
    y = colMeans(out$predictiona),
    yl = colQuantiles(out$predictiona,probs=c(0.025)),
    yh = colQuantiles(out$predictiona,probs=c(0.975)),
    Rt = colMeans(out$Rta),
    Rtl = colQuantiles(out$Rta,probs=c(0.025)),
    Rth = colQuantiles(out$Rta,probs=c(0.975)),
    region = "a"
  )
  
  data_b<-data.frame(
    x=N_simulated_data$t,
    y_real=stan_data$casesb,
    y = colMeans(out$predictionb),
    yl = colQuantiles(out$predictionb,probs=c(0.025)),
    yh = colQuantiles(out$predictionb,probs=c(0.975)),
    Rt = colMeans(out$Rtb),
    Rtl = colQuantiles(out$Rtb,probs=c(0.025)),
    Rth = colQuantiles(out$Rtb,probs=c(0.975)),
    region = "b"
  )
  
  
  data_c<-data.frame(
    x=N_simulated_data$t,
    y_real=stan_data$casesc,
    y = colMeans(out$predictionc),
    yl = colQuantiles(out$predictionc,probs=c(0.025)),
    yh = colQuantiles(out$predictionc,probs=c(0.975)),
    Rt = colMeans(out$Rtc),
    Rtl = colQuantiles(out$Rtc,probs=c(0.025)),
    Rth = colQuantiles(out$Rtc,probs=c(0.975)),
    region = "c"
  )
  
  
  data <- rbind(data_a, data_b, data_c)
  
  data$replicates <- i
  
  data_list[[i_read]] <- data
  
  i_read <- i_read + 1
  
}

data_all <- do.call("rbind", data_list)

# 
# data2_a<-data.frame(
#   x=c(1:stan_data$N)[1:stan_data$N2],
#   ex = colMeans(out$mua)[1:stan_data$N2],
#   exl = colQuantiles(out$mua,probs=c(0.025))[1:stan_data$N2],
#   exh = colQuantiles(out$mua,probs=c(0.975))[1:stan_data$N2],
#   region = "a"
# )
# 
# 
# data2_b<-data.frame(
#   x=c(1:stan_data$N)[1:stan_data$N2],
#   ex = colMeans(out$mub)[1:stan_data$N2],
#   exl = colQuantiles(out$mub,probs=c(0.025))[1:stan_data$N2],
#   exh = colQuantiles(out$mub,probs=c(0.975))[1:stan_data$N2],
#   region = "b"
# )
# 
# data2_c<-data.frame(
#   x=c(1:stan_data$N)[1:stan_data$N2],
#   ex = colMeans(out$muc)[1:stan_data$N2],
#   exl = colQuantiles(out$muc,probs=c(0.025))[1:stan_data$N2],
#   exh = colQuantiles(out$muc,probs=c(0.975))[1:stan_data$N2],
#   region = "c"
# )
# 
# 
# data2 <- rbind(data2_a, data2_b, data2_c)
# data_all_summarise <- data_all %>% group_by(., x) %>% summarise(., y = mean(y), 
#                                                                          yl = quantile(y, 0.025),
#                                                                          yu = quantile(y, 0.975),
#                                                                          Rt = mean(Rt),
#                                                                          Rtl = quantile(Rt, 0.025),
#                                                                          Rtu = quantile(Rt, 0.975))
# 
data_grouped <- data_all %>% group_by(., x, region)

data_all_summarise <- data_grouped  %>% summarise(., y_mean = mean(y),
                                                  y_real = mean(y_real),
                                                  yl = mean(yl),
                                                  yh = mean(yh),
                                                  Rt_mean = mean(Rt),
                                                  Rtl = mean(Rtl),
                                                  Rth = mean(Rth),
                                                  Rt_mean_test = mean(Rt>1))

saveRDS(data_all_summarise, "/projectnb2/hrtgrp/zhenwei/simulation/data_all_summarise_v1.rds")

g1 <- ggplot(data_all_summarise, aes(x = x, color = region)) +
  # geom_smooth(data,mapping=aes(y=y,x=x, color = region),size=1.2) +
  geom_line(data_all_summarise,mapping=aes(y=y_mean,x=x, color = region),size=1.2) +
  # geom_bar(data = data, aes(y = y_real, fill = region), stat='identity', alpha=0.5) +
  geom_ribbon(data_all_summarise,mapping=aes(x=x,ymin=yl,ymax=yh, fill = region),alpha=0.3) + 
  theme_classic() +  xlab("Days") + ylab('Infections')

g2 <- ggplot(data_all_summarise, aes(x = x, color = region)) +
  geom_line(data_all_summarise,mapping=aes(y=Rt_mean,x=x, color = region),size=1.2)+ 
  coord_cartesian(ylim = c(0, 5)) +
  # lims(y = c(0, 5)) +
  # geom_line(data,mapping=aes(y=Rt,x=x, color = region),size=1.2)+lims(y = c(0, 5)) +
  geom_ribbon(data_all_summarise,mapping=aes(x=x,ymin=Rtl,ymax=Rth, fill = region),alpha=0.4) + 
  geom_hline(yintercept = 1,color = "black", size=1,alpha=0.5)  +
  theme_classic() +  xlab("Days") + ylab('Reproduction Number')# + scale_y_continuous(trans = "log10")


# 
# g3 <- ggplot(data2, aes(x = x, color = region)) +
#   geom_line(data2,mapping=aes(y=ex,x=x, color = region),size=1.2) +
#   geom_ribbon(data2,mapping=aes(x=x,ymin=exl,ymax=exh, fill = region),alpha=0.4) + 
#   theme_classic() +  xlab("Days") + ylab('Imported Infections')

grid.arrange(g1,g2#,
             # g3
)


jpeg("/projectnb2/hrtgrp/zhenwei/simulation/rt_v1.jpg", res = 300, w = 3000, h = 2000)
grid.arrange(g1,g2)
dev.off()



jpeg("/projectnb2/hrtgrp/zhenwei/simulation/rt_v1.jpg", res = 300, w = 3000, h = 2000)
grid.arrange(g1,g2)
dev.off()


library(Cairo)
library(ggpubr)


ggsave(ggarrange(g1,g2, labels = c("A", "B"), ncol = 1), file="/projectnb2/hrtgrp/zhenwei/simulation/rt_v1.eps", w = 8, h = 7, units = "in", device = cairo_ps)



# 
# g1 <- ggplot(data_all, aes(y=y, x = x, group = interaction(region, replicates), color = region)) + geom_line()
# geom_line(data,mapping=aes(y=y,x=x, group = interaction(region, replicates), color = region),size=1.2) +
#   # geom_bar(data = data, aes(y = y_real, fill = region), stat='identity', alpha=0.5) +
#   geom_ribbon(data,mapping=aes(x=x,ymin=yl,ymax=yh, fill = region),alpha=0.3) + 
#   theme_classic() +  xlab("Days") + ylab('Infections')