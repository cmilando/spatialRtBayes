
## simulate base on individual data
set.seed(1)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(Cairo)
m <- 3 # number of regions
pop <- rnorm(m, 50000, 10000) %>% round()
initial_inf_rate <- rbeta(m, 1, 500)
initial_inf_num <- round(pop*initial_inf_rate)

initial_inf_num <- c(10, 10, 10)
N_day <- 200
w <- sapply(1:14, function(x){
  pgamma(x, 2, 0.5) - pgamma(x-1, 2, 0.5)
})  

w <- w/sum(w)

P <- matrix(c(0.8, 0.2, 0.1, 
              0.15, 0.6, 0.3,
              0.05, 0.2, 0.6), 3)



x <- 1:214

y1 <- (20*cos(x/500) + ((0.8*x-50))^2 - (x*0.115)^3)/1000 + 0.8

plot(x, y1)

y2 <- (30*sin(x/150) + cos(x/20) - (x/50)^2)/8 - x*0.006

plot(x, y2)

y3 <- (30*cos(x/150) + 2*sin(x/20) + 2*(x/50)^2)/20 - x*0.005

plot(x, y3)


R_matrix_dataframe <- data.frame(a = y1, b = y2, c = y3)

R_matrix <- as.matrix(R_matrix_dataframe)

R_matrix_dataframe$t <- 1:nrow(R_matrix_dataframe)

R_matrix_dataframe_long <- gather(R_matrix_dataframe , region, R_t, a:c, factor_key=TRUE)


g1 <- ggplot(data = R_matrix_dataframe_long, aes(x = t, y = R_t, color = region)) + geom_line()+ xlab("Day") + ylab("Reproductive Number")+ theme_classic()


### theoretical Ns
t_max= 214

N_simulated <- matrix(NA, t_max + 1, m)

N_simulated[1, ] <- initial_inf_num

for(t in 1:t_max){
  
  N_new <- rev(w[1:min(length(w), t)])%*%N_simulated[c(max((t-(length(w)-1)), 1): t), ,drop = F]%*%diag(R_matrix[t, ])%*%P
  
  N_simulated[t+1, ] <- N_new
}
N_simulated <- as.data.frame(N_simulated)
names(N_simulated) <- c("a", "b", "c")
N_simulated$t <- 0:t_max

plot_dat_long <- gather(N_simulated, region, incidence, a:c, factor_key=TRUE)

g2 <- ggplot(data = plot_dat_long, aes(x = t, y = incidence, color = region)) + geom_line() + xlab("Day") + ylab("Incidence") + theme_bw()


grid.arrange(g2,g1)



## theoretical R

N_theoretical <- N_simulated[, -c(4, 5)] %>% as.matrix()

R_est_list <- lapply(1:t_max, function(t){
  (rev(w[1:min(length(w), t)])%*%N_theoretical[c(max(t-length(w) + 1, 1): t), ,drop = F])^(-1)*N_theoretical[t+1,]%*%solve(P)
})

R_est <- do.call("rbind", R_est_list)

R_est <- as.data.frame(R_est)

R_est$t <- 1:t_max

plot_dat_long <- gather(R_est, region, incidence, a:c, factor_key=TRUE)

ggplot(data = plot_dat_long, aes(x = t, y = incidence, color = region)) + geom_line() + theme_bw()

### simulated Ns
t_max= 214

# list of 100 simulated datasets
N_simulated_list <- list()
for(index in 1:100){
  
  N_simulated <- matrix(initial_inf_num, 1)
  
  for(t in 1:t_max){
    
    cases_expectation <- rev(w[1:min(length(w), t)])%*%N_simulated[c(max(t-length(w) + 1, 1): t), ,drop = F]%*%diag(R_matrix[t, ])%*%P
    
    
    cases <- sapply(cases_expectation, function(x){rpois(1, x)})
    
    N_simulated <- rbind(N_simulated, cases)
  }
  N_simulated <- as.data.frame(N_simulated)
  names(N_simulated) <- c("a", "b", "c")
  N_simulated$t <- 0:t_max
  N_simulated$sim <- index
  N_simulated_list[[index]] <- N_simulated
}

# save the first one as example_data.csv
write.csv(N_simulated_list[1], "example_data.csv", row.names = F)
