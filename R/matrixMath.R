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
head(Y)
N <- as.integer(nrow(Y))
J <- as.integer(ncol(Y))
#Y[1, ] <- 0  # Padding the first row as mentioned in your comments

# Transfer matrix, matrix P (random example)
P <- matrix(c(0.8, 0.2, 0.1, 
              0.15, 0.6, 0.3,
              0.05, 0.2, 0.6), 3)

# Serial interval, vector W
S <- as.integer(si_t_max)
## obtain discretized gamma distribution for serial interval
W <- sapply(1:S, function(x){
  pgamma(x, si_shape, si_rate) - pgamma(x-1, si_shape, si_rate)
})  

W <- W/sum(W) 
plot(W, type = 'l')

Ra <- function(t) (20*cos(t/500) + (0.8*t - 50)^2 - (0.115 * t)^3)/1000 + 0.8 
Rb <- function(t) (30*sin(t/150) + cos(t/20) - (t/50)^2)/8 - 0.006*t 
Rc <- function(t) (30*cos(t/150) + 2*sin(t/20) + 2*(t/50)^2)/20 - 0.005*t

plot(x = 1:200, y = Ra(1:200), type = 'l', col = 'red')
lines(x = 1:200, y = Rb(1:200), type = 'l', col = 'green')
lines(x = 1:200, y = Rc(1:200), type = 'l', col = 'blue')

## -----------------------------------------------------------
## ok now get the matrices to solve for R
M <- matrix(integer(), nrow = N, ncol = J)
M[1, ] <- 10
head(M)
head(Y)

R_rev <- matrix(NA, nrow = N, ncol = J)

n_R <- N
Rx <- matrix(c(Ra(1:n_R), Rb(1:n_R), Rc(1:n_R)), byrow = T, nrow = 3, ncol = n_R)
head(Rx)

for(n in 2:N) {
  print(n)
  
  # first, get the M*W dot product
  sum_m_w <- vector("numeric", J)
  start_index = max(1, n - S + 1)
  end_index = n - 1
  start_index
  end_index
  ## sum(m*w)
  for (jx in 1:J) {
    sum_m_w[jx] = M[start_index:end_index, jx] %*% W[1:(end_index - start_index + 1)] 
  }
  sum_m_w
  
  sum_m_w_mat <- matrix(rep(sum_m_w, times = J), byrow = T, nrow = J)
  sum_m_w_mat <- t(sum_m_w_mat)
  sum_m_w_mat
  
  # now the a, b, c ... matrix is really just P * rowwise the above
  c_mat <- t(P) * sum_m_w_mat
  c_mat <- t(c_mat)
  c_mat
  
  ####
  # now mult by R
  M[n, ] <- c_mat %*% Rx[, n]
  M[n, ]
  
  ####
  # and then the classic OLS matrix inverse to get R
  # R = (C'C)^-1C'Y
  R_rev[n, ] <- solve(t(c_mat) %*% c_mat) %*% t(c_mat) %*% M[n, ]
  
}

head(M)

plot(x = 1:200, y = M[1:200, 1], type = 'l', col = 'red')
lines(x = 1:200, y = M[1:200, 2], type = 'l', col = 'green')
lines(x = 1:200, y = M[1:200, 3], type = 'l', col = 'blue')

head(R_rev)
head(t(Rx))
