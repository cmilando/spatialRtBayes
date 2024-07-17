#### INPUT
start_day <- 12 ## day for start estimating
window_length <- 40 ## smoothing window for estimaiton
si_shape <- 2 ## shape parameter for serical interval assuming gamma distribution
si_rate <- 0.5 ## rate parameter for serical interval assuming gamma distribution
si_t_max <- 14 ## maximum number of days with non-zero probability for serial interval

# Observed cases, matrix Y
incidence_data <- readRDS("old_data/hrt-main/data/MA_all_counties_cases_count.rds")
dt1 = min(incidence_data$cdc_report_dt)
dt2 = max(incidence_data$cdc_report_dt)
seqq <- seq.Date(from = as.Date(dt1), to = as.Date(dt2), by = 'day')
length(seqq)
head(incidence_data)
incidence_data <- incidence_data %>% 
  select(county_fips_code, cdc_report_dt, n) %>%
  pivot_wider(id_cols = cdc_report_dt, names_from = county_fips_code,
              values_from = n, values_fill = 0) %>% arrange(cdc_report_dt)
head(incidence_data)
incidence_data <- left_join(data.frame(cdc_report_dt = seqq), incidence_data)
for(j in 2:ncol(incidence_data)) {
  dd <- which(is.na(incidence_data[, j]))
  incidence_data[dd, j] <- 0
}

dim(incidence_data)
head(incidence_data)

plot(incidence_data[, 1], incidence_data[, 2], type = 'l')
for(j in 3:ncol(incidence_data)) {
  lines(incidence_data[, 1], incidence_data[, j], type = 'l')
}

# just get the bottom
incidence_data <- incidence_data[150:362, ]
plot(incidence_data[, 1], incidence_data[, 2], type = 'l')
for(j in 3:ncol(incidence_data)) {
  lines(incidence_data[, 1], incidence_data[, j], type = 'l')
}

Y <- as.matrix(incidence_data[, -1])
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
# P <- matrix(c(0.8, 0.2, 0.1, 
#               0.15, 0.6, 0.3,
#               0.05, 0.2, 0.6), 3)
P <- readRDS("old_data/hrt-main/data/mobility_ma_list.rds")
head(P)
names(P)
length(P)
length(P[[1]])

# Serial interval, vector W
S <- as.integer(si_t_max)
## obtain discretized gamma distribution for serial interval
W <- sapply(1:S, function(x){
  pgamma(x, si_shape, si_rate) - pgamma(x-1, si_shape, si_rate)
})  

W <- W/sum(W) 
plot(W, type = 'l')

## -----------------------------------------------------------
## ok now get the matrices to solve for R
# M <- matrix(integer(), nrow = N, ncol = J)
# M[1, ] <- 10
# head(M)
# head(Y)

R_rev <- matrix(NA, nrow = N, ncol = J)

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
    sum_m_w[jx] = Y[start_index:end_index, jx] %*% W[1:(end_index - start_index + 1)] 
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
  # M[n, ] <- c_mat %*% Rx[, n]
  # M[n, ]
  # 
  ####
  # and then the classic OLS matrix inverse to get R
  # R = (C'C)^-1C'Y
  R_rev[n, ] <- solve(t(c_mat) %*% c_mat) %*% t(c_mat) %*% Y[n, ]
  
}

head(R_rev)

plot(x = 1:200, y = R_rev[1:200, 1], type = 'l', col = 'red')
lines(x = 1:200, y = R_rev[1:200, 2], type = 'l', col = 'green')
lines(x = 1:200, y = R_rev[1:200, 3], type = 'l', col = 'blue')

head(R_rev)

