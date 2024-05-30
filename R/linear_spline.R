

# Define your data points
x <- c(1, 2, 3, 4, 5)[1:KNOTS]
y <- c(2, 3, 5, 4, 6)[1:KNOTS]

# so, what the the euqation is really choosing is Y
# Nah, you give it the number of points, and it creates it afterwards


# Create the design matrix

## so this is true for X in between kn[1] <= x <= kn[2]

# so first solve


# then project
ll <- vector("list", KNOTS - 1)
for(i in 1:length(ll)) {

  X <- cbind(1, x[c(i, i+1)])
  #print(X)
  X
  
  # Convert y to a column vector
  Y <- matrix(y[c(i, i+1)], ncol = 1)
  Y
  
  # Compute the coefficients using the normal equation
  t(X) %*% X
  
  beta <- solve(t(X) %*% X) %*% t(X) %*% Y
  
  # for npoints - 1
  beta
  
  # Extract the coefficients
  a <- beta[1, 1]
  b <- beta[2, 1]
  # cc <- beta[3, 1]
  # d <- beta[4, 1]
  
  # Predict new Y values using the linear model
  new_x <- seq(min(x[c(i, i+1)]), max(x[c(i, i+1)]), length.out = L_OUT)
  new_y <- a + b * new_x #+ cc * new_x^2 + d * new_x^3
  
  ll[[i]] <- data.frame(new_x = round(new_x, 3), new_y = round(new_y, 3))
}

ll_df <- do.call(rbind, ll)
ll_df <- unique(ll_df)


# Plot the original points and the linear fit
plot(x, y, main = "Linear Fit Using Matrix Math",
     xlab = "X", ylab = "Y", pch = 19)
points(ll_df$new_x, ll_df$new_y, col = "red")
#lines(ll_df$new_x, ll_df$new_y, col = "red")
