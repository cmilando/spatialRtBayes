#
library(splines)
library(tidyverse)
source("R/00_generate_data.R")

y <- R_this[, 2]
ytrue = R_rev[, 2]

x <- 1:length(y)
N_KNOTS = 4 ## inversely proportional to the length of y
k <- quantile(x, probs = seq(from = 0, to = 1, length.out  = N_KNOTS))

b1 <- lm(y ~ bs(x, knots = k, degree = 3, intercept = F))
x.pred <- predict(b1, se = T)

head(x.pred$se.fit)
plot(x, y, type = 'l')
lines(x, x.pred$fit, col = 'red')
lines(x, x.pred$fit - 1.96 * x.pred$se.fit, col = 'red')
lines(x, x.pred$fit + 1.96 * x.pred$se.fit, col = 'red')
lines(x, ytrue, col = 'green')

## the thing is, you don't know what the true data generating mechanism is
## so, you would want it to choose the most
## but you could set some params based on your thoughts
## for Number of Knots and Spacing between them

N_KNOTS = 4 ## inversely proportional to the length of y
k <- quantile(x, probs = c(0.1, 0.2, 0.3, 0.4))

b1 <- lm(y ~ bs(x, knots = k, degree = 3, intercept = F))
x.pred <- predict(b1, se = T)

# head(x.pred$se.fit)
# plot(x, y)
lines(x, x.pred$fit, col = 'blue')
lines(x, x.pred$fit - 1.96 * x.pred$se.fit, col = 'blue')
lines(x, x.pred$fit + 1.96 * x.pred$se.fit, col = 'blue')
