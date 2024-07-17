# what are the boundaries of tau
tau_end = min(S, t - 1)

RR <- diag(Rmatrix[t, ])
RR

## MM is m(t - 1), ..., m(1)
## where rows are regions
## and columns are time
if(t == 2) {
  MM <- as.matrix(M[t - 1:tau_end, ])
} else {
  MM <- t(as.matrix(M[t - 1:tau_end, ]))
}
MM
WW <-  as.matrix(w[1:tau_end])

result <- t(P) %*% RR %*% MM %*% WW

