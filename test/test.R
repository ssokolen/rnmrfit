library(rnmrfit)

n <- as.integer(5)
n.peaks <- as.integer(4)
n.baseline <- as.integer(0)
n.phase <- as.integer(0)

x <- as.double(seq(0, n))
y <- as.double(rnorm(n*2))
par <- as.double(runif(n.peaks))
lb <- as.double(rep(0, n.peaks))
ub <- as.double(rep(1, n.peaks))

out <- .Call("fit_1d_wrapper", x = x, y = y, par = par, lb = lb, ub = ub,
             n = n, nl = n.peaks, nb = n.baseline, np = n.phase, basis = x)

print(out)
