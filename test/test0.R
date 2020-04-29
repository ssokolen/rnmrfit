library(rnmrfit)

nmroptions$sf <- 500

# Generating data
r1 <- nmrresonance_1d('1.5 s')
r2 <- nmrresonance_1d('1.6 s')

s <- nmrspecies_1d(list(r1, r2))

n <- 2000
x <- seq(1, 2, length.out = n)
p <- tibble(direct.shift = x,
            intensity = values(r1, x) +
                        values(r2, x) +
                        cmplx1(r = rnorm(n, 0, .005), i = rnorm(n, 0, .005)))

d <- new("NMRData1D")
d@processed <- p

# Fitting data
r1 <- nmrresonance_1d('1.48 s')
r2 <- nmrresonance_1d('1.58 s')
s <- nmrspecies_1d(c(r1, r2))

f <- nmrfit_1d(s, d)

