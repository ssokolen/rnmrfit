library(rnmrfit)

nmroptions$sf <- 500

# Generating data
r1 <- nmrresonance_1d('1 d 10')
r2 <- nmrresonance_1d('1.5 t 10')

s <- nmrspecies_1d(list(r1, r2), areas = c(1, 2))

n <- 2000
x <- seq(0.5, 2, length.out = n)
p <- tibble(direct.shift = x,
            intensity = values(r1, x) +
                        2*values(r2, x) +
                        cmplx1(r = rnorm(n, 0, .005), i = rnorm(n, 0, .005)))

d <- new("NMRData1D")
d@processed <- p

# Fitting data
r1 <- nmrresonance_1d('0.995 d 8', position.leeway = 0.4)
r2 <- nmrresonance_1d('1.491 t 8', position.leeway = 0.4)

s <- nmrspecies_1d(list(r1, r2), areas = c(1, 2))

f <- nmrfit_1d(s, d)

