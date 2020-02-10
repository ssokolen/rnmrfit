library(rnmrfit)

nmroptions$direct$sf <- 500
nmroptions$indirect$sf <- 500

r1 <- nmrresonance_2d('2 d 10', '2 d 10')
r2 <- nmrresonance_2d('2.5 d 10', '2 s')

s <- nmrspecies_2d(list(r1, r2))

x1 <- seq(1.5, 3, length.out = 100)
x2 <- seq(1.5, 3, length.out = 100)
p <- expand.grid(x1 = x1, x2 = x2)

p <- tibble(direct.shift = p$x1, indirect.shift = p$x2,
            intensity = values(s, p$x1, p$x2))

print(summary(p$intensity))


d <- new("NMRData2D")
d@processed <- p

plot(d)

print(s)
