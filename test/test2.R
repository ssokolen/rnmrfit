library(rnmrfit)

nmroptions$direct$sf <- 500
nmroptions$indirect$sf <- 500

r1 <- nmrresonance_2d('1.2 d 10', '1.2 d 10', width = 4)
r2 <- nmrresonance_2d('1.4 d 10', '1.4 s', width = 4)

s <- nmrspecies_2d(list(r1, r2))

x1 <- seq(1, 1.5, length.out = 100)
x2 <- seq(1, 1.5, length.out = 100)
p <- expand.grid(x1 = x1, x2 = x2)

intensity <- values(s, p$x1, p$x2, components = 'rr/ri/ir/ii')

p <- tibble(direct.shift = p$x1, indirect.shift = p$x2,
            intensity = intensity)

#print(summary(p$intensity))


d <- new("NMRData2D")
d@processed <- p

plot(d)

#print(s)
