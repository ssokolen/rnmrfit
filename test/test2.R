library(rnmrfit)

nmroptions$direct$sf <- 500
nmroptions$indirect$sf <- 500

r1 <- nmrresonance_2d('2 d 1', '2 d 2')
r2 <- nmrresonance_2d('3 d 1', '3 d 2')

s <- nmrspecies_2d(list(r1, r2))

print(s)
