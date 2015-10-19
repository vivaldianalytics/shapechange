library(abind)
data(hoverflies)
A = hoverflies

B = abind(A-sqrt(A),A + sqrt(A),along = 3)
B[,15,1] = A$year
B[,15,2] = A$year
B[B<1] = 1
show(B)
hoverflies_INT = B
save(hoverflies_INT, file = "data/hoverflies_INT.rda")