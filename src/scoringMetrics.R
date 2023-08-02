require(SID)
require(qgraph)


G <- matrix(data = c(0, 1, 0,
                     0, 0, 0,
                     1, 1, 0),
            nrow = 3,
            byrow = T)

E0 <- matrix(data = c(0, 0, 0,
                      0, 0, 0,
                      0, 0, 0),
             nrow = 3,
             byrow = T)

E1 <- matrix(data = c(0, 1, 1,
                     0, 0, 0,
                     0, 1, 0),
            nrow = 3,
            byrow = T)

E2 <- matrix(data = c(0, 1, 1,
                     0, 0, 0,
                     1, 1, 0),
            nrow = 3,
            byrow = T)

E3 <- matrix(data = c(0, 1, 0,
                     0, 0, 0,
                     0, 0, 0),
            nrow = 3,
            byrow = T)

print("SHD")
for (e in list(G,E0,E1,E2,E3, t(G))) {
  print(hammingDist(G,e, allMistakesOne = T))
}

print("SID")
for (e in list(G,E0, E1,E2,E3, t(G))) {
  print("upper")
  print(structIntervDist(G,e)$sidUpperBound)
  print("lower")
  print(structIntervDist(G,e)$sidLowerBound)
}



####explaining SID

Tr <- t(matrix(data = c(0, 1, 0, 0 ,1,
                     0, 0, 1, 0 ,0,
                     0, 0, 0, 0, 0,
                     0, 0, 1, 0, 0,
                     0, 0, 0, 0, 0
                     ),
            nrow = 5,ncol=5,
            byrow = T))
qgraph(Tr)

G1 <- t(matrix(data = c(0, 1, 0, 0 ,0,
                           0, 0, 1, 0 ,0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0,
                           1, 0, 0, 0, 0
),
nrow = 5,
byrow = T))
qgraph(G1)

s <- structIntervDist(Tr, G1)

par(mfrow=c(1,3))
a <- qgraph(Tr)
b <- qgraph(G1, layout=a$layout)
c <- qgraph(s$incorrectMat, layout=a$layout)
