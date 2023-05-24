# R code to generate the figures of the manuscript
#
# Birth of immunity in early life: stepwise emergence of resistance and
# immunity to complex molecular parasites via progressively degenerating
# hyperparasites
#
# Copyright (c) 2023, Christian Iseli, Bernard Conrad, Joseph A. Curran and Magnus Pirovino
#

# Model 1(1-): transition from hyperparasitism to antiparasitism: now pure antiparasitism

dt <- 0.25
m <- 5000
a <- 0.05
b <- 0.06
c <- 0.07
d <- 0.00 # protoreplicase F2 loses its ability to reproduce prototemplate R1 upon encounter with P
e <- 0.1
theta_P <- 0.04
theta_F2 <- 0.04
alpha <- 0.005 #  inflow rate of parasite P into habitat
beta <- 0.1 # trigger parameter of replicase R2 producing F1 upon encounter with P
delta_1 <- 0.07
delta_2 <- 0.0005
delta_3 <- 0.00 # protoreplicase F2 loses its ability to parasitize template P upon encounter with P
d_R2 <- 0.01 
d_R1 <- 0.01 
d_P <- 0.01
d_F1 <- 0.01
d_F2 <- 0.01   
q_P <- 2
q_F1 <- 10.0
q_F2 <- 10.0
K <- 10.0

# initial conditions setting:
R2_0 <- 1.0
R1_0 <- 3.0
P_0 <- 0.0
F1_0 <- 0.0
F2_0 <- 0.0

pdf('Fig4B.pdf',paper="a4r",width=0,height=0)
op <- par(mfrow = c(1,1), mgp = c(2,.8,0), mar = 0.1+c(3,3,3,1))

# Construct first row of results

r <- as.data.frame(matrix(
  c(0.0,
    R2_0,
    R1_0,
    P_0,
    F1_0,
    F2_0),
  nrow=1, ncol=6, byrow=TRUE))

colnames(r) <- c("t","R2","R1","P","F1","F2")

# Construct more rows, based on the last one

for (k in 1:m) {
  Pi_uppercase <- (1 - (r$R1[k] + r$R2[k] + r$P[k] / q_P + r$F1[k] / q_F1 + r$F2[k] / q_F2) / K)
  R2 <- r$R2[k] + dt * ((a * r$R1[k] *r$R1[k] - delta_2 * r$R1[k] *r$F1[k]) * Pi_uppercase                                       - d_R2 * r$R2[k])
  R1 <- r$R1[k] + dt * ((b * r$R2[k] *r$R1[k] - delta_1 * r$R2[k] *r$P[k])  * Pi_uppercase                                       - d_R1 * r$R1[k])
  P  <- r$P[k]  + dt * ((c * r$R2[k] *r$P[k]  - delta_3 * r$F2[k] *r$P[k])  * Pi_uppercase + alpha - theta_P  * r$F2[k] * r$P[k] - d_P  * r$P[k])
  F1 <- r$F1[k] + dt * ((d * r$F2[k] *r$P[k]  + beta    * r$R2[k] *r$P[k])  * Pi_uppercase                                       - d_F1 * r$F1[k])
  F2 <- r$F2[k] + dt * ( e * r$R1[k] *r$F1[k]                               * Pi_uppercase         - theta_F2 * r$F2[k] * r$P[k] - d_F2 * r$F2[k])
  t  <- r$t[k] + dt

  # create the full row and add
  y <- cbind(t, R2, R1, P, F1, F2, deparse.level = 0)
  colnames(y) <- c("t","R2","R1","P","F1","F2")
  r <- rbind(r, y)
}

# plot the computed data
plot(r$t,r$R1,col="green", main="Model 1(1-)", xlab="time", ylab="population size", pch=".", ylim=c(0,6))
points(r$t,r$R2,col="cyan", pch=".")
points(r$t,r$P,col="blue", pch=".")
points(r$t,r$F1,col="red", pch=".")
points(r$t,r$F2,col="pink", pch=".")
legend("topright", c("R1 host ribozyme", "R2 host replicase", "P parasite", "F1 prototemplate", "F2 antibody"), lty = 1, col = c("green","cyan","blue","red","pink"), box.lwd = 0)

dev.off()
#
