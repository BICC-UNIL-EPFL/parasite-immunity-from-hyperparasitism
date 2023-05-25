# R code to generate the figures of the manuscript
#
# Birth of immunity in early life: stepwise emergence of resistance and
# immunity to complex molecular parasites via progressively degenerating
# hyperparasites
#
# Copyright (c) 2023, Christian Iseli, Bernard Conrad, Joseph A. Curran and Magnus Pirovino
#

# Model 1(0) habitat with primordial autocatalytic cycle, exposed to template invasion P
# and invasion of protoenzyme E_1

dt    <-    0.2
m     <- 5000
a     <-    0.1
b     <-    0.02
alpha <-    0.005 # inflow rate of RNA-template P into habitat
beta  <-    0.01  # inflow rate of protoenzyme/catalyst E1 into habitat
d_E1  <-    0.01 
d_E2  <-    0.01 
d_P   <-    0.01
d_F1  <-    0.01   
q_F1  <-   10.0
K     <-   10.0

# initial conditions setting:
E2_0 <- 0.5
E1_0 <- 0.0
P_0  <- 0.0
F1_0 <- 0.0

pdf('Fig2.pdf',paper="a4r",width=0,height=0)
op <- par(mfrow = c(1,1), mgp = c(2,.8,0), mar = 0.1+c(3,3,3,1))

# Construct first row of results

r <- as.data.frame(matrix(
  c(0.0,
    E2_0,
    E1_0,
    P_0,
    F1_0),
  nrow=1, ncol=5, byrow=TRUE))

colnames(r) <- c("t","E2","E1","P","F1")

# Construct more rows, based on the last one

for (k in 1:m) {
  Pi_uppercase <- (1 - (r$E2[k] + r$F1[k] / q_F1) /K)
  P  <- r$P[k]  + dt * (alpha - d_P * r$P[k])
  E1 <- r$E1[k] + dt * (beta - d_E1 * r$E1[k])
  E2 <- r$E2[k] + dt * (a * r$E1[k] * r$F1[k] * Pi_uppercase - d_E2 * r$E2[k])
  F1 <- r$F1[k] + dt * (b * r$E2[k] * r$P[k]  * Pi_uppercase - d_F1 * r$F1[k])
  t  <- r$t[k] + dt

  # create the full row and add
  y <- cbind(t, E2, E1, P, F1, deparse.level = 0)
  colnames(y) <- c("t","E2","E1","P","F1")
  r <- rbind(r, y)
}

# plot the computed data
plot(r$t,r$E1,col="green", main="Model 1(0)", xlab="time", ylab="population size", pch=".", ylim=c(0,7))
points(r$t,r$E2,col="cyan", pch=".")
points(r$t,r$P,col="blue", pch=".")
points(r$t,r$F1,col="grey", pch=".")
legend("topleft", c("E1 invaded protoenzyme/catalyst", "P invaded RNA-template", "E2 protoreplicase", "F1 prototemplate"), lty = 1, col = c("green","blue","cyan","grey"), box.lwd = 0)

dev.off()
#
