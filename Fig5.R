# R code to generate the figures of the manuscript
#
# Birth of immunity in early life: stepwise emergence of resistance and
# immunity to complex molecular parasites via progressively degenerating
# hyperparasites
#
# Copyright (c) 2023, Christian Iseli, Bernard Conrad, Joseph A. Curran and Magnus Pirovino
#

# From Model 1(1+) to 1(1-) ∶ transition from hyperparasitism to antiparasitism: show stepwise transition in one (or more) pictures

dt <- 0.25
m <- 4000
a <- 0.05
b <- 0.06
c <- 0.07
d_0 <- 0.05
e <- 0.1
theta_P_0 <- 0.04
theta_F2_0 <- 0.04
alpha <- 0.005 #  inflow rate of parasite P into habitat
beta <- 0.1 # trigger parameter of replicase R2 producing F1 upon encounter with P
delta_1 <- 0.07
delta_2 <- 0.0005
delta_3_0 <- 0.05
transition_steps <- 4
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

pdf('Fig5.pdf',paper="a4",width=0,height=0)
op <- par(mfrow = c(transition_steps+1,1), mgp = c(2,.8,0), mar = 0.1+c(3,3,3,1))

for (j in 0:transition_steps) {
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

  # transition parameters
  eta      <- j / transition_steps
  d        <- d_0        * (1 - eta)
  delta_3  <- delta_3_0  * (1 - eta)
  theta_P  <- theta_P_0  * eta
  theta_F2 <- theta_F2_0 * eta

  # Construct more rows, based on the last one

  for (k in 1:m) {
    Pi_uppercase <- (1 - (r$R1[k] + r$R2[k] + r$P[k] / q_P + r$F1[k] / q_F1 + r$F2[k] / q_F2) / K)
    R2 <- r$R2[k] + dt * ((a * r$R1[k] * r$R1[k] - delta_2 * r$R1[k] *r$F1[k]) * Pi_uppercase                                       - d_R2 * r$R2[k])
    R1 <- r$R1[k] + dt * ((b * r$R2[k] * r$R1[k] - delta_1 * r$R2[k] *r$P[k])  * Pi_uppercase                                       - d_R1 * r$R1[k])
    P  <- r$P[k]  + dt * ((c * r$R2[k] * r$P[k]  - delta_3 * r$F2[k] *r$P[k])  * Pi_uppercase + alpha - theta_P  * r$F2[k] * r$P[k] - d_P  * r$P[k])
    F1 <- r$F1[k] + dt * ((d * r$F2[k] * r$P[k]  + beta    * r$R2[k] *r$P[k])  * Pi_uppercase                                       - d_F1 * r$F1[k])
    F2 <- r$F2[k] + dt * ( e * r$R1[k] * r$F1[k]                               * Pi_uppercase         - theta_F2 * r$F2[k] * r$P[k] - d_F2 * r$F2[k])
    t  <- r$t[k] + dt

    # create the full row and add
    y <- cbind(t, R2, R1, P, F1, F2, deparse.level = 0)
    colnames(y) <- c("t","R2","R1","P","F1","F2")
    r <- rbind(r, y)
  }
  print(paste("R1 is", R1, "R2 is", R2))

  # plot the computed data
  plot(r$t,r$R1,col="green", main=paste0("Model 1(1+) to 1(1-) step ",j," (",eta,")"), xlab="time", ylab="population size", pch=".", ylim=c(0,10), xlim=c(0,1300))
  points(r$t,r$R2,col="cyan", pch=".")
  points(r$t,r$P,col="blue", pch=".")
  points(r$t,r$F1,col="red", pch=".")
  points(r$t,r$F2,col="pink", pch=".")
  lines(x=c(0,1000), y=rep(4.26,2), col="grey", lty=2, lwd=1)
  legend("topright", c("R1 host ribozyme", "R2 host replicase", "P parasite", "F1 prototemplate", "F2 protoreplicase/antibody", "population size = 4.26"),
         lty = c(1,1,1,1,1,2), col = c("green","cyan","blue","red","pink","grey"), box.lwd = 0)
}

dev.off()
#
