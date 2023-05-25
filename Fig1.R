# R code to generate the figures of the manuscript
#
# Birth of immunity in early life: stepwise emergence of resistance and
# immunity to complex molecular parasites via progressively degenerating
# hyperparasites
#
# Copyright (c) 2023, Christian Iseli, Bernard Conrad, Joseph A. Curran and Magnus Pirovino
#

# Strong and weak parasite invasion

dt <- 0.3
m <- 6000
a <- 0.08
b <- 0.11
c <- 0.05
delta <- 0.001
d_R2 <- 0.01 
d_R1 <- 0.01 
d_P <- 0.01   
q_P <- 2.0
K <- 10.0

pdf('Fig1.pdf',paper="a4r",width=0,height=0)
op <- par(mfrow = c(2,1), mgp = c(2,.8,0), mar = 0.1+c(3,3,3,1))

for (alpha in c(0.1, 0.01)) {   # high and low parasite invasion constant 
  # state of alpha
  state <- ifelse(alpha == 0.1, "strong", "weak")

  # initial conditions

  R2_0 <- 0.5
  R1_0 <- 1.0
  P_0  <- 0.0

  # Construct first row of results

  r <- as.data.frame(matrix(
    c(0.0,
      R2_0,
      R1_0,
      P_0),
    nrow=1, ncol=4, byrow=TRUE))

  colnames(r) <- c("t","R2","R1","P")

  # Construct more rows, based on the last one

  for (k in 1:m) {
    Pi_uppercase <- (1 - (r$R1[k] + r$R2[k] + r$P[k] / q_P) / K)
    R2 <- r$R2[k] + dt * ( a * r$R1[k] * r$R1[k]                             * Pi_uppercase         - d_R2 * r$R2[k])
    R1 <- r$R1[k] + dt * ((b * r$R2[k] * r$R1[k] - delta * r$R2[k] * r$P[k]) * Pi_uppercase         - d_R1 * r$R1[k])
    P  <- r$P[k]  + dt * ( c * r$R2[k] * r$P[k]                              * Pi_uppercase + alpha - d_P  * r$P[k])
    t <- r$t[k] + dt

    # create the full row and add
    y <- cbind(t, R2, R1, P, deparse.level = 0)
    colnames(y) <- c("t","R2","R1","P")
    r <- rbind(r, y)
  }

  # plot the computed data
  plot(r$t, r$R2,col="cyan", main=paste0("Model 1(1) ",state," parasites (alpha = ",alpha,")"), xlab="time", ylab="population size", pch=".", ylim=c(0,16))
  points(r$t, r$R1,col="green", pch=".")
  points(r$t, r$P,col="blue", pch=".")
  legend("topright", c("R2 host replicase", "R1 host ribozyme", "P parasite"), lty = 1, col = c("cyan","green","blue"), box.lwd = 0)
}
dev.off()
#
