library("mvtnorm")

# Enter parameters of the model.

M <- 4
C <- c(-Inf,-0.1,0,0.1,Inf)

X <- c(-1,1)
G <- c(-Inf,0,Inf)

Z <- c(-1,-7/9,-5/9,-3/9,-1/9,1/9,3/9,5/9,7/9,1)

alpha_1 <- 1
alpha_2 <- 0.5

Mu <- c(0,0)
Sigma <- matrix(c(1,0.5,0.5,1),nrow=2,ncol=2)

# Enter grid of values over which you wish to search for the identified set.

a_0 <- seq(-1 , 1 , length.out = 5)
a_1 <- seq(0 , 2 , length.out = 5)

# The following code will create a data.frame of all possible combinations of (a_0 , a_1 , Y , X , Z). 

Y <- seq(1 , M , 1)
Pr.frame <- expand.grid(Y , X , Z , a_0 , a_1 )
colnames(Pr.frame) <- c("Y" , "X" , "Z" , "a_0" , "a_1")

# We now need to calculate the probability of the event (Y=m and X=x | Z=z) using the mvtnorm package. 

Pr <- vector(length = nrow(Pr.frame))
for(i in seq(1,nrow(Pr.frame),1)){
  Pr[i] <- pmvnorm(lower = c(C[match(Pr.frame$Y[i] , Y)] - alpha_1*Pr.frame$X[i] , G[match(Pr.frame$X[i] , X)] - 0.5*Pr.frame$Z[i]) , upper = c(C[match(Pr.frame$Y[i] , Y) + 1] - alpha_1*Pr.frame$X[i] , G[match(Pr.frame$X[i] , X) + 1] - 0.5*Pr.frame$Z[i]) , mean = Mu , sigma = Sigma) 
}

# Combine these probabilities with the data.frame.

Pr.frame$Pr <- Pr

# Cumulate these probabilities for Y<=m.

cum.Pr <- vector(length = nrow(Pr.frame))
for(i in seq(1,nrow(Pr.frame),1)){
  cum.Pr[i] <- sum(Pr.frame[Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$X[i] == Pr.frame$X & Pr.frame$Y[i] >= Pr.frame$Y & Pr.frame$a_1[i] == Pr.frame$a_1 & Pr.frame$a_0[i] == Pr.frame$a_0 , "Pr"])
}

# Combine these probabilities with the data.frame.

Pr.frame$cum.Pr <- cum.Pr

# Compute parameter values.

Par_L <- vector(length = nrow(Pr.frame)) 
Par_U <- vector(length = nrow(Pr.frame))
for(i in seq(1,nrow(Pr.frame),1)){
  Par_L[i] <- C[match(Pr.frame$Y[i] , Y)] - Pr.frame$a_0[i] - Pr.frame$a_1[i]*Pr.frame$X[i]
  Par_U[i] <- C[match(Pr.frame$Y[i] , Y) + 1] - Pr.frame$a_0[i] - Pr.frame$a_1[i]*Pr.frame$X[i]
}

# Combine these parameters with the data.frame.

Pr.frame$Par_L <- Par_L
Pr.frame$Par_U <- Par_U

# Normalise parameters such that they fall in the [0,1] interval.

Pr.frame <- transform(Pr.frame , Par_L=pnorm(Par_L))
Pr.frame <- transform(Pr.frame , Par_U=pnorm(Par_U))

# The first condition for a combination (a_0 , a_1) lying inside the identified set to be tested.

output.1 <- vector(length = nrow(Pr.frame))
for(i in seq(1,nrow(Pr.frame),1)){
  output.1[i] <- sum(Pr.frame[Pr.frame$Par_U[i] >= Pr.frame$Par_U & Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$a_0[i] == Pr.frame$a_0 & Pr.frame$a_1[i] == Pr.frame$a_1 & Pr.frame$Y != M, "Pr"])
} 
condition.1 <- ifelse(output.1 <= Pr.frame$Par_U , TRUE , FALSE)

# The second condition for a combination (a_0 , a_1) lying inside the identified set to be tested.

output.2 <- vector(length = nrow(Pr.frame))
for(i in seq(1,nrow(Pr.frame),1)){
  output.2[i] <- sum(Pr.frame[Pr.frame$Par_U[i] > Pr.frame$Par_L & Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$a_0[i] == Pr.frame$a_0 & Pr.frame$a_1[i] == Pr.frame$a_1 , "Pr"])
} 
condition.2 <- ifelse(output.2 >= Pr.frame$Par_U | Pr.frame$Y == M , TRUE , FALSE)

# The third condition for a combination (a_0 , a_1) lying inside the identified set to be tested.

output.3 <- vector(length = nrow(Pr.frame))
for(i in seq(1,nrow(Pr.frame),1)){
  output.3[i] <- ifelse(Pr.frame$Par_U[i] - Pr.frame$cum.Pr[i] >= 0, TRUE , FALSE)
}

output.4 <- vector(length = nrow(Pr.frame))
for(i in seq(1,nrow(Pr.frame),1)){
  output.4[i] <- ifelse(any(ifelse(
    (Pr.frame$Par_U[i] - min(Pr.frame[Pr.frame$Y[i] > Pr.frame$Y & Pr.frame$X[i] == Pr.frame$X & Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$a_0[i] == Pr.frame$a_0 & Pr.frame$a_1[i] == Pr.frame$a_1 , "Par_U"] , 0)) - 
      (Pr.frame$cum.Pr[i] - min(Pr.frame[Pr.frame$Y[i] > Pr.frame$Y & Pr.frame$X[i] == Pr.frame$X & Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$a_0[i] == Pr.frame$a_0 & Pr.frame$a_1[i] == Pr.frame$a_1 , "cum.Pr"] , 0)) >= 0 | Pr.frame$Y[i] == 1 , 
    TRUE , FALSE)) ,
                        TRUE , FALSE)
}

condition.3 <- as.logical(output.3 * output.4)

# Combine conditions into one identification condition.

Identification <- as.logical(condition.1 * condition.2 * condition.3)

# Attach Identification results to the data.frame.

Pr.frame$Identification <- Identification

# We now need to combine the identification conditions for each individual parameter into a vector that can be checked.

Results <- expand.grid(a_0 , a_1)
colnames(Results) <- c("a_0" , "a_1")
output.final <- vector(length = nrow(Results))
for(i in seq(1,nrow(Results),1)){
  output.final[i] <- ifelse(
    all(
      Pr.frame[Pr.frame$a_0 == Results$a_0[i] & Pr.frame$a_1 == Results$a_1[i] , "Identification"]) == TRUE ,
    TRUE , FALSE)
}

# Final results can be compiled in the following data.frame.

Results$Identified <- output.final

# Restrict attention to the identified set.

library(ggplot2)

Graphical <- subset(Results , Identified == TRUE , select = c(a_0 , a_1))
p <- ggplot(data = Graphical , aes(x = a_0 , y = a_1 ))
p + geom_point(alpha = 0.3)