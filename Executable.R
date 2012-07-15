library("mvtnorm")

# Enter parameters of the model.

M <- 6
C <- c(-Inf,-0.3,-0.1,0,0.1,0.3,Inf)

X <- c(-1,1)
G <- c(-Inf,0,Inf)

Z <- c(-1,-7/9,-5/9,-3/9,-1/9,1/9,3/9,5/9,7/9,1)

alpha_1 <- 1
alpha_2 <- 0.5

Mu <- c(0,0)
Sigma <- matrix(c(1,0.5,0.5,1),nrow=2,ncol=2)

# Enter grid of values over which you wish to search for the identified set.

a_0 <- seq(0,2,0.5)
a_1 <- seq(0,2,0.5)

# The following code will create a data.frame of all possible combinations of (a_0 , a_1 , Y , X , Z). 

Y <- seq(1 , M , 1)
frame <- expand.grid(a_0 , a_1 , Y , X , Z)
colnames(frame) <- c("a_0" , "a_1" , "Y" , "X" , "Z")

# We now need to calculate the probability of the event (Y=m and X=x | Z=z) using the mvtnorm package. 

Pr <- NULL
for(i in seq(1,nrow(frame),1)){
  Pr <- c(Pr , pmvnorm(lower = c(C[match(frame$Y[i] , Y)] - alpha_1*frame$X[i] , G[match(frame$X[i] , X)] - 0.5*frame$Z[i]) , upper = c(C[match(frame$Y[i] , Y) + 1] - alpha_1*frame$X[i] , G[match(frame$X[i] , X) + 1] - 0.5*frame$Z[i]) , mean = Mu , sigma = Sigma)) 
}

# Combine these probabilities with the data.frame.

Pr.frame <- data.frame(cbind(frame , Pr))
colnames(Pr.frame) <- c("a_0" , "a_1" , "Y" , "X" , "Z" , "Pr")

# Cumulate these probabilities for Y<=m.

cum.Pr <- NULL
for(i in seq(1,nrow(frame),1)){
  cum.Pr <- c(cum.Pr , with(Pr.frame , sum(Pr.frame[Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$X[i] == Pr.frame$X & Pr.frame$Y[i] >= Pr.frame$Y & Pr.frame$a_1[i] == Pr.frame$a_1 & Pr.frame$a_0[i] == Pr.frame$a_0 , "Pr"]))) 
}

# Combine these probabilities with the data.frame.

frame <- data.frame(cbind(Pr.frame , cum.Pr))
colnames(frame) <- c("a_0" , "a_1" , "Y" , "X" , "Z" , "Pr" , "cum.Pr")

# Compute parameter values.

Par_L <- NULL
Par_U <- NULL
for(i in seq(1,nrow(frame),1)){
  Par_L <- c(Par_L , C[match(frame$Y[i] , Y)] - frame$a_0[i] - frame$a_1[i]*frame$X[i])
  Par_U <- c(Par_U , C[match(frame$Y[i] , Y) + 1] - frame$a_0[i] - frame$a_1[i]*frame$X[i])
}

# Combine these parameters with the data.frame.

Pr.frame <- data.frame(cbind(frame , Par_L , Par_U))
colnames(Pr.frame) <- c("a_0" , "a_1" , "Y" , "X" , "Z" , "Pr" , "cum.Pr" , "Par_L" , "Par_U")

# Normalise parameters such that they fall in the [0,1] interval.

Pr.frame <- transform(Pr.frame , Par_L=pnorm(Par_L))
Pr.frame <- transform(Pr.frame , Par_U=pnorm(Par_U))

# The first condition for a combination (a_0 , a_1) lying inside the identified set to be tested.

output.1 <- NULL
for(i in seq(1,nrow(frame),1)){
  output.1 <- c(output.1 , with(Pr.frame , sum(Pr.frame[Pr.frame$Par_L[i] >= Pr.frame$Par_L & Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$a_0[i] == Pr.frame$a_0 & Pr.frame$a_1[i] == Pr.frame$a_1 & Pr.frame$Y != M, "Pr"])))
} 
condition.1 <- ifelse(output.1 <= Pr.frame$Par_L , 1 , 0)

# The second condition for a combination (a_0 , a_1) lying inside the identified set to be tested.

output.2 <- NULL
for(i in seq(1,nrow(frame),1)){
  output.2 <- c(output.1 , with(Pr.frame , sum(Pr.frame[Pr.frame$Par_L[i] > Pr.frame$Par_L & Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$a_0[i] == Pr.frame$a_0 & Pr.frame$a_1[i] == Pr.frame$a_1 , "Pr"])))
} 
condition.2 <- ifelse(output.1 >= Pr.frame$Par_L , 1 , 0)

# The third condition for a combination (a_0 , a_1) lying inside the identified set to be tested.
