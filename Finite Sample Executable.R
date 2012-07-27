rm(list=ls(all=TRUE))

library("mvtnorm")

# Enter parameters of the model.

M <- 4
C <- c(-Inf,-0.1,0,0.1,Inf)

X <- c(-1,1)
G <- c(-Inf,0,Inf)

Z <- c(-1,-7/9,-5/9,-3/9,-1/9,1/9,3/9,5/9,7/9,1)

alpha.1 <- 1
alpha.2 <- 0.5

Mu <- c(0,0)
Sigma <- matrix(c(1,0.5,0.5,1),nrow=2,ncol=2)

# Enter range of values over which you wish to search for the identified set.

a.0 <- 0
a.1 <- seq(0.2 , 1.6 , by = 0.01)

# Enter number of observations and number of times to repeat.

n <- c(100 , 200 , 300 , 400 , 500 , 1000 , 2000 , 3000 , 5000 , 8000 , 13000 , 22000 , 36000 , 60000 , 100000 , 160000 , 270000 , 440000 , 700000 , 1000000)
r <- 25

if(any(n%%length(Z) != 0))  stop("n is not a multiple of Z") 

# The following code will create a data.frame of all possible combinations of (a_0 , a_1 , Y , X , Z). 

Y <- seq(1 , M , 1)
Pr.frame <- expand.grid(Y , X , Z , a.1 )
colnames(Pr.frame) <- c("Y" , "X" , "Z" , "a.1")

# Compute parameter values.

Par.L <- vector(length = nrow(Pr.frame)) 
Par.U <- vector(length = nrow(Pr.frame))
for(i in seq(1,nrow(Pr.frame),1)){
  Par.L[i] <- C[match(Pr.frame$Y[i] , Y)] - a.0 - Pr.frame$a.1[i]*Pr.frame$X[i]
  Par.U[i] <- C[match(Pr.frame$Y[i] , Y) + 1] - a.0 - Pr.frame$a.1[i]*Pr.frame$X[i]
}

# Combine these parameters with the data.frame.

Pr.frame$Par.L <- Par.L
Pr.frame$Par.U <- Par.U
rm(list=c("Par.L" , "Par.U"))

# Normalise parameters such that they fall in the [0,1] interval.

Pr.frame <- transform(Pr.frame , Par.L=pnorm(Par.L))
Pr.frame <- transform(Pr.frame , Par.U=pnorm(Par.U))

# Run loop.

k <- 1
Store <- numeric(length = r * length(n))
Sample.size <- numeric(length = r * length(n))

for(m in 1:length(n)){
  repeat{
    
    k <- k+1
    
    Sample.data <- data.frame(rmvnorm(n[m] , mean = Mu , sigma = Sigma))
    colnames(Sample.data) <- c("w" , "v")
    
    Sample.data$z <- c(rep(Z , each = n[m]/length(Z)))
    
    Sample.data$x <- cut(alpha.2 * Sample.data$z + Sample.data$v , breaks = G , labels = X)
    Sample.data <- transform(Sample.data , x = as.numeric(x))
    Sample.data <- transform(Sample.data , x = X[x])
    
    Sample.data$y <- cut(alpha.1 * Sample.data$x + Sample.data$w , breaks = C , labels = c(1:M))
    Sample.data <- transform(Sample.data , y = as.numeric(y))
    
    Sample.space <- data.frame(expand.grid(Y , X , Z))
    colnames(Sample.space) <- c("Y" , "X" , "Z")
    
    Pr <- vector(length = nrow(Sample.space))	
    for(i in seq(1 , nrow(Sample.space) , 1)){
      Pr[i] <- (nrow(subset(Sample.data , Sample.space$Y[i] == Sample.data$y & Sample.space$X[i] == Sample.data$x & Sample.space$Z[i] == Sample.data$z)))/
        (nrow(subset(Sample.data , Sample.space$Z[i] == Sample.data$z)))
    }
    
    Sample.space$Pr <- Pr
    rm(Pr)
    
    Pr <- vector(length = nrow(Pr.frame))
    for(i in seq(1 , nrow(Pr.frame) , 1)){
      Pr[i] <- Sample.space[Pr.frame$Y[i] == Sample.space$Y & Pr.frame$X[i] == Sample.space$X & Pr.frame$Z[i] == Sample.space$Z , "Pr"]
    } 
    
    Pr.frame$Pr <- Pr
    rm(list = c("Sample.data" , "Sample.space"))
    
    cum.Pr <- vector(length = nrow(Pr.frame))
    for(i in seq(1,nrow(Pr.frame),1)){
      cum.Pr[i] <- sum(Pr.frame[Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$X[i] == Pr.frame$X & Pr.frame$Y[i] >= Pr.frame$Y & Pr.frame$a.1[i] == Pr.frame$a.1 , "Pr"])
    }
    
    Pr.frame$cum.Pr <- cum.Pr
    rm(cum.Pr)
    
    output.1 <- vector(length = nrow(Pr.frame))
    for(i in seq(1,nrow(Pr.frame),1)){
      output.1[i] <- sum(Pr.frame[Pr.frame$Par.U[i] >= Pr.frame$Par.U & Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$a.1[i] == Pr.frame$a.1 & Pr.frame$Y != M, "Pr"])
    } 
    condition.1 <- ifelse(output.1 <= Pr.frame$Par.U , TRUE , FALSE)
    rm(output.1)
    
    output.2 <- vector(length = nrow(Pr.frame))
    for(i in seq(1,nrow(Pr.frame),1)){
      output.2[i] <- sum(Pr.frame[Pr.frame$Par.U[i] > Pr.frame$Par.L & Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$a.1[i] == Pr.frame$a.1 , "Pr"])
    } 
    condition.2 <- ifelse(output.2 >= Pr.frame$Par.U | Pr.frame$Y == M , TRUE , FALSE)
    rm(output.2)
    
    output.3 <- vector(length = nrow(Pr.frame))
    for(i in seq(1,nrow(Pr.frame),1)){
      output.3[i] <- ifelse(Pr.frame$Par.U[i] - Pr.frame$cum.Pr[i] >= 0, TRUE , FALSE)
    }
    
    output.4 <- vector(length = nrow(Pr.frame))
    for(i in seq(1,nrow(Pr.frame),1)){
      output.4[i] <- ifelse(any(ifelse(
        (Pr.frame$Par.U[i] - min(Pr.frame[Pr.frame$Y[i] > Pr.frame$Y & Pr.frame$X[i] == Pr.frame$X & Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$a.1[i] == Pr.frame$a.1 , "Par.U"] , 0)) - 
          (Pr.frame$cum.Pr[i] - min(Pr.frame[Pr.frame$Y[i] > Pr.frame$Y & Pr.frame$X[i] == Pr.frame$X & Pr.frame$Z[i] == Pr.frame$Z & Pr.frame$a.1[i] == Pr.frame$a.1 , "cum.Pr"] , 0)) >= 0 | Pr.frame$Y[i] == 1 , 
        TRUE , FALSE)) ,
                            TRUE , FALSE)
    }
    
    condition.3 <- as.logical(output.3 * output.4)
    rm(list=c("output.3" , "output.4"))
    
    Identification <- as.logical(condition.1 * condition.2 * condition.3)
    rm(list=c("condition.1" , "condition.2" , "condition.3"))
    
    Pr.frame$Identification <- Identification
    
    Results <- expand.grid(a.1)
    colnames(Results) <- c("a.1")
    output.final <- vector(length = nrow(Results))
    for(i in seq(1,nrow(Results),1)){
      output.final[i] <- ifelse(
        all(
          Pr.frame[Pr.frame$a.1 == Results$a.1[i] , "Identification"]) == TRUE ,
        TRUE , FALSE)
    }
    
    Results$Identified <- output.final
    rm(output.final)
    
    Store[k-1] <- subset(Results, Identified == TRUE , select = a.1)
    Sample.size[k-1] <- n[m]
    
    print(k-1)
    
    if(k==r*m+ 1)
      break 
  }
}

Finite.identification <- matrix(nrow = length(Store) , ncol = 3)
for(i in 1:length(Store)){
  Finite.identification[i,1] <- ifelse(length(Store[[i]])==0, NA, max(Store[[i]]))
  Finite.identification[i,2] <- ifelse(length(Store[[i]])==0, NA, min(Store[[i]]))	
}
Finite.identification[,3] <- Finite.identification[,1] - Finite.identification[,2]

Results.frame <- data.frame(Finite.identification , Sample.size)
colnames(Results.frame) <- c("Maximum" , "Minimum" , "Range" , "Observations")

Results.frame <- Results.frame[!is.na(Results.frame$Range) , ]

Summary.stats <- data.frame(n)

Mean.range <- vector(length = nrow(Summary.stats))
Maximum.range <- vector(length = nrow(Summary.stats))
Maximum.maximum <- vector(length = nrow(Summary.stats))
Minimum.minimum <- vector(length = nrow(Summary.stats))
Simulations <- vector(length = nrow(Summary.stats))
for(i in 1:nrow(Summary.stats)){
  Mean.range[i] <- mean(Results.frame[Results.frame$Observations == Summary.stats$n[i] , "Range"])
  Maximum.range[i] <- max(Results.frame[Results.frame$Observations == Summary.stats$n[i] , "Range"])
  Maximum.maximum[i] <- max(Results.frame[Results.frame$Observations == Summary.stats$n[i] , "Maximum"])
  Minimum.minimum[i] <- min(Results.frame[Results.frame$Observations == Summary.stats$n[i] , "Minimum"])
  Simulations[i] <- nrow(Results.frame[Results.frame$Observations == Summary.stats$n[i] , ])
}

Summary.stats$Supremum <- Maximum.maximum
Summary.stats$Infimum <- Minimum.minimum
Summary.stats$Mean.range <- Mean.range
Summary.stats$Observed.range <- Maximum.range
Summary.stats$Covered.range <- Summary.stats$Supremum - Summary.stats$Infimum
Summary.stats$Simulatons <- Simulations
rm(list = c("Mean.range" , "Maximum.range" , "Maximum.maximum" , "Minimum.minimum" , "Simulations"))