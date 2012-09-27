rm(list=ls(all=TRUE))

library("mvtnorm")

# ========================================

# ----------
# Specify the true model.
# ----------
  
# Y = f(X , W) , X = g(Z , V)
# Presume that f is linear and specify the parameter values on X and Z
alpha.1 <- 1       # Coefficient on X
alpha.2 <- 0.5     # Coefficient on Z

# ----------
# Specify a true distribution.
# ----------

# Variance, Covariance matrix.
Variance.W <- 1
Variance.V <- 1    # Always set equal to 1 as a normalisation.
Covariance.W.V <- 0.5

# ----------
# Specify the number of outcomes M and thresholds C.
# ----------

M <- 4
C <- c(-Inf,-0.1,0,0.1,Inf)

# ----------
# Specify the (ordered, increasing) support of X and thresholds G.
# ----------

X <- c(-1,1)
G <- c(-Inf,0,Inf)

# ----------
# Specify the support of Z.
# ----------
  
Z <- c(-1,-7/9,-5/9,-3/9,-1/9,1/9,3/9,5/9,7/9,1)

# ----------
# Specify the number of observations in a finite sample setting.
# ----------

n <- 200            # Note that this must be an integer multiple of length(Z).   

if(any(n%%length(Z) != 0))  stop("n is not a multiple of Z")

# ----------
# Specification is now complete. Run the rest of the code to obtain results
# for maximum likelihood (control function approach) and SEIV estimator
# results. 
# ----------

# ========================================

# ----------
# Generate finite sample data.
# ----------

Sigma <- matrix(c(Variance.W , Covariance.W.V , Covariance.W.V , Variance.V) ,
                nrow = 2 ,
                ncol = 2 ,
                byrow = T)
Mu <- c(0 , 0)
Y <- seq(1 , M , 1)

Generation <- function(n){
  
  Sample.data <- data.frame(rmvnorm(n , mean = Mu , sigma = Sigma))
  colnames(Sample.data) <- c("w" , "v")
  Sample.data$z <- c(rep(Z , each = n/length(Z)))
  Sample.data$x <- cut(alpha.2 * Sample.data$z + Sample.data$v , breaks = G , labels = X)
  Sample.data <- transform(Sample.data , x = as.numeric(x))
  Sample.data <- transform(Sample.data , x = X[x])
  Sample.data$y <- cut(alpha.1 * Sample.data$x + Sample.data$w , breaks = C , labels = c(1:M))
  Sample.data <- transform(Sample.data , y = as.numeric(y))
  
  return(data.frame(Y = Sample.data$y , X = Sample.data$x , Z = Sample.data$z))
}

Sample <- function(t){
  
  Sample.space <- data.frame(expand.grid(Y , X , Z))
  colnames(Sample.space) <- c("Y" , "X" , "Z")
  
  Pr <- vector(length = nrow(Sample.space))  
  for(i in 1:nrow(Sample.space)){
    Pr[i] <- (nrow(subset(t , Sample.space$Y[i] == t$Y & Sample.space$X[i] == t$X & Sample.space$Z[i] == t$Z)))/
      (nrow(subset(t , Sample.space$Z[i] == t$Z)))
  }
  
  Sample.space$Pr <- Pr
  
  return(data.frame(Y = Sample.space$Y , X = Sample.space$X , Z = Sample.space$Z , Pr = Sample.space$Pr))
}

# If the data is such that the SEIV estimator returns false for all parameter
# values in the search area then we want the data to be regenerated.

repeat{
  
Distribution <- Generation(n)

# ----------
# Maximum likelihood estimation for triangular model.
# ----------

Log.likelihood <- function(x){
  a.0 <- x[1]
  a.1 <- x[2]
  b <- x[3]
  Var <- x[4]
  Cov <- x[5]
  output <- vector(length = nrow(Distribution))
  for(i in 1:nrow(Distribution)){
    output[i] <- log(
      pmvnorm(lower = c(C[match(Distribution$Y[i] , Y)] - a.0 - a.1*Distribution$X[i] , G[match(Distribution$X[i] , X)] - b*Distribution$Z[i]) , 
              upper = c(C[match(Distribution$Y[i] , Y) + 1] - a.0 - a.1*Distribution$X[i] , G[match(Distribution$X[i] , X) + 1] - b*Distribution$Z[i]) , 
              mean = c(0 , 0) , 
              sigma = matrix(c(Var , Cov , Cov , 1) , nrow = 2 , ncol = 2 , byrow = T))
                              )
  }
  answer <- sum(output)
  return(-1 * answer)
}
# Optimisation procedure with starting values at true parameter values.
MLE <- optim(c(0 , alpha.1 , alpha.2 , Variance.W , Covariance.W.V) , Log.likelihood , hessian = TRUE)
print(MLE)
# ----------
# Single Equation Instrumental Variable method.
# ----------

# Define conditions for lying in the outer region of the identified set.  
  # Output 1  
    Output.1 <- function(q){
  output.1 <- vector(length = nrow(q))
  for(i in seq(1,nrow(q),1)){
    output.1[i] <- sum(q[q$Par.U[i] >= q$Par.U & q$Z[i] == q$Z & q$Y != M, "Pr"])
  } 
  return(ifelse(output.1 <= q$Par.U , TRUE , FALSE))
  }
  # Output 2
    Output.2 <- function(q){
  output.2 <- vector(length = nrow(q))
  for(i in seq(1,nrow(q),1)){
    output.2[i] <- sum(q[q$Par.U[i] > q$Par.L & q$Z[i] == q$Z , "Pr"])
  } 
  return(ifelse(output.2 >= q$Par.U | q$Y == M , TRUE , FALSE))
  }
  # Output 3
    Output.3 <- function(q){
  output.3 <- vector(length = nrow(q))
  for(i in seq(1,nrow(q),1)){
    output.3[i] <- ifelse(q$Par.U[i] - q$cum.Pr[i] >= 0, TRUE , FALSE)
  }
  output.4 <- vector(length = nrow(q))
  for(i in seq(1,nrow(q),1)){
    output.4[i] <- ifelse(any(ifelse(
      (q$Par.U[i] - min(q[q$Y[i] > q$Y & q$X[i] == q$X & q$Z[i] == q$Z , "Par.U"] , 0)) - 
        (q$cum.Pr[i] - min(q[q$Y[i] > q$Y & q$X[i] == q$X & q$Z[i] == q$Z , "cum.Pr"] , 0)) >= 0 | q$Y[i] == 1 , 
      TRUE , FALSE)) ,
                          TRUE , FALSE)
  }
  
  return(as.logical(output.3 * output.4))
  }
  # Organise probabilities and parameter values into a table.
    Organiser <- function(e){
    q <- Sample(e)
    
    cum.Pr <- vector(length = nrow(q))
    for(i in seq(1,nrow(q),1)){
      cum.Pr[i] <- sum(q[q$Z[i] == q$Z & q$X[i] == q$X & q$Y[i] >= q$Y , "Pr"])
    }
    
    q$cum.Pr <- cum.Pr
    
    Par.L <- vector(length = nrow(q)) 
    Par.U <- vector(length = nrow(q))
    for(i in seq(1,nrow(q),1)){
      Par.L[i] <- C[match(q$Y[i] , Y)] - a.0 - a.1*q$X[i]
      Par.U[i] <- C[match(q$Y[i] , Y) + 1] - a.0 - a.1*q$X[i]
    }
    
    q$Par.L <- Par.L
    q$Par.U <- Par.U
    
    q <- transform(q , Par.L=pnorm(Par.L))
    q <- transform(q , Par.U=pnorm(Par.U))
    
    return(q)
  }
  # Combine the three output conditions into a single logical test.  
    Output <- function(q){
      return(if(all(Output.1(q) * Output.2(q) * Output.3(q)))
      {TRUE}
             else {FALSE})
                }
# Computation using the MLE results as starting values.
  a.0 <- MLE$par[1]
  a.1 <- MLE$par[2]
  Result <- Output(Organiser(Distribution))
# If the MLE results are estimated to be in the identified set then break
# the loop and search for other values in the set. If the MLE results are
# not estimated in the set then resample and begin loop again.
if(Result == TRUE)
  break
}
# The MLE result is estimated in the set. Use polar coordinates to search 
# for other values (a.0,a.1) in the set. First, define a function that 
# validates belonging using polar coordinates (to do so, we will also 
# need to redefine the Organiser function (Pr.function) so that it works for polar
# coordinates).

Nested.engine <- function(r){
  
  h <- MLE$par[1] + r * sin(theta*pi/180)
  v <- MLE$par[2] + r * cos(theta*pi/180)
  
  Workhorse <- Sample(Distribution)
  
    Par.L <- vector(length = nrow(Workhorse)) 
  Par.U <- vector(length = nrow(Workhorse))
  for(i in seq(1,nrow(Workhorse),1)){
    Par.L[i] <- C[match(Workhorse$Y[i] , Y)] - h - v*Workhorse$X[i]
    Par.U[i] <- C[match(Workhorse$Y[i] , Y) + 1] - h - v*Workhorse$X[i]
  }
  
  Workhorse$Par.L <- Par.L
  Workhorse$Par.U <- Par.U
  
  Workhorse <- transform(Workhorse , Par.L=pnorm(Par.L))
  Workhorse <- transform(Workhorse , Par.U=pnorm(Par.U))
  
  cum.Pr <- vector(length = nrow(Workhorse))
  for(i in seq(1,nrow(Workhorse),1)){
    cum.Pr[i] <- sum(Workhorse[Workhorse$Z[i] == Workhorse$Z & Workhorse$X[i] == Workhorse$X & Workhorse$Y[i] >= Workhorse$Y , "Pr"])
  }
  
  Workhorse$cum.Pr <- cum.Pr
  
  return(Workhorse)
}

theta <- 0
k <- 0
Out.r <- vector(length = 360)
Out.theta <- vector(length = 360)
Out.function <- vector(length = 360)
r.increment <- 0.001
theta.increment <- 1
repeat{
  theta <- theta + theta.increment
  r <- 0
  while(Output(Nested.engine(r)) == T){
    k <- k + 1
    Out.r[k] <- r
    Out.function[k] <- Output(Nested.engine(r))
    r <- r + r.increment
    Out.r[k] <- r
    Out.theta[k] <- theta
    print(c(k,theta))
  }
if(theta >= 360)
  break
}

Out <- data.frame(Out.r , Out.theta)
colnames(Out) <- c("r" , "theta")
Out$a.0 <- MLE$par[1] + Out$r * sin(Out$theta * pi / 180)
Out$a.1 <- MLE$par[2] + Out$r * cos(Out$theta * pi / 180) 
Out$Logical <- Out.function
rm(list=c("Out.r" , "Out.theta" , "Out.function"))

library(ggplot2)
p <- ggplot(Out , aes(x = a.0 , y = a.1))
p + geom_point(alpha = 0.5)

Hull.problem <- chull(Out$a.0 , Out$a.1)
Hull.problem <- append(Hull.problem , Hull.problem[1])
Hull <- Out[Hull.problem , c("a.0" , "a.1")]
h <- ggplot(Hull , aes(x = a.0 , y = a.1))
h + geom_polygon(alpha = 0.5) + geom_point(aes(x = MLE$par[1] , y = MLE$par[2]) , colour = "white")