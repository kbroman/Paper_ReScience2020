######################################################################
# Is it a surprise to see as many as 18 TA sites observed twice?

# Assume no genes (all 65649 sites) are non-essential; sample 1183 mutants
# 65649 = sum(data100stop[,1] + data100stop[,3])
n.sim <- 10000
ge2 <- ge3 <- 1:n.sim
for(i in 1:n.sim) {
  x <- table(sample(1:65649,1183,repl=TRUE))
  ge2[i] <- sum(x>1) # mean = 10.5    sd = 3.1    Pr(>=22) = 1/1000
  ge3[i] <- sum(x>2) # mean = 0.065
}

# Suppose instead only 42672 eligible sites (=65%)
ge2b <- ge3b <- 1:n.sim
for(i in 1:n.sim) {
  x <- table(sample(1:42672,1183,repl=TRUE))
  ge2b[i] <- sum(x>1) # ave = 15.9     sd = 3.9   Pr(>=22) = 7%
  ge3b[i] <- sum(x>2) # ave =  0.1
}

# Suppose instead only 39389 sites (=60%)
ge2c <- ge3c <- 1:n.sim
for(i in 1:n.sim) {
  x <- table(sample(1:39389,1183,repl=TRUE))
  ge2c[i] <- sum(x>1) # ave = 17.4     sd = 4.1   Pr(>=22) = 15%
  ge3c[i] <- sum(x>2) # ave =  0.2
}

