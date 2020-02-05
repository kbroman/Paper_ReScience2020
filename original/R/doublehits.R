######################################################################
# Is it a surprise to see as many as 18 TA sites observed twice?

# Assume no genes (all 65649 sites) are non-essential; sample 986 mutants
n.sim <- 10000
ge2 <- ge3 <- 1:n.sim
for(i in 1:n.sim) {
  x <- table(sample(1:65649,986,repl=TRUE))
  ge2[i] <- sum(x>1) # mean = 7.3    sd = 2.6    Pr(>=16) = 3/1000
  ge3[i] <- sum(x>2) # mean = 0.04
}

# Suppose instead only 36763 sites (=56%)
ge2b <- ge3b <- 1:n.sim
for(i in 1:n.sim) {
  x <- table(sample(1:36763,986,repl=TRUE))
  ge2b[i] <- sum(x>1) # ave = 13     sd = 3.5   Pr(>=16) = 23%
  ge3b[i] <- sum(x>2) # ave =  0.1
}


######################################################################
# do it all again with 1189 mutants

# Assume no genes (all 65649 sites) are non-essential; sample 1189 mutants
n.sim <- 10000
ge2 <- ge3 <- 1:n.sim
for(i in 1:n.sim) {
  x <- table(sample(1:65649,1189,repl=TRUE))
  ge2[i] <- sum(x>1) # mean = 10.6    sd = 3.2    Pr(>=16) = 7%
  ge3[i] <- sum(x>2) # mean =  0.06
}

# Suppose instead only 36763 sites (=56%)
ge2b <- ge3b <- 1:n.sim
for(i in 1:n.sim) {
  x <- table(sample(1:36763,1189,repl=TRUE))
  ge2b[i] <- sum(x>1) # ave = 18.7     sd = 4.2   Pr(>=16) = 77%;  Pr(<=16) = 30%
  ge3b[i] <- sum(x>2) # ave =  0.2
}


rm(i,x,ge2,ge3,ge2b,ge3b,n.sim)

