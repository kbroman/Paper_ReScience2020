
######################################################################
# run gibbs samplers
out100 <- negenes(data100[,1], data100[,2], data100[,3], data100[,4],
                  n.mcmc=100000,skip=49)
out90 <- negenes(data90[,1], data90[,2], data90[,3], data90[,4],
                 n.mcmc=100000,skip=49)
out80 <- negenes(data80[,1], data80[,2], data80[,3], data80[,4],
                 n.mcmc=100000,skip=49)
out70 <- negenes(data70[,1], data70[,2], data70[,3], data70[,4],
                 n.mcmc=100000,skip=49)
out60 <- negenes(data60[,1], data60[,2], data60[,3], data60[,4],
                 n.mcmc=100000,skip=49)

######################################################################
# Final Results

# 3 minutes, 50 seconds
temp <- negenes(data80[,1], data80[,2], data80[,3], data80[,4],
                n.mcmc=50000, skip=49, return=TRUE, trace=FALSE)

mymat <- matrix(ncol=109,nrow=nrow(temp$output))
for(i in 1:109) 
  if(any(geneclasses==i)) 
    mymat[,i] <- apply(1-temp$output[,geneclasses==i,drop=FALSE],1,sum)

finalres <- list(total = temp$n.essential,
                 bygene = apply(1-temp$output,2,sum),
                 byfam = mymat)


# 3 minutes, 50 seconds
for(j in 1:3) {
cat(j,"\n")
temp <- negenes(data80[,1], data80[,2], data80[,3], data80[,4],
                n.mcmc=50000, skip=49, return=TRUE, trace=FALSE)

mymat <- matrix(ncol=109,nrow=nrow(temp$output))
for(i in 1:109) 
  if(any(geneclasses==i)) 
    mymat[,i] <- apply(1-temp$output[,geneclasses==i,drop=FALSE],1,sum)

finalres <- list(total = c(finalres$total, temp$n.essential),
                 bygene = finalres$bygene + apply(1-temp$output,2,sum),
                 byfam = rbind(finalres$byfam, mymat))
}

rm(mymat,temp)
famprob <- 1:109
for(i in 1:109)
  famprob[i] <- mean(finalres$byfam[,i]/ngeneperfam[i] > finalres$total/4207)
famprob <- famprob*100

######################################################################
# test
#n.sites <- rpois(4200,5)+1
#viable <- rep(0,4200)
#viable[sample(1:4200,2500)] <- 1
#dat <- sim.mutants(n.sites,viable,1200)
#output <- negenes(n.sites,dat)
#rm(dat,n.sites,viable,output)

# Genes with high posterior probability of being essential
temp <- cbind(info[match(rownames(data80),info[,1]),],
              numta=numTAs,prob.ess=finalres[[2]]/100)
