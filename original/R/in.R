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
for(j in 1:9) {
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

rm(temp,mymat)

