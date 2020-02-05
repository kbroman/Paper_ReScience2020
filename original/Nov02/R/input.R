temp <- negenes(mydata[,1], mydata[,2], mydata[,3], mydata[,4],
                n.mcmc=50000, skip=49, return=TRUE, trace=FALSE)

mymat <- matrix(ncol=109,nrow=nrow(temp$output))
for(i in 1:109) 
  if(any(geneclasses==i)) 
    mymat[,i] <- apply(1-temp$output[,geneclasses==i,drop=FALSE],1,sum)

finalres <- list(total = temp$n.essential,
                 bygene = temp$geneprob,
                 byfam = mymat)


for(j in 1:9) {
  cat(j,"\n")
  temp <- negenes(mydata[,1], mydata[,2], mydata[,3], mydata[,4],
                  n.mcmc=50000, skip=49, return=TRUE, trace=FALSE)

  mymat <- matrix(ncol=109,nrow=nrow(temp$output))
  for(i in 1:109) 
    if(any(geneclasses==i)) 
      mymat[,i] <- apply(1-temp$output[,geneclasses==i,drop=FALSE],1,sum)

  finalres <- list(total = c(finalres$total, temp$n.essential),
                   bygene = finalres$bygene + temp$geneprob,
                   byfam = rbind(finalres$byfam, mymat))
}

rm(mymat,temp)
famprob <- 1:109
for(i in 1:109)
  famprob[i] <- mean(finalres$byfam[,i]/ngeneperfam[i] > finalres$total/4207)
famprob <- famprob*100

finalres$bygene <- finalres$bygene/10
