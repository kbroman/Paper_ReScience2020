######################################################################
# Final Results
######################################################################

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


######################################################################
# Summary information for Gyanu
######################################################################
# Overall summary
mean(finalres$bygene) # 35%
quantile(finalres$total/4204, c(0.025, 0.975)) # 28 - 41%

# genes with >= 75% prob'y of being essential
wh <- cbind(numTAs,finalres$bygene)[finalres$bygene >= 0.749,]
wh <- wh[rev(order(wh[,2])),]
wh[,2] <- round(wh[,2],3)
colnames(wh)[2] <- "Prob.essential"

# families with high or low prob'y to be enriched in ess'l genes
whfam <- (1:109)[!is.na(famprob) & (famprob <= 10 | famprob >= 75)]
famme <- tapply(finalres$bygene, geneclasses, mean,na.rm=TRUE)[as.character(whfam)]
famci <- t(apply(finalres$byfam[,whfam],2,quantile,c(0.025,0.975)))/ngeneperfam[whfam]
fams <- cbind(whfam,famprob[whfam],famme*100,famci*100)
fams <- round(fams)
fams <- cbind(fams,ngeneperfam[fams[,1]])
colnames(fams) <- c("family","prob.enriched","percent.essential","2.5%","97.5%","n.genes")
fams <- fams[fams[,6] > 3,]
fams <- fams[rev(order(fams[,2])),]
