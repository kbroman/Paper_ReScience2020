######################################################################
# Results with different cutoffs
######################################################################
res100nostop <- negenes(data100nostop[,1], data100nostop[,2],
                        data100nostop[,3], data100nostop[,4],
                        n.mcmc=500000, skip=49, return=FALSE,
                        trace=FALSE)
res100stop <- negenes(data100stop[,1], data100stop[,2],
                      data100stop[,3], data100stop[,4],
                      n.mcmc=500000, skip=49, return=FALSE,
                      trace=FALSE)
res90 <- negenes(data90[,1], data90[,2],
                 data90[,3], data90[,4],
                 n.mcmc=500000, skip=49, return=FALSE,
                 trace=FALSE)
res70 <- negenes(data70[,1], data70[,2],
                 data70[,3], data70[,4],
                 n.mcmc=500000, skip=49, return=FALSE,
                 trace=FALSE)
res60 <- negenes(data60[,1], data60[,2],
                 data60[,3], data60[,4],
                 n.mcmc=500000, skip=49, return=FALSE,
                 trace=FALSE)


######################################################################
# Final Results
######################################################################

temp <- negenes(data80[,1], data80[,2], data80[,3], data80[,4],
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
  temp <- negenes(data80[,1], data80[,2], data80[,3], data80[,4],
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
mean(finalres$bygene) # 40%
quantile(finalres$total/4204, c(0.025, 0.975)) # 33 - 46%

# Plot of Pr(essential) vs no. TAs
postscript("prob_essential.ps", horiz=FALSE, height=4.5,width=6.5)
par(las=1,mar=c(5.1,6.1,4.1,2.1))
p <- finalres$bygene
p[p==0] <- p[p==0] + runif(sum(p==0),-0.01,0.01)
plot(numTAs, p, ylim=c(-0.01,1), xlab="Number of TAs in proximal 80% of gene",
     ylab="Probability gene is essential",type="n")
abline(h=seq(0,1,by=0.1),lty=2)
points(numTAs, p, cex=0.6)
dev.off()

# genes with >= 75% prob'y of being essential
wh <- cbind(numTAs,finalres$bygene)[numTAs >=50 & finalres$bygene >0,]
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
#fams[rev(order(fams[,2])),]
