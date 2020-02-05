######################################################################
# Picture of results
postscript("Figs/pr_essential.ps", horiz=TRUE,height=7.5,width=10)
par(las=1)
p <- finalres[[2]]/length(finalres[[1]])
p[p==0] <- p[p==0] + runif(sum(p==0),-0.03, 0.03)
plot(numTAs, p, ylim=c(-0.03,1), xlab="no. TA sites",
     ylab="Pr(gene is essential | data)",type="n")
abline(h=seq(0,1,by=0.1),lty=2)
points(numTAs, p, cex=0.6)
dev.off()
rm(p)

######################################################################
# study of gene families
temp <- matrix(1,nrow=1000,ncol=nrow(data80))
temp[,data80[,2]==0] <- data80b$output
output <- temp
genes <- unique(taloc[taloc[,2]/taloc[,3] <= 0.8,1])
subclasses <- info[match(genes,info[,1]),3]
genefam.prop <- matrix(ncol=111,nrow=1000)
for(i in 1:111) {
  if(any(subclasses==i)) 
    genefam.prop[,i] <-
      apply(output[,subclasses==i,drop=FALSE],1,mean)
  cat(i,"\n")
}
genefam.summary <- cbind(mean = 1 - apply(genefam.prop,2,mean),
                         lo = 1 - apply(genefam.prop,2,quantile,0.975,na.rm=TRUE),
                         hi = 1 - apply(genefam.prop,2,quantile,0.025,na.rm=TRUE),
                         n = table(factor(subclasses,1:111)))
rm(genes,output,temp,subclasses,i)

# estimated proportion of essential genes in each class
postscript("Figs/prop_essential_byclass.ps",horiz=TRUE,height=7.5,width=10)
par(las=1)
plot(1:111,genefam.summary[,1],ylim=c(0,1),type="n",
     xlab="Gene family",ylab="Proportion of essential genes")
abline(h=1-out80b$summary[1]/nrow(data80),lty=2)
segments(1:111,genefam.summary[,2],1:111,genefam.summary[,3])
abline(v=c(0.5,seq(10.5,110.5,by=10)),lty=2,col="gray")
points(1:111,genefam.summary[,1])
x <- c(13,101,55)
segments(x,genefam.summary[x,2],x,genefam.summary[x,3],col="blue")
points(x,genefam.summary[x,1],col="blue")
x <- 56
segments(x,genefam.summary[x,2],x,genefam.summary[x,3],col="green")
points(x,genefam.summary[x,1],col="green")
x <- c(104,15,111)
segments(x,genefam.summary[x,2],x,genefam.summary[x,3],col="red")
points(x,genefam.summary[x,1],col="red")
dev.off()

# Posterior probability that a class has more than average essential genes
overall <- apply(output,1,mean)
probs <- 1:111
for(i in 1:111)
  probs[i] <- mean(genefam.prop[,i] < overall)

postscript("Figs/pr_essential_family.ps",horiz=TRUE,height=7.5,width=10)
par(las=1)
plot(1:111, probs, ylim=c(0,1), xlab="Gene family", type="n",
     ylab="Pr(family has a more than ave no. essential genes | data)")
abline(v=c(0.5,seq(10.5,110.5,by=10)),lty=2,col="gray")
abline(v=c(13,101,55),col="blue")
abline(v=56,col="green")
abline(v=c(104,15,111),col="red")
points(1:111, probs)
dev.off()

