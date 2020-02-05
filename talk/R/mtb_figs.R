postscript("../Figs/mtb_mcmc.ps",horiz=TRUE,height=6.5,width=10)
par(las=1,mar=c(5.1,5.1,0.1,5.1))
plot(0:1000,c(3435,out$n.essential)/42.04,type="l",ylim=c(0,100),
     ylab="Percent essential", xlab="MCMC steps",lwd=2)
text(500,80,"Overall:  39%   (28 - 41%)",adj=c(0,0.5),cex=1.8,font=2,col="blue")
dev.off()


postscript("../Figs/prob_essential.ps", horiz=TRUE,height=6.5,width=10)
par(las=1,mar=c(5.1,5.1,0.1,5.1))
p <- finalres[[2]]
p[p==0] <- p[p==0] + runif(sum(p==0),-0.01, 0.01)
plot(numTAs, p, ylim=c(-0.01,1), xlab="Number of TAs in proximal portion of gene",
     ylab="Probability gene is essential",type="n")
abline(h=seq(0,1,by=0.1),lty=2)
points(numTAs, p, cex=0.6)
dev.off()
rm(p)

# histogram of number of ta sites
postscript("../Figs/numTAs.ps",horiz=TRUE,width=8,height=3)
par(mar=c(4.1,0.1,0.1,0.1))
hist(numTAs,breaks=30,ylab="",yaxt="n",xlab="Number of TA sites per gene",main="",xaxt="n")
u <- par("usr")
x <- numTAs + runif(length(numTAs),-0.4,0.4)
segments(x,u[3]-diff(u[3:4])*0.03,x,0,xpd=TRUE)
segments(0,u[3]-diff(u[3:4])*0.03,u[2],u[3]-diff(u[3:4])*0.03,xpd=TRUE)
segments(0,0,u[2],0,xpd=TRUE)
x <- seq(0,200,by=25)
segments(x,u[3]-diff(u[3:4])*0.03,x,2*(u[3]-diff(u[3:4])*0.03),xpd=TRUE)
text(x,3*(u[3]-diff(u[3:4])*0.03),paste(x),xpd=TRUE)
dev.off()
