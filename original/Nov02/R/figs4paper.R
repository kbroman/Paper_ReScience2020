postscript("Figs/fig1.ps", horiz=FALSE, height=5,width=5)
#jpeg("Figs/fig1.jpg", height=500, width=500, quality=100, pointsize=14)
par(las=1,mar=c(5.1,6.1,4.1,2.1))
z <- taloc[taloc[,4]>0,]
mtnum <- unique(rownames(data100nostop))
y <- as.numeric(mtnum[match(z[,1],mtnum)])
x <- z[,2]/z[,3]*100
o <- (z[,2] == z[,3]-1); x <- x[!o]; y <- y[!o]
plot(x,y,cex=0.7,xlab="Site of insertion",ylab="",
     main="")
u <- par("usr")
text(u[1]-diff(u[1:2])*0.25,mean(u[3:4]),"Gene number",xpd=TRUE,srt=90)
dev.off()

postscript("Figs/fig2.ps", horiz=FALSE,height=4.5,width=6.5)
#jpeg("Figs/fig2.jpg", height=500, width=750, pointsize=14, quality=100)
par(las=1,mar=c(5.1,5.1,4.1,2.1))
p <- finalres[[2]]
p[p==0] <- p[p==0] + runif(sum(p==0),-0.01, 0.01)
plot(numTAs, p, ylim=c(-0.01,1), xlab="Number of TAs in proximal portion of gene",
     ylab="Probability gene is essential",type="n")
abline(h=seq(0,1,by=0.1),lty=2)
points(numTAs, p, cex=0.6)
dev.off()
rm(p)

#postscript("Figs/fig3.ps", horiz=FALSE,height=4.5,width=6.5)
##jpeg("Figs/fig3.jpg", height=500, width=750, pointsize=14, quality=100) 
#par(las=1,mar=c(5.1,5.1,4.1,2.1))
#attach("TApostprob2.RData")
#x <- c(750,1000,1500,2000,3000,4000,5000,6000)
##TAlim75m <- c(TAlim75m[1:2],25.2,TAlim75m[-(1:2)])
##TAlim90m <- c(TAlim90m[1:2],50.3,TAlim90m[-(1:2)])
#yl <- range(c(TAlim90m,TAlim75m));yl[1] <- 0
#plot(x,TAlim90m,ylim=yl,xlab="Number of intragenic mutants",
#     ylab="Minimum number of TAs in gene", 
#     xlim=c(500,6000),type="n")
#abline(h=seq(0,90,by=10),lty=2)
#points(x,TAlim90m)
#points(x,TAlim75m,pch=0)
#rm(TAlim75m,TAlim90m)
#detach(2)
#dev.off()
#rm(x,y,yl,u,z,o)
