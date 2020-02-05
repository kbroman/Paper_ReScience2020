
# create data on locations of insertions around circular genome
#mygeneloc <- read.csv("../Data/mygeneloc.csv")
L <- 4403836 # genome location

# determine locations of all TA insertions in (1..4403836)
a <- taloc[taloc$obs>0,]
insertloc <- 1:nrow(a)
for(i in 1:nrow(a)) {
  g <- a[i,1]
  loc <- a[i,2]
  x <- mygeneloc[mygeneloc[,1]==g,]
  if(x[[4]] == 1) # forward orientation
    insertloc[i] <- x[[2]] + loc - 1
  else
    insertloc[i] <- x[[3]] - loc + 1
}
rm(a,g,loc,x)

postscript("../Figs/circlefig.ps",horiz=FALSE,height=5,width=5)
par(pty="s",bty="n")
plot(0,0,type="l",xlab="",ylab="",xlim=c(-1,1),ylim=c(-1,1),xaxt="n",yaxt="n")
x <- seq(0,L,length=501)
lines(0.9*sin(x/L*2*pi),0.9*cos(x/L*2*pi),lwd=2)
X <- sin(insertloc/L*2*pi)
Y <- cos(insertloc/L*2*pi)
b <- taloc$obs[taloc$obs>0]
b[b==1] <- 0.95
b[b==2] <- 1
segments(0.9*X,0.9*Y,b*X,b*Y)
segments(0.9*sin(0),0.9*cos(0),0.85*sin(0),0.85*cos(0),lwd=2)
text(0.8*sin(0),0.8*cos(0),"0",cex=1.1,adj=c(0.5,1))
dev.off()
