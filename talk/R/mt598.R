mt598 <- scan("mt598.txt",what=character(0))
temp <- vector("list",length(mt598))
for(i in 1:3) {
  temp[[i]] <- rev(strsplit(mt598[i],"")[[1]])
  for(j in 1:length(temp[[i]])) {
    if(temp[[i]][j] == "A") temp[[i]][j] <- "T"
    else if(temp[[i]][j] == "T") temp[[i]][j] <- "A"
    else if(temp[[i]][j] == "G") temp[[i]][j] <- "C"
    else temp[[i]][j] <- "G"
  }
}
mt598[1] <- paste(temp[[3]],collapse="")
mt598[2] <- paste(temp[[2]],collapse="")
mt598[3] <- paste(temp[[1]],collapse="")
rm(temp)

postscript("../Figs/mt598.ps", horiz=FALSE, height=2.5, width=6.5)
par(mar=c(0.1,0.1,0.1,0.1))
plot(0,0,xlim=c(0,100),ylim=c(0,100),type="n",
     xlab="",ylab="",bty="n",xaxt="n",yaxt="n")
u <- par("usr")
first <- u[1]+diff(u[1:2])*0.05
last <- u[2]-diff(u[1:2])*0.05
mult <- (last-first)/51
cury <- u[3]+diff(u[3:4])*5/6
for(i in 1:5) 
  text(first+(i-1)*mult,cury,substr(mt598[1],5+i,5+i),
       cex=0.8,adj=c(0.5,0.5))
arrows(first+4.5*mult,cury+10,first+4.5*mult,cury+3,length=0.1,xpd=TRUE,col="red",lwd=2)
segments(first+c(4.5,4.5,7.5)*mult,cury-c(5,8,5),
         first+c(4.5,7.5,7.5)*mult,cury-c(8,8,8))
text(first+6*mult,cury-14,"start",xpd=TRUE,cex=0.8)
for(i in 1:41) 
  text(first+(i-1+5)*mult,cury,substr(mt598[2],i,i),
       cex=0.8,adj=c(0.5,0.5))
points(first+seq(-2,-1,len=3)*mult,rep(cury,3),xpd=TRUE,cex=0.25,pch=16)
segments(first+(c(10,20,30,40)+4)*mult,cury-5,
         first+(c(10,20,30,40)+4)*mult,cury-8)
text(first+(c(10,20,30,40)+4)*mult,cury-14,as.character(c(10,20,30,40)),cex=0.8)
cury <- u[3]+diff(u[3:4])*0.5
for(i in 1:41)
  text(first+(i-1+5)*mult,cury,substr(mt598[2],i+41,i+41),
       cex=0.8,adj=c(0.5,0.5))
arrows(first+35.5*mult,cury+10,first+35.5*mult,cury+3,length=0.1,col="red",lwd=2)
segments(first+(c(50,60,70,80)-37)*mult,cury-5,
         first+(c(50,60,70,80)-37)*mult,cury-8)
text(first+(c(50,60,70,80)-37)*mult,cury-14,as.character(c(50,60,70,80)),cex=0.8)
cury <- u[3]+diff(u[3:4])/6
for(i in 1:41)
  text(first+((i-1)+5)*mult,cury,substr(mt598[2],i+41*2,i+41*2),
       cex=0.8,adj=c(0.5,0.5))
for(i in 1:5)
  text(first+(i-1+41+5)*mult,cury,substr(mt598[3],i,i),
       cex=0.8,adj=c(0.5,0.5))
arrows(first+43.5*mult,cury+10,first+43.5*mult,cury+3,length=0.1,col="red",lwd=2)
segments(first+c(42.7,42.7,45.5)*mult,cury-c(5,8,5),
         first+c(42.7,45.5,45.5)*mult,cury-c(8,8,8))
text(first+44*mult,cury-14,"stop",xpd=TRUE,cex=0.8)
points(first+seq(51,52,len=3)*mult,rep(cury,3),xpd=TRUE,cex=0.25,pch=16)
segments(first+(c(90,100,110)-78)*mult,cury-5,
         first+(c(90,100,110)-78)*mult,cury-8)
text(first+(c(90,100,110)-78)*mult,cury-14,as.character(c(90,100,110)),cex=0.8)
dev.off()
rm(i,cury,first,last,mult,u,j)
