
# schematic of Himar1 transposon
postscript("../Figs/himar.ps",horiz=TRUE,height=5.5,width=9)
plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n",
     xlim=c(0,100),ylim=c(30,100))
arrows(15,80,85,80,code=3,lwd=3,length=0.2,col="blue")
segments(c(25,75),c(80,80),c(25,75),c(90,90),lwd=2,col="blue")
arrows(c(25,75),c(90,90),c(15,85),c(90,90),lwd=2,length=0.1,col="blue")
text(40,84,"KanR",cex=1.6,font=2)
text(60,84,"Ori R6K",cex=1.6,font=2)
text(15,73,"Inverted end",cex=1.4,font=2)
text(85,73,"Inverted end",cex=1.4,font=2)
text(15,95,"T7 Promoter",cex=1.4,font=2)
text(85,95,"T7 Promoter",cex=1.4,font=2)
text(50,74,"2.1 kbp",font=2,cex=1.6,col="darkviolet")
segments(c(5,95),c(69,69),c(50,50),c(40,40),lwd=2,col="red")
arrows(50,40,50,30,lwd=2,col="red")
rect(0,69,100,100,lwd=2)
dev.off()



loc <- c(0.39,0.55,0.89,1.94,2.05,2.93,3.05,3.25,
         3.90,4.05,4.19,4.9,5.72,5.89)
loc2 <- c(1.33,1.66,2.44,2.62,5.22,6.22)

plot.genome <-
function(x,y,r=10,LOC=loc,LOC2=loc2)
{
  z <- seq(0,2*pi,len=501)
  lines(r*sin(z)+x,r*cos(z)+y)

  zp <- z[z>0.1 & z<1.2]
  lines(r*sin(zp)+x,r*cos(zp)+y,lwd=5,col="red")
  zp <- z[z>1.8 & z<2.2]
  lines(r*sin(zp)+x,r*cos(zp)+y,lwd=5,col="green")
  zp <- z[z>2.8 & z<3.4]
  lines(r*sin(zp)+x,r*cos(zp)+y,lwd=5,col="green")
  zp <- z[z>3.7 & z<4.4]
  lines(r*sin(zp)+x,r*cos(zp)+y,lwd=5,col="red")
  zp <- z[z>4.6 & z<5.0]
  lines(r*sin(zp)+x,r*cos(zp)+y,lwd=5,col="green")
  zp <- z[z>5.5 & z<6.0]
  lines(r*sin(zp)+x,r*cos(zp)+y,lwd=5,col="red")
  segments((r*1.05)*sin(LOC)+x,(r*1.05)*cos(LOC)+y,
           (r*1.25)*sin(LOC)+x,(r*1.25)*cos(LOC)+y,lwd=2)
  segments(r*sin(LOC2)+x,r*cos(LOC2)+y,
           (r*1.25)*sin(LOC2)+x,(r*1.25)*cos(LOC2)+y,lwd=2,col="orange")
}


postscript("../Figs/mut1.ps",horiz=TRUE,height=6.5,width=10)
par(mar=c(0.1,0.1,0.1,0.1))
plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n",
     xlim=c(0,100),ylim=c(-1,69))
plot.genome(5,34,r=5)

text(0,16,"Red = essential", col="red", cex=1.4, font=2, adj=c(0,0.5))
text(0,12,"Green = non-essential", col="green", cex=1.4, font=2, adj=c(0,0.5))
dev.off()

loc3 <- c(loc,loc2)
wh <- c(6,8,1,16,4,3,9,19,13,8)

postscript("../Figs/mut2.ps",horiz=TRUE,height=6.5,width=10)
par(mar=c(0.1,0.1,0.1,0.1))
plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n",
     xlim=c(0,100),ylim=c(-1,69))
plot.genome(5,34,r=5)

text(0,16,"Red = essential", col="red", cex=1.4, font=2, adj=c(0,0.5))
text(0,12,"Green = non-essential", col="green", cex=1.4, font=2, adj=c(0,0.5))

arrows(15,34,25,34,lwd=3,col="purple3",length=0.2)
text(18,43,"Create random\nsingle-insertion\nmutants",cex=1.4,col="purple3",font=2)
dev.off()


postscript("../Figs/mut3.ps",horiz=TRUE,height=6.5,width=10)
par(mar=c(0.1,0.1,0.1,0.1))
plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n",
     xlim=c(0,100),ylim=c(-1,69))
plot.genome(5,34,r=5)

text(0,16,"Red = essential", col="red", cex=1.4, font=2, adj=c(0,0.5))
text(0,12,"Green = non-essential", col="green", cex=1.4, font=2, adj=c(0,0.5))

arrows(15,34,25,34,lwd=3,col="purple3",length=0.2)
text(18,43,"Create random\nsingle-insertion\nmutants",cex=1.4,col="purple3",font=2)

xloc <- c(rep(35,5),rep(50,5))
yloc <- rep(c(4,19,34,49,64),2)
for(i in 1:10) {
  plot.genome(xloc[i],yloc[i],r=5)
  temp <- loc3[wh[i]];
  arrows(8*sin(temp)+xloc[i],8*cos(temp)+yloc[i],
         6.5*sin(temp)+xloc[i],6.5*cos(temp)+yloc[i],
         col="blue",len=0.1,lwd=2)
}
dev.off()


postscript("../Figs/mut4.ps",horiz=TRUE,height=6.5,width=10)
par(mar=c(0.1,0.1,0.1,0.1))
plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n",
     xlim=c(0,100),ylim=c(-1,69))
plot.genome(5,34,r=5)

text(0,16,"Red = essential", col="red", cex=1.4, font=2, adj=c(0,0.5))
text(0,12,"Green = non-essential", col="green", cex=1.4, font=2, adj=c(0,0.5))

arrows(15,34,25,34,lwd=3,col="purple3",length=0.2)
text(18,43,"Create random\nsingle-insertion\nmutants",cex=1.4,col="purple3",font=2)

xloc <- c(rep(35,5),rep(50,5))
yloc <- rep(c(4,19,34,49,64),2)
for(i in 1:10) {
  plot.genome(xloc[i],yloc[i],r=5)
  temp <- loc3[wh[i]];
  arrows(8*sin(temp)+xloc[i],8*cos(temp)+yloc[i],
         6.5*sin(temp)+xloc[i],6.5*cos(temp)+yloc[i],
         col="blue",len=0.1,lwd=2)
}

arrows(60,34,70,34,lwd=3,col="purple3",length=0.2)
text(65,43,"Plate\nand select\nviable mutants",cex=1.4,col="purple3",font=2)
dev.off()



postscript("../Figs/mut5.ps",horiz=TRUE,height=6.5,width=10)
par(mar=c(0.1,0.1,0.1,0.1))
plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n",
     xlim=c(0,100),ylim=c(-1,69))
plot.genome(5,34,r=5)

text(0,16,"Red = essential", col="red", cex=1.4, font=2, adj=c(0,0.5))
text(0,12,"Green = non-essential", col="green", cex=1.4, font=2, adj=c(0,0.5))

arrows(15,34,25,34,lwd=3,col="purple3",length=0.2)
text(18,43,"Create random\nsingle-insertion\nmutants",cex=1.4,col="purple3",font=2)

xloc <- c(rep(35,5),rep(50,5))
yloc <- rep(c(4,19,34,49,64),2)
for(i in 1:10) {
  plot.genome(xloc[i],yloc[i],r=5)
  temp <- loc3[wh[i]];
  arrows(8*sin(temp)+xloc[i],8*cos(temp)+yloc[i],
         6.5*sin(temp)+xloc[i],6.5*cos(temp)+yloc[i],
         col="blue",len=0.1,lwd=2)
}

arrows(60,34,70,34,lwd=3,col="purple3",length=0.2)
text(65,43,"Plate\nand select\nviable mutants",cex=1.4,col="purple3",font=2)

xloc <- c(rep(80,5),rep(95,5))
yloc <- rep(c(4,19,34,49,64),2)
for(i in c(1,2,4,5,8,10)) {
  plot.genome(xloc[i],yloc[i],r=5)
  temp <- loc3[wh[i]];
  arrows(8*sin(temp)+xloc[i],8*cos(temp)+yloc[i],
         6.5*sin(temp)+xloc[i],6.5*cos(temp)+yloc[i],
         col="blue",len=0.1,lwd=2)
}
dev.off()


postscript("../Figs/mut6.ps",horiz=TRUE,height=6.5,width=10)
par(mar=c(0.1,0.1,0.1,0.1))
plot(0,0,type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n",
     xlim=c(0,100),ylim=c(-1,69))
plot.genome(5,34,r=5)

text(0,16,"Red = essential", col="red", cex=1.4, font=2, adj=c(0,0.5))
text(0,12,"Green = non-essential", col="green", cex=1.4, font=2, adj=c(0,0.5))

arrows(15,34,25,34,lwd=3,col="purple3",length=0.2)
text(18,43,"Create random\nsingle-insertion\nmutants",cex=1.4,col="purple3",font=2)

xloc <- c(rep(35,5),rep(50,5))
yloc <- rep(c(4,19,34,49,64),2)
for(i in 1:10) {
  plot.genome(xloc[i],yloc[i],r=5)
  temp <- loc3[wh[i]];
  arrows(8*sin(temp)+xloc[i],8*cos(temp)+yloc[i],
         6.5*sin(temp)+xloc[i],6.5*cos(temp)+yloc[i],
         col="blue",len=0.1,lwd=2)
}

arrows(60,34,70,34,lwd=3,col="purple3",length=0.2)
text(65,43,"Plate\nand select\nviable mutants",cex=1.4,col="purple3",font=2)

xloc <- c(rep(80,5),rep(95,5))
yloc <- rep(c(4,19,34,49,64),2)
for(i in c(1,2,4,5,8,10)) {
  plot.genome(xloc[i],yloc[i],r=5)
  temp <- loc3[wh[i]];
  arrows(8*sin(temp)+xloc[i],8*cos(temp)+yloc[i],
         6.5*sin(temp)+xloc[i],6.5*cos(temp)+yloc[i],
         col="blue",len=0.1,lwd=2)
}

r <- 5
segments(c(xloc[4]-r,xloc[4]+r),rep(yloc[4]+r,2),
         c(xloc[4]+r,xloc[4]-r),rep(yloc[4]-r,2),lwd=4,col="purple3")
segments(c(xloc[8]-r,xloc[8]+r),rep(yloc[8]+r,2),
         c(xloc[8]+r,xloc[8]-r),rep(yloc[8]-r,2),lwd=4,col="purple3")

dev.off()
