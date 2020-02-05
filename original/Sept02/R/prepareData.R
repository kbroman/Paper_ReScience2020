######################################################################
# Read data
taloc <- cbind(read.csv("../Data/TAloc.csv"),obs=0)
info <- read.csv("../Data/geneinfo.csv")
classes <- levels(info$class)
info$class <- as.numeric(info$class)
mut <- read.csv("../Rawdata/GenomicData-final4_rev.csv")
mut <- mut[,-c(3,5)]
colnames(mut) <- c("ORFsize","POI","MT")
mut[,2] <- mut[,2] + 2

newmut <- read.csv("../Rawdata/Phase2-FinalData_rev.csv")
newmut[,2] <- newmut[,2]+2
colnames(newmut) <- colnames(mut)

mut <- rbind(mut,newmut)
rm(newmut)

######################################################################
# try to find nearest TA site
lendiff <- mindiff <- mindiff2 <- 1:nrow(mut)
wh <- (1:nrow(taloc))
for(i in 1:nrow(mut)) {
  a <- wh[taloc[,1]==mut[i,3]]
  differ <- abs(mut[i,2]-taloc[a,2])
  differ2 <- mut[i,2]-taloc[a,2]
  mindiff[i] <- min(differ)
  mindiff2[i] <- differ2[differ==min(differ)][1]
  lendiff[i] <- abs(mut[i,1]-info[info[,1]==mut[i,3],2])
  taloc[a[differ==min(differ)][1],4] <- 
    taloc[a[differ==min(differ)][1],4] + 1
  cat(i,mindiff[i],mindiff2[i],lendiff[i],"\n")
}
rm(i,a,wh,differ,mindiff,lendiff,differ2,mindiff2)

######################################################################$
# Look at double-counted TA sites
doubleta <- cbind(read.csv("../Data/doubleTA.csv"),hits1=0,hits2=0)
doubleta <- cbind(doubleta[,1:2],len1=0,
                  doubleta[,3:4],len2=0, hits=0)
temp1 <- paste(taloc[,1],taloc[,2],sep=":")
temp2a <- paste(doubleta[,1],doubleta[,2],sep=":")
temp2b <- paste(doubleta[,4],doubleta[,5],sep=":")
doubleta[,3] <- taloc[match(temp2a,temp1),3]
doubleta[,6] <- taloc[match(temp2b,temp1),3]
doubleta[,7] <- taloc[match(temp2a,temp1),4] + taloc[match(temp2b,temp1),4]
rm(temp1,temp2a,temp2b)


######################################################################
# Create data sets for use with the gibbs sampler
#

prepData <- 
function(ratio=0.8, tas=taloc, doubles=doubleta, discard.stop=TRUE)
{
  doubles <- doubles[doubles[,2]/doubles[,3] <= ratio |
                     doubles[,5]/doubles[,6] <= ratio, ]

  temp1 <- paste(tas[,1],tas[,2],sep=":")
  temp2a <- paste(doubles[,1],doubles[,2],sep=":")
  temp2b <- paste(doubles[,4],doubles[,5],sep=":")

  # remove double-counted TAs 
  tas <- tas[-c(match(temp2a,temp1),match(temp2b,temp1)),]

  # remove TAs not in early part of gene or in stop codon
  if(!discard.stop)
    tas <- tas[tas[,2]/tas[,3] <= ratio, ]
  else
    tas <- tas[tas[,2]/tas[,3] <= ratio & tas[,2] < tas[,3]-1, ]

  # use only the double-counted TAs that are within
  #     the first 80% of *both* genes; discard those in stop codon
  if(!discard.stop)
    doubles <- doubles[doubles[,2]/doubles[,3] <= ratio &
                       doubles[,5]/doubles[,6] <= ratio, ]
  else
    doubles <- doubles[doubles[,2]/doubles[,3] <= ratio &
                       doubles[,5]/doubles[,6] <= ratio &
                       doubles[,2] < doubles[,3]-1 &
                       doubles[,5] < doubles[,6]-1, ]

  mydata1 <- cbind(tapply(tas[,4],tas[,1],length),
                   tapply(tas[,4],tas[,1],sum))

  mydata2 <- cbind(tapply(doubles[,7],doubles[,1],length),
                   tapply(doubles[,7],doubles[,1],sum))

  n <- nrow(info)
  mydata <- matrix(0,ncol=4,nrow=n)
  dimnames(mydata) <- list(paste(info[,1]),c("n.sites","counts","n.sites2","counts2"))
  mydata[rownames(mydata1),1:2] <- mydata1
  mydata[rownames(mydata2),3:4] <- mydata2

  o <- (mydata[,1]>0 | mydata[,3]>0 | c(mydata[n,3],mydata[-n,3])>0)
  mydata[o,]
}
  
data100stop <- prepData(1,discard.stop=FALSE)
data100nostop <- prepData(1,discard.stop=TRUE)
data90 <- prepData(0.9)
data80 <- prepData(0.8)
data70 <- prepData(0.7)
data60 <- prepData(0.6)

# number of TA sites for each gene
numTAs <- data80[,1]+data80[,3]+ c(0,data80[-nrow(data80),3])

#c(sum(data100[,1])+sum(data100[,3]), sum(data100[,3]),
#  sum(data100[,2])+sum(data100[,4]), sum(data100[,4]), nrow(data100),
#  sum(data100[,2]>0 | data100[,4]>0 | c(0,data100[-nrow(data100),4])>0))

# gene classes
geneclasses <- info[match(rownames(data80),info[,1]),3]
ngeneperfam <- table(factor(geneclasses,levels=1:109))
