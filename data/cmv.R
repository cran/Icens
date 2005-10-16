 bet<- matrix(scan("CMVdata", quiet=TRUE),nc=5,byr=TRUE)

 cmv1 <- rep(bet[,1],bet[,5])
 cmv2 <- rep(bet[,2],bet[,5])

 mac1<-rep(bet[,3],bet[,5])
 mac2<-rep(bet[,4],bet[,5])

 cmv <- cbind(cmv1,cmv2, mac1, mac2)
 rm(bet, cmv1, cmv2, mac1, mac2) 

 dimnames(cmv) <- list(NULL, c("cmvL", "cmvR", "macL", "macR"))

