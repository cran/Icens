#some random data
library(Icens)

intvlx <- matrix(c(
0.8820387, 10.666764,
15.2923703, 18.390665,
10.0710104,
18.9,
7.9796946, 10.964210,
5.2703924, 11.267734,
18.7,
19.875977,
5.9667531, 19.886629,
9.7729062, 13.055671,
3.1947369,
7.482414,
4.2636605,  7.216566,
5.3197158, 15.686208,
0.2885009,
11.463272,
0.2885009,
11.463272),ncol=2,byrow=TRUE)

intvly <- matrix(c(
8.431484, 11.324923,
9.6,
18.739108,
1.438516,  3.232738,
10.6, 11.711857,
14.298833,
16.752745,
9.431221, 16.958045,
2.396955,  7.541405,
12.334413,
21.932913,
7.0, 19.268005,
9.342461, 13.843589,
14.717762,
22.361883,
16.983453, 20.541734,
7.918273, 10.),ncol=2,byrow=TRUE)

#find the cliques
BVcliques(intvlx,intvly)

#find the support
BVsupport(intvlx,intvly)

#find the clique matrix
clmat <- BVclmat(BVcliques(intvlx,intvly))

#the matrix is rank deficient
 clmat[4,]+clmat[7,]-clmat[3,]-clmat[8,]

#should be the zero vector

#now for some estimation

 p1 <- VEM(clmat)

 p2 <- PGM(clmat)

 #p3 seems to be different from p1 and p2!
 p3 <- EMICM(clmat)

 # so is the est unique?

  w<-clmat%*%t(clmat)
  b<-eigen(w)

   b$values
   # one zero eigenvalue

   ev1 <- b$vectors[,10]

   #but the estimator is unique since we cannot move in the direction of
   #recesion
