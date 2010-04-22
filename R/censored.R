# S-plus/R functions to determine the NPMLE of the distribution for interval-
# censored event time.
# Copyright 1998 Alain Vandal and Robert Gentleman
#
# These functions should work, but their intent is to illustrate the
# concepts involved.  Functions are provided as is, with no guarantee.
# Redistribute freely, without modification & without charge.
# Send comments, criticism & flak to vandal@stat.auckland.ac.nz or
# rgentlem@hsph.harvard.edu# Returns TRUE if s1 is a subset of s2, F otherwise
# These functions are aimed at solving the bivariate/multivariate problem


Subset<-function(s1,s2) {
if (is.null(s1))
	return(TRUE)
if (length(unique(s1))>length(unique(s2)))
	return(FALSE)
else
	for (i in 1:length(s1))
		if (sum(s2==s1[i])==0) return(FALSE)
	return(TRUE)
}

# Returns the cliques for box data
#arguments: intvlx- can be either a matrix or the macs for the 1st component
#	    intvly- can be either a 2 x n matrix of interval endpoints
#	 	    or the macs for the second coordinate
BVcliques<-function(intvlx, intvly, Lxopen=TRUE, Rxopen=FALSE,
                    Lyopen=TRUE, Ryopen=FALSE )
{
    if( is.matrix(intvlx) || is.data.frame(intvlx) )
        intvlx<-Maclist(intvlx, Lopen=Lxopen, Ropen=Rxopen)
    if( is.matrix(intvly) || is.data.frame(intvlx) )
        intvly<-Maclist(intvly, Lopen=Lyopen, Ropen=Ryopen)
    lencand<-0
    cliques<-vector("list",length=0)
    for (i in 1:length(intvlx))
        for (j in 1:length(intvly)) {
            curcand<-Intersection(intvlx[[i]],intvly[[j]])
            if (is.null(curcand)) next
            found<-FALSE
            superset<-FALSE
            k<-1
            while (k<=lencand) {
                if (Subset(curcand,cliques[[k]])) {
                    found<-TRUE
                    break
                }
                if (Subset(cliques[[k]],curcand)) {
                    found<-TRUE
                    superset<-TRUE
                    cliques[[k]]<-curcand
                    k<-k+1
                    break
                }
                k<-k+1
            }
            while(superset && k<=lencand) {
                if (Subset(cliques[[k]],curcand)) {
                    k1<-k
                    while (k1<=(lencand-1)) {
                        cliques[[k1]]<-cliques[[k1+1]]
                        k1<-k1+1
                    }
                    cliques[[lencand]]<-NULL
                    lencand<-lencand-1
                } else
                k<-k+1
            }
            if (!found) {
                cliques<-c(cliques,list(curcand))
                lencand<-lencand+1
            }
        }
    cliques
}

# Returns the box coordinates (xlow,xhi,ylo,yhi) of the clique intersections
# in the plane
BVsupport<-function(intvlx,intvly,cliques=BVcliques(intvlx,intvly)) {
	m<-length(cliques)
	boxes<-matrix(0,nrow=m,ncol=4)
	boxes[,c(2,4)]<-Inf
	for (i in 1:m)
		for (j in 1:length(cliques[[i]])) {
			boxes[i,1]<-max(boxes[i,1],intvlx[cliques[[i]][j],1])
			boxes[i,2]<-min(boxes[i,2],intvlx[cliques[[i]][j],2])
			boxes[i,3]<-max(boxes[i,3],intvly[cliques[[i]][j],1])
			boxes[i,4]<-min(boxes[i,4],intvly[cliques[[i]][j],2])
		}
	dimnames(boxes)<-list(1:m,c("xlo","xhi","ylo","yhi"))
	boxes
}

# Plots the boxes given to matrices of intervals; optionally (textp=TRUE) number
# the boxes in their upper left corner & optionally (showsupp=TRUE) displays the
# support boxes
Plotboxes<-
function(int1,int2,textp=FALSE,showmac=FALSE,showsupp=FALSE,showmp=FALSE,cliques=NULL,
    macprod=NULL,density=c(2,8,20),col=c(2,3,4),offsetx=0.02,offsety=0.03)
{
	plot(c(0, max(int1)), c(0, max(int2)), type = "n",xlab="",ylab="")
	segments(int1[, 1], int2[, 1], int1[, 2], int2[, 1])
	segments(int1[, 2], int2[, 1], int1[, 2], int2[, 2])
	segments(int1[, 2], int2[, 2], int1[, 1], int2[, 2])
	segments(int1[, 1], int2[, 2], int1[, 1], int2[, 1])
	if (is.null(density))
		density<-c(2,8,10)
	else if (length(density)==1)
		density<-rep(density,3)
	else
		density<-rep(density,2)[1:3]
	if (is.null(col))
		col<-c(2,3,4)
	else if (length(col)==1)
		density<-rep(col,3)
	else
		density<-rep(density,2)[1:3]

	if (textp) {
		delx<-offsetx*max(int1)
		dely<-offsety*max(int2)
		text(int1[,1]+delx,int2[,2]-dely,paste(1:dim(int1)[1]))
	}
	if (showmac) {
		temp<-MLEintvl(int1)
		for (i in 1:dim(temp)[1])
			segments(temp[i,1],0,temp[i,2],0,lwd=density[1],col=col[1])
		temp<-MLEintvl(int2)
		for (i in 1:dim(temp)[1])
			segments(0,temp[i,1],0,temp[i,2],lwd=density[1],col=col[1])
	}
	if (showsupp) {
		if (is.null(cliques)) cliques<-BVcliques(int1,int2)
		temp<-BVsupport(int1,int2,cliques=cliques)
		for (i in 1:dim(temp)[1])
            polygon(temp[i,c(1,2,2,1)],temp[i,c(3,3,4,4)],density=density[2],
                    col=col[2])
	}
	if (showmp) {
		if (is.null(macprod)) macprod<-BVmacprod(int1,int2)$mpcoor
		for (i in 1:dim(macprod)[1])
            polygon(macprod[i,c(1,2,2,1)],macprod[i,c(3,3,4,4)],
                    density=density[3],col=col[3],angle=-45)
	}
}

# Returns the clique matrix for the data
BVclmat<-function(cliques) {
	m<-length(cliques)
	temp<-NULL
	for (i in 1:m) temp<-c(temp,cliques[[i]])
	n<-length(unique(temp))
	ret<-matrix(0,m,n)
	for (i in 1:m)
		for (j in 1:length(cliques[[i]]))
			ret[i,cliques[[i]][j]]<-1
	ret
}

# Returns the list containing the product of the maximal antichains as it
# intersects the data
BVmacprod<-function(intvlx,intvly) {
	biclq<-list(Maclist(intvlx),Maclist(intvly))
	m<-c(length(biclq[[1]]),length(biclq[[2]]))
	hmx<-MLEintvl(intvlx,biclq[[1]])
	hmy<-MLEintvl(intvly,biclq[[2]])
	macprod<-vector("list",length=m[1])
	for (i2 in 1:m[2]) macprod[[i2]]<-vector("list",length=m[2])
	mpcoor<-NULL
	for (i1 in 1:m[1])
		for (i2 in 1:m[2]) {
			temp<-Intersection(biclq[[1]][[i1]],biclq[[2]][[i2]])
			if (!is.null(temp)) mpcoor<-rbind(mpcoor,c(hmx[i1,],hmy[i2,]))
			macprod[[i1]][[i2]]<-temp
		}
	dimnames(mpcoor)<-list(NULL,c("xlo","xhi","ylo","yhi"))
	ret<-NULL
	ret$macprod<-macprod
	ret$mpcoor<-mpcoor
	ret
}

Intersection<-function(vec1,vec2) {
	rvec<-NULL
	for(i in vec1)
		for(j in vec2)
			if(i==j) rvec<-c(j,rvec)
	return(rvec)
}


