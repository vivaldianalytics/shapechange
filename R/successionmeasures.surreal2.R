library(Matrix)
require(Matrix)
#functions for calculating and plotting rate of succession using surreal numbers

#input data frame data (rows are time points, columns are species abundances) and times (vector of time points)
#return list containing:
#rbar: one row per time point. Columns are coefficients of omega^(1/omega) and omega^0 in among-species mean growth rate
#sr: one row per time point. Columns are coefficients of omega^(1/omega) and omega^0 in among-species st dev of growth rates (there is also an infinitesimal part that isn't calculated)
#rnorm: one row per time point. Columns are coefficients of omega^(1/omega) and omega^0 in norm of growth rates (ignores infinitesimal part)
#cossqtheta: vector of real numbers: proportion of change that is size change (ignores infinitesimal part if there are colonizations and extinctions)
#cossqthetarank: ranks of cossqtheta (including real and infinitesimal parts
#rnormrank: rank of rnorm (including all parts)
#numcossqtheta, denomcossqtheta: numerator and denominator of cos^2 theta: each three columns, coefficients of omega^{2/omega}, omega^{1/omega}, omega^0.
getsurr <- function(data,times){
  deltat <- diff(times)
  s <- dim(data)[2]
  ntimes <- dim(data)[1]
  sr <- array(dim=c(ntimes,2))
  rbar <- array(dim=c(ntimes,2))
  rnorm <- array(dim=c(ntimes,2))
  cossqtheta <- array(dim=c(ntimes,1))
  numcossqtheta <- array(dim=c(ntimes,3))
  denomcossqtheta <- array(dim=c(ntimes,3)) 
  lnd <- log(data)

  lnb <- lnd[1:ntimes-1,,drop=FALSE]
  lna <- lnd[2:ntimes,,drop=FALSE]

  for(tt in 1:(length(times)-1)){
    a <- lna[tt,]
    b <- lnb[tt,]
    S1 <- a!= -Inf & b!= -Inf
    S2 <- a== -Inf & b!= -Inf
    S3 <- a!= -Inf & b== -Inf
    S4 <- a== -Inf & b== -Inf
    k4 <- sum(S4)
    k3 <- sum(S3)
    k2 <- sum(S2)
    as <- sum(a[S1 | S3])
    bs <- sum(b[S1 | S2])

    rbar[tt+1,1] <- (k3-k2)/(s*deltat[tt]) #coefficient of omega^(1/omega)
    rbar[tt+1,2] <- (as-bs)/(s*deltat[tt]) #coefficient of omega^0

    if((k3==0) & (k2==0)){#no colonizations or extinctions
      sr[tt+1,1] <- 0
      meanr <- c((a[S1]-b[S1])/deltat[tt],rep(0,k4)) #append zeros for species absent everywhere
      sr[tt+1,2] <- sd(meanr)

      #if zeros appended to meanr as above, the following is identical:
#      gamma <- sum((a[S1]-b[S1])^2)-(as-bs)^2/s
#      sr[tt+1,2] <- sqrt(gamma)/(sqrt(s-1)*deltat[tt])
      gammap <- sum((a[S1]-b[S1])^2)
      rnorm[tt+1,1] <- 0
      rnorm[tt+1,2] <- sqrt(gammap)/deltat[tt]
      cossqtheta[tt+1] <- s*mean(meanr)^2/gammap
      numcossqtheta[tt+1,] <-c(0,0,s*mean(meanr)^2)
      denomcossqtheta[tt+1,] <- c(0,0,gammap)

    } else {
      #standard deviation of growth rates
        alpha <- k2+k3-1/s*(k3-k2)^2
        beta <- 2*(sum(b[S2])+sum(a[S3])-1/s*(as-bs)*(k3-k2))
        sr[tt+1,1] <- sqrt(alpha)/(sqrt(s-1)*deltat[tt])
        sr[tt+1,2] <- beta/(2*sqrt(alpha)*sqrt(s-1)*deltat[tt])

        #norm of growth rates divided by number of species
        alphap <- k2+k3
        betap <- 2*(sum(b[S2])+sum(a[S3]))
        gammap <- sum((a[S1]-b[S1])^2)+sum(b[S2]^2)+sum(a[S3]^2
                                                        )
        rnorm[tt+1,1] <- sqrt(alphap)/(deltat[tt])
        rnorm[tt+1,2] <- betap/(2*sqrt(alphap)*deltat[tt])

        #what proportion of change is size change? ignores infinitesimal parts
        cossqtheta[tt+1] <- (k3-k2)^2/(s*(k2+k3))
        numcossqtheta[tt+1,] <-c((k3-k2)^2,2*(as-bs)*(k3-k2),as-bs)
        denomcossqtheta[tt+1,] <- s*c(alphap,betap,gammap)
      }
  }
  csrat <- list(num=numcossqtheta[-1,],denom=denomcossqtheta[-1,])
  cossqthetarank <- c(NA,rankfromsortsurrr(bubblesortsurrrat(csrat)))
  rnormrank <- c(NA,ranksurr(denomcossqtheta[-1,]))
  return(list(rbar=rbar,sr=sr,rnorm=rnorm,cossqtheta=cossqtheta,cossqthetarank=cossqthetarank,rnormrank=rnormrank,numcossqtheta=numcossqtheta,denomcossqtheta=denomcossqtheta))
}

#plot a time series of infinite and real parts of a size or shape measure of rate of change
#surr is an array: rows are times, first column is infinite part (coefficient of omega^(1/omega)), second column is real part.
#time is a vector of time points
#x.lab and y.lab are axis labels
#panel.label is a label for top left corner
plotsurr <- function(surr,time,x.lab,y.lab,panel.label){
  



    par(mgp=c(2.5,0.5,0)) #axis title, labels, line
  ascale <- 0.1
  ry <- c(range(surr[,1],na.rm=TRUE),range(surr[,1]+ascale*surr[,2],na.rm=TRUE))
  rx <- c(range(time),range(time))
  plot(rx,ry,type="n",xlab="",ylab="",cex.axis=1,yaxt="n")
  axl <- par("usr")


    mtext(x.lab,side=1,line=2,cex=0.75)
    mtext(y.lab,side=2,line=3.5,cex=0.75)


    mytick <- pretty(range(surr[,1],na.rm=TRUE),2) #not too many ticks  
  mytick[mytick<axl[3]] <- NA
  mytick[mytick>axl[4]] <- NA
  ytp <- axis(2,col.axis="transparent",at=mytick) 
  labs <- lapply(ytp,function(x) bquote(.(x)~psi))

#  labs <- lapply(ytp,function(x) bquote(.(x)~"\U0001D713"))  #unicode italic psi: displays fine but doesn't export to eps
  
  axis(2,at=ytp,labels=do.call(expression,labs))
  lines(time,surr[,1],col=colors()[320])
  points(time,surr[,1],pch=16,cex=0.5)
  arrows(time,surr[,1],time,surr[,1]+ascale*surr[,2],length=0.05)

  par(xpd=TRUE)

  arx <- axl[1]-0.04*(axl[2]-axl[1])
  artick <- max(abs(surr[,2]),na.rm=TRUE)
  ary <- ascale*artick
  yt <- pretty(c(-ary,ary),n=3) #tick marks for scale of real part
  arrows(arx,0,arx,max(yt),length=0.05)
  arrows(arx,0,arx,min(yt),length=0.05)  
  ytrim <- yt[c(-1,-length(yt))]
  points(rep(arx,length(ytrim)),ytrim,pch=3)
  text(rep(arx,length(ytrim)),ytrim,ytrim/ascale,cex=1,pos=2,offset=0.25)

#  xtx <- min(time)+0.5
#  ryy <- range(ry)
#  text(xtx,ryy[2]-0.025*diff(ryy),panel.label,cex=1.5)
  axl <- par("usr")
  text(axl[1]+0.01*(axl[2]-axl[1]),axl[4]-0.07*(axl[4]-axl[3]),panel.label,cex=1.5,pos=4)  
  
}

#scatter plot shape against size, either real or infinite parts
#x and y are coordinates
#x.lab and y.lab are axis labels
#panel.label is label for panel
#if isinf is TRUE, label tick marks in terms of psi
#if jit is true, jitter points (useful for infinite parts)
plotpair <- function(x,y,x.lab,y.lab,panel.label,isinf,jit){

  if(jit){
    x <- jitter(x)
    y <- jitter(y)
  }
  if(isinf){#plotting infinite part
    plot(x,y,type="n",xlab="",ylab="",cex.axis=1,xaxt="n",yaxt="n")
    mtext(x.lab,side=1,line=3,cex=1)
    mtext(y.lab,side=2,line=4,cex=1)
    ytp <- axis(2,col.axis="transparent") #get current tick mark locations
    ylabs <- lapply(ytp,function(x) bquote(.(x)~psi))
    axis(2,at=ytp,labels=do.call(expression,ylabs))
    xtp <- axis(1,col.axis="transparent") #get current tick mark locations
    xlabs <- lapply(xtp,function(x) bquote(.(x)~psi))
    axis(1,at=xtp,labels=do.call(expression,xlabs))
    points(x,y)
  }
  else {#plotting finite part
    plot(x,y,xlab=x.lab,ylab=y.lab,cex.axis=1)
  }

#  xx <- range(x,na.rm=TRUE)
#  yy <- range(y,na.rm=TRUE)
#  text(xx[1]+0.1*diff(xx),yy[2]-0.1*diff(yy),panel.label,cex=1.5)
  axl <- par("usr")
  text(axl[1]+0.01*(axl[2]-axl[1]),axl[4]-0.07*(axl[4]-axl[3]),panel.label,cex=1.5,pos=4)  
  
}

plotsurrray <- function(rbar,sr,nspp,panel.label){
#ray diagram of changes in rates of change of size and shape
#rbar is array of rates of change in size (first column infinite part, second column real part)
#sr is array of rates of change in shape (same arrangement)
#nspp is number of species included
#panel.label is label for top LH corner

  par(mgp=c(3.2,0.5,0)) #axis title, labels, line
  
  lr <- dim(rbar)[1]
  x <- diff(rbar[,1]*sqrt(nspp)) #infinite part
  y <- diff(sr[,1]*sqrt(nspp-1))
  xr <- diff(rbar[,2]*sqrt(nspp)) #real part
  yr <- diff(sr[,2]*sqrt(nspp-1))

  ascale <- 0.1

  rx <- c(range(x,na.rm=TRUE),range(x+ascale*xr,na.rm=TRUE)) #good axis limits
  ry <- c(range(y,na.rm=TRUE),range(y+ascale*yr,na.rm=TRUE))
  
  xl <- c(-max(abs(rx)),max(abs(rx)))
  yl <- c(-max(abs(ry)),max(abs(ry)))
  
  plot(rx,ry,type="n",asp=1,xlim=xl,ylim=yl,cex.lab=1,cex.axis=1,xlab=bquote(Delta(italic(n)^{1/2}~bar(italic(r)))),ylab=bquote(Delta(italic(n)-1)^{1/2}~italic(s[r])),xaxt="n",yaxt="n")

  axl <- par("usr")
  mytick <- pretty(range(y,na.rm=TRUE),2) #not too many ticks  
  mytick[mytick<axl[3]] <- NA
  mytick[mytick>axl[4]] <- NA
  ytp <- axis(2,col.axis="transparent",at=mytick) 
  ylabs <- lapply(ytp,function(x) bquote(.(x)~psi))
  axis(2,at=ytp,labels=do.call(expression,ylabs))
  
  mytick <- pretty(range(x,na.rm=TRUE),2) #not too many ticks  
  mytick[mytick<axl[1]] <- NA
  mytick[mytick>axl[2]] <- NA
  xtp <- axis(1,col.axis="transparent",at=mytick)
  xlabs <- lapply(xtp,function(x) bquote(.(x)~psi))
  axis(1,at=xtp,labels=do.call(expression,xlabs))
  for(i in 1:(lr-1)){
    lines(x=c(0,x[i]),y=c(0,y[i]),col=colors()[320])
    arrows(x[i],y[i],x[i]+ascale*xr[i],y[i]+ascale*yr[i],length=0.05)    
  }
  points(x,y,pch=16,cex=0.5)


  
  par(xpd=TRUE) #scale for real part
  arx <- axl[1]-0.12*(axl[2]-axl[1])
  artick <- 4*max(abs(ry),na.rm=TRUE)
  ary <- ascale*artick
  yt <- pretty(c(-ary,ary),n=3) #tick marks for scale of real part
  arrows(arx,0,arx,max(yt),length=0.05)
  arrows(arx,0,arx,min(yt),length=0.05)  
  ytrim <- yt[c(-1,-length(yt))]
  points(rep(arx,length(ytrim)),ytrim,pch=3)
  text(rep(arx,length(ytrim)),ytrim,ytrim/ascale,cex=1,pos=2,offset=0.3)


  ary <- axl[3]-0.13*(axl[4]-axl[3])
  artick <- 4*max(abs(rx),na.rm=TRUE)
  arx <- ascale*artick
  xt <- pretty(c(-arx,arx),n=3) #tick marks for scale of real part
  arrows(0,ary,max(xt),ary,length=0.05)
  arrows(0,ary,min(xt),ary,length=0.05)  
  xtrim <- xt[c(-1,-length(xt))]
  points(xtrim,rep(ary,length(xtrim)),pch=3)
  text(xtrim,rep(ary,length(xtrim)),xtrim/ascale,cex=1,pos=1,offset=0.5)

#  text(axl[1]+0.05*(axl[2]-axl[1]),axl[4]-0.06*(axl[4]-axl[3]),panel.label,cex=1.5)
  axl <- par("usr")
  text(axl[1]+0.01*(axl[2]-axl[1]),axl[4]-0.07*(axl[4]-axl[3]),panel.label,cex=1.5,pos=4)  
  
}

plotsurrcircle <- function(rbar,sr,nspp,panel.label){
#circle diagram of rates of change of size and shape
#rbar is array of rates of change in size (first column infinite part, second column real part)
#sr is array of rates of change in shape (same arrangement)
#nspp is number of species included
#panel.label is label for top LH corner

  par(mgp=c(3.2,0.5,0)) #axis title, labels, line

  colors <- brewer.pal(5,"Dark2")[2:3]
  xs=rbar[,1]*sqrt(nspp)
  xr=rbar[,2]*sqrt(nspp)
  ys=sr[,1]*sqrt(nspp-1)
  yr=sr[,2]*sqrt(nspp-1)

  ascale <- 0.1
  rx <- c(range(xs,na.rm=TRUE),range(xs+ascale*xr,na.rm=TRUE)) #good axis limits
  ry <- c(range(ys,na.rm=TRUE),range(ys+ascale*yr,na.rm=TRUE))
  
  xl <- c(-max(abs(rx)),max(abs(rx)))
  yl <- c(0,max(abs(ry)))

  plot(rx,ry,type="n",asp=1,xlim=xl,ylim=yl,cex.lab=1,cex.axis=1,xlab=bquote(italic(n)^{1/2}~bar(italic(r))),ylab=bquote((italic(n)-1)^{1/2}~italic(s[r])),xaxt="n",yaxt="n")

  par(xpd=FALSE)
  theta <- c(seq(from=0,to=pi,length.out=1000))
  for(r in 1:3){
    x <- r*cos(theta)
    y <- r*sin(theta)
    lines(t(x),t(y),lty="dashed",col=colors()[320])
  }

#blue when infinite part of rbar decreasing, orange when increasing
  lr <- dim(rbar)[1]
  for(i in 1:(lr-1)){
    if(xs[i+1] < xs[i]) col=colors[2]
    else if (xs[i+1]>xs[i]) col=colors[1]
    else col=colors()[320]
    lines(x=xs[i:(i+1)],y=ys[i:(i+1)],col=col,lwd=1)
    arrows(xs[i],ys[i],xs[i]+ascale*xr[i],ys[i]+ascale*yr[i],length=0.05)    
    if(i==1) points(xs[1],ys[1],pch=16,col=col)
  }
  points(xs[lr],ys[lr],pch=1,col=col)
  arrows(xs[lr],ys[lr],xs[lr]+ascale*xr[lr],ys[lr]+ascale*yr[lr],length=0.05)    
  
  axl <- par("usr")
  mytick <- pretty(range(ys,na.rm=TRUE),2) #not too many ticks  
  mytick[mytick<axl[3]] <- NA
  mytick[mytick>axl[4]] <- NA
  ytp <- axis(2,col.axis="transparent",at=mytick) 
  ylabs <- lapply(ytp,function(x) bquote(.(x)~psi))
  axis(2,at=ytp,labels=do.call(expression,ylabs))
  
  mytick <- pretty(range(xs,na.rm=TRUE),2) #not too many ticks  
  mytick[mytick<axl[1]] <- NA
  mytick[mytick>axl[2]] <- NA
  xtp <- axis(1,col.axis="transparent",at=mytick)
  xlabs <- lapply(xtp,function(x) bquote(.(x)~psi))
  axis(1,at=xtp,labels=do.call(expression,xlabs))

  par(xpd=TRUE) #scale for real part
  arx <- axl[1]-0.1*(axl[2]-axl[1])
  artick <- 2*max(abs(ry),na.rm=TRUE)
  ary <- ascale*artick
  yt <- pretty(c(-ary,ary),n=3) #tick marks for scale of real part
  arrows(arx,0,arx,max(yt),length=0.05)
  arrows(arx,0,arx,min(yt),length=0.05)  
  ytrim <- yt[c(-1,-length(yt))]
  points(rep(arx,length(ytrim)),ytrim,pch=3)
  text(rep(arx,length(ytrim)),ytrim,ytrim/ascale,cex=1,pos=2,offset=0.3)

  ary <- axl[3]-0.13*(axl[4]-axl[3])
  artick <- 5*max(abs(rx),na.rm=TRUE)
  arx <- ascale*artick
  xt <- pretty(c(-arx,arx),n=3) #tick marks for scale of real part
  arrows(0,ary,max(xt),ary,length=0.05)
  arrows(0,ary,min(xt),ary,length=0.05)  
  xtrim <- xt[c(-1,-length(xt))]
  points(xtrim,rep(ary,length(xtrim)),pch=3)
  text(xtrim,rep(ary,length(xtrim)),xtrim/ascale,cex=1,pos=1,offset=0.5)

#  text(axl[1]+0.05*(axl[2]-axl[1]),axl[4]-0.06*(axl[4]-axl[3]),panel.label,cex=1.5)
  axl <- par("usr")
  text(axl[1]+0.01*(axl[2]-axl[1]),axl[4]-0.07*(axl[4]-axl[3]),panel.label,cex=1.5,pos=4)  
 
}

plotlogcount <- function(year,data,n,panel.label,blackline,xoff=0.01,cexax=1,cexlab=0.75,col=colors()[320]){
#plot log counts+1 against time for Surtsey data
#year is survey year
#data is a data frame, in which first n columns are counts for each species (rows are years)
#panel.label is label for top left corner
#blackline is index of species to highlight with a black line
#xoff is offset of label from left edge
#cexax is size for axis tick labels
#cexlab is size for axis labels
#col is colour for lines  
  yl <- range(log(data[,1:n]+1))
  plot(year,log(data[,1]),type="n",ylim=yl,xlab="",ylab="",cex.axis=cexax)
  mtext("Year",side=1,line=3,cex=cexlab) #use mtext so size matches surreal number plots
  mtext(bquote(log(italic(c[i])+1)),side=2,line=3,cex=cexlab)
  for(i in 1:n){
    lines(year,log(data[,i]+1),col=col)
  }
  if(!is.na(blackline)) lines(year,log(data[,blackline]+1)) #discussed by name in text
  axl <- par("usr")
  text(axl[1]+xoff*(axl[2]-axl[1]),axl[4]-0.07*(axl[4]-axl[3]),panel.label,cex=1.5,pos=4)
}

#activity vs proportion of change that is "size" change
#ssur is output from getsurr
#nspp is number of species
#panel.label is text for top left-hand corner
plottempvpropsurr <- function(ssur,nspp,panel.label){
  par(mgp=c(3.2,0.5,0)) #axis title, labels, line
  ascale <- 0.1

  xs <- 1/sqrt(nspp)*ssur$rnorm[,1]
  xr <- 1/sqrt(nspp)*ssur$rnorm[,2]
  rx <- c(range(xs,na.rm=TRUE),range(xs+ascale*xr,na.rm=TRUE)) #good axis limits
  ry <- c(0,1,0,1)
  xl <- c(0,max(abs(rx)))

  plot(rx,ry,type="n",xlim=xl,ylim=c(0,1),xaxt="n",yaxp=c(0,1,2),cex.lab=0.83,cex.axis=1,ylab=bquote(Proportion~'"size"'~change),xlab=bquote(Scaled~activity~(years^-1))) #0.83 because par(mfrow) reduces fonts by this amount
  axl <- par("usr")
  mytick <- pretty(range(rx,na.rm=TRUE),4) #not too many ticks  
  mytick[mytick<axl[1]] <- NA
  mytick[mytick>axl[2]] <- NA
  xtp <- axis(1,col.axis="transparent",at=mytick)
  xlabs <- lapply(xtp,function(x) bquote(.(x)~psi))
  axis(1,at=xtp,labels=do.call(expression,xlabs))

  ary <- axl[3]-0.13*(axl[4]-axl[3])
  artick <- max(abs(rx),na.rm=TRUE)
  arx <- ascale*artick
  xt <- pretty(c(-arx,arx),n=3) #tick marks for scale of real part
  arrows(0,ary,max(xt),ary,length=0.05)
  arrows(0,ary,min(xt),ary,length=0.05)  
  xtrim <- xt[c(-1,-length(xt))]
  points(xtrim,rep(ary,length(xtrim)),pch=3)
  text(xtrim,rep(ary,length(xtrim)),xtrim/ascale,cex=1,pos=1,offset=0.5)
  points(xs,ssur$cossqtheta)
  arrows(xs,ssur$cossqtheta,xs+ascale*xr,ssur$cossqtheta,length=0.05)    
  text(axl[1]+0.01*(axl[2]-axl[1]),axl[4]-0.07*(axl[4]-axl[3]),panel.label,cex=1.5,pos=4)
  
}

#rank a set of surreal numbers in the form \sum_i a_i psi^i
#surr is an array: columns are coefficients for powers of psi in descending order from left to right
ranksurr <- function(surr){
    rx <- rank(do.call(interaction,args=list(data.frame(surr),lex.order=TRUE)))
    return(rx)
}

#compare two ratios of surreal numbers
largerr <- function(pairnum,pairdenom) {
    a1 <- pairnum[1,1]
    b1 <- pairnum[1,2]
    c1 <- pairnum[1,3]
    d1 <- pairdenom[1,1]
    e1 <- pairdenom[1,2]
    f1 <- pairdenom[1,3]
    a2 <- pairnum[2,1]
    b2 <- pairnum[2,2]
    c2 <- pairnum[2,3]
    d2 <- pairdenom[2,1]
    e2 <- pairdenom[2,2]
    f2 <- pairdenom[2,3]
    tC <- cbind(a2*d1,a2*e1+b2*d1,a2*f1+b2*e1+c2*d1,b2*f1+c2*e1,c2*f1)
    tD <- cbind(a1*d2,a1*e2+b1*d2,a1*f2+b1*e2+c1*d2,b1*f2+c1*e2,c1*f2)
    cdr <- ranksurr(rbind(tC,tD))

    #sign of product of denominators
    sgs <- c(d1*d2,d1*e2+e1*d2,d1*f2+e1*e2+f1*d2,e1*f2+f1*e2,f1*f2)
    ts <- which(!sgs==0)
    if(length(ts)>0) xms <- sign(sgs[ts[1]]) else xms <- 0 #at least one of the ratios has a zero denominator: need to catch this at entry, not here
    
    if(xms*cdr[1] < xms*cdr[2]) return(TRUE) else return(FALSE)
}

#used in bubble sort
swapiflargerr <- function(pairnum,pairdenom,pairindex) {
    if(largerr(pairnum,pairdenom)) {
        return(list(num=pairnum[c(2,1),],denom=pairdenom[c(2,1),],index=pairindex[c(2,1)]))
    } else {
        return(list(num=pairnum,denom=pairdenom,index=pairindex))
    }
}

#used in bubble sort
swappassr <- function(surrr) {
    for(i in seq(1, dim(surrr$num)[1]-1)) {
        sswap <- swapiflargerr(surrr$num[i:(i+1),],surrr$denom[i:(i+1),],surrr$index[i:(i+1)])
        surrr$num[i:(i+1),] <- sswap$num
        surrr$denom[i:(i+1),] <- sswap$denom
        surrr$index[i:(i+1)] <- sswap$index
    }
    return(surrr)
}

#call this via bubblesortsurrrat()
##based on http://www.numbertheory.nl/2013/05/10/bubble-sort-implemented-in-pure-r/
bubblesortsurrr <- function(surrr){
    newsurrr <- swappassr(surrr)
    if(isTRUE(all.equal(surrr, newsurrr))) {
        return(newsurrr)
    } else {
        return(bubblesortsurrr(newsurrr))
    }
}


#surrr is a list of two surreal number arrays num, denom
#each one consists of 3 columns representing descending powers of omega
#each row is the coefficients for a single surreal number
#want to sort the ratios num/denom
#return a list containing num and denom sorted, and the sort index
#unless some ratios have zero denominator, in which case prints warning and returns NA
bubblesortsurrrat <- function(surrr) {
    if(any(rowSums(abs(surrr$denom))==0)){
        print("some ratios have zero denominator")
        return(NA)
    }
    surrr$index <- seq(from=1,to=dim(surrr$num)[1])
    return(bubblesortsurrr(surrr))
}

#get ranks from sorted set of surreal ratios 
#input a list (output from bubblesortsurrrat) containing
#num and denom (three columns each, descending powers of omega)
# andindex: index for sort order
#return ranks (using mean rank for ties)
rankfromsortsurrr <- function(surrr){
  library(Matrix)
  require(Matrix)
    n <- dim(surrr$num)[1]
    c <- 1
    ranks <- array(dim=c(n,1))
    ranks[1] <- 1
    for(i in 2:n){
        if(largerr(surrr$num[c(i,i-1),],surrr$denom[c(i,i-1),])) c <- c+1
        ranks[i] <- c
    }
    ranks <- rank(ranks) #fractional ranks to ties
    ranks <- ranks[invPerm(surrr$index)] #order ranks to match original data
}
