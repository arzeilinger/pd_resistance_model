if (!"deSolve" %in% installed.packages()) install.packages("deSolve")
library(deSolve)

one.run<-function(propn,ff=2,harvest=50,c1=1,c2=1,plot.it=T){
  
  params2<-c(
    T=1,
    muC=0.001,
    muI=0.01,
    aC=1,
    aI=1,
    pC=0.001,
    pI=0.01,
    bC=1,
    bI=1,
    eta=0.1,
    delta=0.1,
    gamma=0.1,
    lambda=0.083,
    fudge=ff,
    tau=harvest
  )
  
  times<-seq(0,harvest,length=101)
  times.end<-length(times)
  xstart2<-c(S1=100*(1-propn),E1=0,Hc1=0,Hi1=0,S2=100*(propn),E2=0,Hc2=0,Hi2=0,U=99,Vc=1,Vi=0)
  
  xyl.2<-ode(xstart2,times,xylella.ode2,params2)
  
  if (plot.it){
    par(mfrow=c(3,2),mar=c(2,3,1,1),
        oma=c(4,0.1,2,1),las=1,font=3,
        tcl=-0.2,pty="m",cex=1,las=1,cex.axis=0.6,cex.lab=0.8,mgp=c(2,0.5,0))
    matplot(times,cbind(xyl.2[,2]),type="l",ylim=c(0,100),xlab="",ylab="S (S)")
    matplot(times,cbind(xyl.2[,3:5]),type="l",ylim=c(0,100),xlab="",ylab="E, Hc, Hi, E+Hi (S)")
    lines(times,xyl.2[,3]+xyl.2[,5],type="l",lwd=2,col=1)
    matplot(times,cbind(xyl.2[,6]),type="l",ylim=c(0,100),xlab="",ylab="S (R)")
    matplot(times,cbind(xyl.2[,7:9]),type="l",ylim=c(0,100),xlab="",ylab="E, Hc, Hi, E+Hi (R)")
    lines(times,xyl.2[,7]+xyl.2[,9],type="l",lwd=2,col=1)
    matplot(times,cbind(xyl.2[,10]),type="l",ylim=c(0,100),xlab="Time",ylab="U")
    matplot(times,cbind(xyl.2[,11:12]),type="l",ylim=c(0,100),xlab="Time",ylab="Vc, Vi")
    fig.end()
    profit.x<-(c1*(xyl.2[,2])+c2*(xyl.2[,6]))[times.end]
    mtext(paste("Total profit:",round(profit.x,2)),side=3,cex=1,font=1,line=2)
    mtext("Time",side=1,cex=1,font=1,line=1)
  }
  return(profit.x)
}

# harvest time
harvest<-200
# latent period scaling
ff<-10
# price for each variety; healthy and latent, not infected
c1<-1
c2<-0.1

one.run(0.5,ff,harvest,c1,c2,plot.it=T)

seed.propn<-seq(0,10,by=1)/10
seed.profit<-rep(NA,length=length(seed.propn))
for (i in 1:length(seed.propn)){
  seed.profit[i]<-one.run(seed.propn[i],ff,harvest,c1,c2,plot.it=T)
  cat(seed.propn[i],round(seed.profit[i],1),"\n")
}

par(mfrow=c(1,1),mar=c(3,3,1,1),
    mgp=c(2.7,0.5,0),
    oma=c(5,0,2,1),las=1,font=3,tcl=-0.2,pty="m",
    cex=1,las=1,cex.axis=0.6,cex.lab=0.8,mgp=c(2,0.5,0))
plot(seed.propn,seed.profit,type="b",pch=16,cex=1,cex.axis=1,
     xlab="Proportion of (R)esilient crop",
     ylab="Profit")
