semi.hemi.inv=function(species="N2O", lplot=F,delt=.025,t.range=c(1996+.042,2009-delt), S.scale=1.1,
l.forward=F, params=c(0.15,0.8,0.55,4),c.init=NULL,saveK=F,l.orig=F ) {  
#                      ^    ^    ^  ^    optimized values
# require to run: 5_box_model_initialize.r to make non-parameter objects; comments are in this file and in "make_halocarbon_boxdata.r"
if(!exists("Mtrop"))stop("Must first run 5_box_model_initialize.r")
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## tau values, exchange times for mass between the boxes
tau.intrahemis=params[1];tau.interhemis=params[2] ;strat.nh.fraction=params[3];taustrat=params[4]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c.obs=ghg.gmd[,c("Year", paste(species,"box",1:4,sep="."))] 
tglobal=tglobal.all[species]
Mol.wt = Mol.wt.all[species]
L.freq= loss.monthly[,paste(species,"box",1:5,sep=".")];rownames(L.freq)=1:12  ##monthly loss freq in s^-1
L.freq.ann=apply(L.freq,2,mean) ## annual mean loss frequencies  in s^-1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#      
######################## definition of exch rates kg/yr ################################
    kih=Mtrop[1]/tau.intrahemis               #Midlat-Tropics exchange flux kg/yr
    kh=(Mtrop[1] + Mtrop[2])/tau.interhemis               #N/S exchange flux kg/yr

#The advection flux is defined by taustrat and strat.nh.fraction; FT2S=total mass flux Mstrat/taustrat
    FT2S=Mstrat/taustrat
    kadv.tr2strat=FT2S/2  # 2: equal each hemisph tropics
    kadv.strat2nh=FT2S*strat.nh.fraction
    kadv.strat2sh=FT2S*(1 - strat.nh.fraction )
    kadv.ntrop2strop=FT2S*(2*strat.nh.fraction - 1)/2  #= ( kadv.strat2nh - kadv.strat2sh )/2 to balance the advection
#
#
  alpha=(1/tglobal + 1/taustrat)/(1/taustrat -1/tglobal*.842/.158)
   tstrat= tglobal*.158/(.842*alpha+.158) #  chemical-loss-lifetime-in-strat
   header=(paste(species,"Mol.wt",Mol.wt,", tglobal=",tglobal,"taustrat=",taustrat,"strat.nh.fraction=",
   strat.nh.fraction,"tau.intrahemis=",tau.intrahemis,"tau.interhemis=",tau.interhemis,"alpha=",alpha,"tstrat=",tstrat))
#
##################################Start definition of K matrix #########################################
# equations: (Mass of air / .029 )= moles of air; so e.g. box 1  d(M1/.029*c1)/dt = -kadv/.029*c1 No. moles/time: 0.029 in comm.
#1 NH.hi d(M1*c1) = (-1)*(kadv.strat2nh + kih)*c1 +                                          kih*c2                                                             + kadv.strat2nh*c5
#2 NH.tr d(M2*c2) =      (kadv.strat2nh + kih)*c1 +(-1)*(kih + kadv.tr2strat + kadv.ntrop2strop + kh)*c2 +                             kh *c3
#3 SH.tr d(M2*c3) =                                                          +(kadv.ntrop2strop + kh)*c2 + (-1)*(kih + kadv.tr2strat + kh)*c3 + (kadv.strat2sh+kih)*c4
#4 SH.hi d(M4*c4) =                                        + kih*c3                                                            + (-1)*(kadv.strat2sh + kih)*c4  + kadv.strat2sh*c5
#5 Stratosphere   =                                                                     kadv.tr2strat*c2 +                   kadv.tr2strat*c3     -(kadv.strat2nh + kadv.strat2sh)*c5
##              c1            c2    c3 c4    c5
## add diagonal, losst Mtrop[j]/tautrop <<=
K1j=c((-1)*(kadv.strat2nh + kih),kih,0, 0, kadv.strat2nh )
K2j=c(kadv.strat2nh + kih,(-1)*(kih + kadv.tr2strat + kadv.ntrop2strop + kh),kh,0,0)
K3j=c(0,(kadv.ntrop2strop + kh),(-1)*(kih + kadv.tr2strat + kh),(kadv.strat2sh+kih),0)
K4j=c(0,0,kih,(-1)*(kadv.strat2sh + kih), kadv.strat2sh)
K5j=c(0, kadv.tr2strat, kadv.tr2strat,0,-(kadv.strat2nh + kadv.strat2sh) ) ##chemical loss to be added below
K0=rbind(K1j,K2j,K3j,K4j,K5j)
KK=rbind(K1j/Mtrop[1],K2j/Mtrop[2],K3j/Mtrop[3],K4j/Mtrop[4],K5j/Mstrat)
## holdover fo checking
K5j.orig=c(0, kadv.tr2strat, kadv.tr2strat,0,-(kadv.strat2nh + kadv.strat2sh+Mstrat/tstrat) )
KK.orig=rbind(K1j/Mtrop[1],K2j/Mtrop[2],K3j/Mtrop[3],K4j/Mtrop[4],K5j.orig/Mstrat)
# KK=K0/Mglobal ##units change from kg/yr to yr-1; we prefer not to solve for masses although maybe we should.....
## static: done in test.sf6.r : e.EE=eigen(KK); lambda=e.EE$values; Xe=e.EE$vectors; Xe1=solve(Xe)
#
##simple integration, year by year #sources, annual, in ppt/yr
#Gg2ppt= 1e9/1000/Mglobal*29/Mol.wt.all[species]*1e12
Gg2ppt= 1e9/Mol.wt.all[species]/(Mglobal/.029)*1e12

Yrs=seq(t.range[1],t.range[2],delt)  ## time period to be inverted
l.start=c.obs[,"Year"]>=Yrs[1]-mean(diff(c.obs[,"Year"]))
l.end=c.obs[,"Year"]<=t.range[2]+mean(diff(c.obs[,"Year"])) ## ensure span in approx interpolation
#Initial conditions
if(is.null(c.init)){
## estimate init strat as 3 years later .... (not optimum)
i.start=min(which(l.start))
c0=as.numeric(c.obs[i.start,paste(species,"box",1:4,sep=".")])
c4=as.numeric(c.obs[i.start+38,paste(species,"box",1:4,sep=".")]) ##based on eigenvalue, choose 3.2 yr delay
C0=c(c0,2*mean(c0)-mean(c4))
}else{
C0=c.init }
if(C0[5]<0)C0[5]=0
#
CC=matrix(C0,ncol=1,nrow=5)## ; AA=matrix(rep(0,5),nrow=5,ncol=1) ; CCi=CC
#
# ------------------------------------------------------------------------------------------------------
# interpolate concentrations to integration grid
# we have no data for stratosphere.....
CCC=cbind(C0,matrix(NA,nrow=5,ncol=length(Yrs)))
for(i in 1:4){
CCC[i,1+1:length(Yrs)]=approx(x=c.obs[l.start&l.end,"Year"],y=c.obs[l.start&l.end,paste(species,"box",i,sep=".")],xout=Yrs)$y
}
colnames(CCC)=c(Yrs[1]-delt,Yrs)

PP=NULL
for(i in 1+1:length(Yrs)){
#solve for stratosphere with requirement P5 = 0
#    use annual loss rate, form global lifetime, or input monthly losses derived from Prabir
if(l.orig){K1 = KK.orig} else {
# the next three lines add the monthly loss terms to the diagonal of the KK matrix
MM=trunc( (Yrs[i] - trunc(Yrs[i]) )*12) +1
L.freq.i=L.freq[MM,]*3.15e7  ## monthly loss frequencies converted into yr-1
K1=KK-diag(L.freq.i)}
#
CCC[5,i]=CCC[5,i-1]+ c(K1%*%CCC[,i-1])[5]*delt
## The inverse model equation: P(t.i) = dC.i/dt - K %*% C
ppi=(CCC[,i]-CCC[,i-1])/delt - K1%*%CCC[,i]
PP=cbind(PP,ppi)
}
PP.Tg.yr=PP*c(Mtrop,Mstrat)/Mglobal/Gg2ppt  ## units were ppt/yr in each box; switch to ppb *1000  then to Tg /1000
#
#
if(lplot){ 
dev.new()
UU="Tg/yr";if(species=="N2O")UU="TgN/yr";if(species%in%c("SF6","CFC11","CFC12","HCFC22"))UU="Mg/yr"
matplot(Yrs,t(PP.Tg.yr),type="o",lty=1,lwd=3,pch=16,cex=.4,xlab="Year",ylab=paste(species,UU) )
legend("topleft",legend=paste("Box",1:5,sep="."),col=1:5,lty=1,lwd=3,pch=16) 
abline(h=0,lty=2)
dev.copy(png,paste("Fig.",species,"_sources.png",sep=""));dev.off()
dev.new()
PP.yr=aggregate(t(PP.Tg.yr),list(trunc(Yrs)),mean,na.rm=T)
matplot(PP.yr[,1],(PP.yr[,2:5]),type="o",lty=1,lwd=3,pch=16,cex=.4,xlab="Year",ylab=paste(species,UU) )
legend("topright",legend=paste("Box",1:5,sep="."),col=1:5,lty=1,lwd=3,pch=16) 
abline(h=0,lty=2)
if(species%in%c("CFC12","CFC11") ){ axis(side=2,at=seq(-625,625,25),tck=.015,labels=F)
 axis(side=4,at=seq(-625,625,25),tck=.015,labels=F)
 axis(side=4,at=seq(-600,600,100),tck=.025,labels=F)}
dev.copy(png,paste("Fig.",species,"_sources_annual.png",sep=""));dev.off()
PP.seas=aggregate(t(PP.Tg.yr),list(Yrs-trunc(Yrs)),mean,na.rm=T)
matplot(PP.seas[,1],(PP.seas[,2:5]),type="o",lty=1,lwd=3,pch=16,cex=.4,xlab="Season of Year",ylab=paste(species,UU) )
legend("topright",legend=paste("Box",1:5,sep="."),col=1:5,lty=1,lwd=3,pch=16) 
abline(h=0,lty=2)
dev.copy(png,paste("Fig.",species,"_sources_seasonal.png",sep=""));dev.off()


#
}
colnames(PP)=Yrs
colnames(PP.Tg.yr)=Yrs;rownames(PP.Tg.yr)=1:5
return(PP.Tg.yr)
}


