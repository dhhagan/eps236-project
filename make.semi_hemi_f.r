semi.hemi.forward=function(species="SF6",tglobal=tglobal.all[species] , lplot=F,delt=.025,t.range=c(1995,2009-delt), S.scale=1.1,
#taustrat=4,tau.intrahemis=1/2,tau.interhemis=2 ,strat.nh.fraction=.55,tglobal=tglobal.all[species] ,
params=c(.5,2.3,.55,4),emissions=NULL,c.init=NULL,saveK=F,l.orig=F)
{
tau.intrahemis=params[1];tau.interhemis=params[2] ;strat.nh.fraction=params[3];taustrat=params[4]
# # require to run: 5_box_model_initialize.r to make non-parameter objects; comments are in this file and in "make_halocarbon_boxdata.r"
if(!exists("Mtrop"))stop("You need to run first; 5_box_model_initialize.r")
#      
c.obs=ghg.gmd[,c("Year", paste(species,"box",1:4,sep="."))] 
tglobal=tglobal.all[species]
Mol.wt = Mol.wt.all[species]
L.freq= loss.monthly[,paste(species,"box",1:5,sep=".")];rownames(L.freq)=1:12  ##monthly loss freq
L.freq.ann=apply(L.freq,2,mean) ## annual loss frequencies
#print(L.freq)
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
    kadv.ntrop2strop=FT2S*(2*strat.nh.fraction - 1)/2  #=kadv.strat2nh - kadv.strat2sh
#
#
  alpha=(1/tglobal + 1/taustrat)/(1/taustrat -1/tglobal*.842/.158)
   tstrat= tglobal*.158/(.842*alpha+.158) #  chemical-loss-lifetime-in-strat
   header=(paste(species,"Mol.wt",Mol.wt,"; tglobal=",tglobal,"taustrat=",taustrat,"strat.nh.fraction=",
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
K1j=c((-1)*(kadv.strat2nh + kih),kih,0, 0, kadv.strat2nh )
K2j=c(kadv.strat2nh + kih,(-1)*(kih + kadv.tr2strat + kadv.ntrop2strop + kh),kh,0,0)
K3j=c(0,(kadv.ntrop2strop + kh),(-1)*(kih + kadv.tr2strat + kh),(kadv.strat2sh+kih),0)
K4j=c(0,0,kih,(-1)*(kadv.strat2sh + kih), kadv.strat2sh)
#K5j=c(0, kadv.tr2strat, kadv.tr2strat,0,-(kadv.strat2nh + kadv.strat2sh+Mstrat/tstrat) )
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
Gg2ppt= 1e9/Mol.wt.all[species]/(Mglobal/.029)*1e12
if(is.null(emissions)){
## if species!="SF6" then specify emissions via emissions=<obj> in arg list; sf6.sources give the sources for SF6 from EDGAR
# Format is annual in each tropospheric box : > sf6.sources
##     yr     box.1     box.2     box.3     box.4
##1  1970 0.7212501 0.0019824 0.0003115 0.0012474
##2  1971 0.8182210 0.0340240 0.0017166 0.0020934
##3  1972 0.8543246 0.0601936 0.0017468 0.0018853
##4  ......
S.SPECIES=t(cbind(sf6.sources[,2:5],rep(0,nrow(sf6.sources)) ))*Gg2ppt*S.scale ; dimnames(S.SPECIES)=list(1:5,sf6.sources[,1])
#note EDGAR adds 2.9 ppt bet 1995 and end 2008, but atm adds 3.23, so we should scale up by 10% +-
} else {
S.SPECIES=emissions
}
#
Yrs=seq(t.range[1],t.range[2],delt)  ## time period to be considered for this run
l.start=c.obs[,"Year"]>=Yrs[1]
l.end=c.obs[,"Year"]<=t.range[2]
#Initial conditions
if(is.null(c.init)){
## estimate init strat as 3 years later .... (not optimum)
i.start=min(which(l.start))
c0=as.numeric(c.obs[i.start,paste(species,"box",1:4,sep=".")])
c4=as.numeric(c.obs[i.start+38,paste(species,"box",1:4,sep=".")]) ##based on eigenvalue, choose 3.2 yr delay
C0=c(c0,2*mean(c0)-mean(c4))
}else{
C0=c.init }
#

# Integrate the equations using the simple time stepping

CC=matrix(C0,ncol=1,nrow=5)## ; AA=matrix(rep(0,5),nrow=5,ncol=1) ; CCi=CC
for(yr in Yrs){
if(l.orig){K1 = KK.orig} else {
##  The next 3 lines add the chemical loss term to the diagonal elements, can vary monthly
MM=trunc( (yr - trunc(yr) )*12) +1
L.freq.i=L.freq[MM,]*3.15e7  ## monthly loss frequencies converted into yr-1
K1=KK-diag(L.freq.i)}
#
## The forward model equation:  Delta.C = ( K %*% C  + P  )*Delta.t
del.CC=(K1%*%CC[,ncol(CC)]+S.SPECIES[,as.character(trunc(yr))]/c(Mtrop,Mstrat)*Mglobal)  *delt
CC.new=CC[,ncol(CC)]+del.CC
CC=cbind(CC,CC.new)
}
dimnames(CC)=list(paste("Box",1:5,sep="."),c(Yrs[1]-delt,Yrs) )


if(lplot){
 matplot(c(Yrs[1]-delt,Yrs) ,t(CC),type="l",lwd=4,xlab="Year",ylab=(paste(species," (ppt)")),lty=1,
                 main=paste("Model and Observations",species))
yy=ghg.gmd[,"Year"]
matpoints(yy,ghg.gmd[,paste(species,"box",1:4,sep=".")],pch=3,cex=.4)
 legend("topleft",legend=c(paste("Box",1:5,sep="."),paste("Obs",1:4,sep="-")),lty=c(1,1,1,1,1,-1,-1,-1,-1),
pch=c(rep(-1,5),rep(3,4)),col=1:5)

}
if(saveK)write.table(KK,file=paste("KK",params[1],params[2],params[3],params[4],".tbl",sep="_"),row.names=F,col.names=F,quote=F)

return(CC)

}


