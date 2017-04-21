## seasonal variations:  where we make serial correlation the focus of our science
sf6.global=read.table("HATS_global_SF6.txt",header=T,na.strings=c("NA","nan"))
#what is in the data
head(sf6.global)
##what is in a time series
help(ts)
sf6.ts=ts(sf6.global,start=1995-1/24,frequency=12)
plot(sf6.ts[,"HATS_NH_SF6"])

library(feather)
library(data.table)

sf6.df <- data.frame(sf6.ts)

feather::write_feather(sf6.df, "HATS_global_SF6.feather")

sf6.df.feather <- feather::read_feather("HATS_global_SF6.feather")

##are there NAs ?
sum(is.na(sf6.ts[,"HATS_NH_SF6"]))
help(stl)
library(stlplus)
## is there a seasonal cycle, how does it change with time ??
#
zum=stl(sf6.ts[,"HATS_NH_SF6"],s.window=15)
plot(zum)
tail(sf6.ts)
lok=sf6.ts[,1]%in%1995:2015
zam=tapply(sf6.ts[lok,"HATS_NH_SF6"],sf6.ts[lok,2],mean)
dev.new();plot(zam,type="o")
t7=lowess(sf6.ts[lok,"HATS_NH_SF6"],f=.1)
dev.new();plot(t7)
zam=tapply(sf6.ts[lok,"HATS_NH_SF6"]-t7$y[lok],sf6.ts[lok,2],mean)
plot(zam,type="o")
## NOW a short workshop using the famous BRW and MLO data -- how does the amplitude of the seasonal cycle change with time?
## Thoning method .....
## stl or stlplus (with na's).   note:  na.action=na.exclude may give an answer (e.g. in ar()) where na.omit gives an error.
dev.new()
par(mfrow=c(2,2))
plot(sf6.ts[,"HATS_NH_SF6"]-sf6.ts[,"HATS_SH_SF6"],main="Northern - Southern Hemisphere",ylab=expression(paste(Delta,SF[6])))## !!!!!
plot(sf6.ts[,"HATS_mlo_SF6"]-sf6.ts[,"HATS_spo_SF6"],main="Mauna Loa - South Pole",ylab=expression(paste(Delta,SF[6])))## !!!!!
plot(sf6.ts[,"HATS_mlo_SF6"]-sf6.ts[,"HATS_cgo_SF6"],main="Mauna Loa - Cape Grim",ylab=expression(paste(Delta,SF[6])))## !!!!!
plot(sf6.ts[,"HATS_alt_SF6"]-sf6.ts[,"HATS_spo_SF6"],main="Alert - South Pole",ylab=expression(paste(Delta,SF[6])))## !!!!!

#South Pole (SPO, 90°S, 2837 m asl)
#Palmer Station, Antarctica (PSA, 64.6°S, 64.0°W, 10 m asl)
#Cape Grim, Australia (CGO, 40.682°S, 144.688°E, 164 m asl; inlet is 70 m agl)
#American Samoa (SMO, 14.247°S, 170.564°W, 77 m asl)
#Mauna Loa, USA (MLO, 19.5362°N, 155.5763°W, 3397 m asl)
#Cape Kumukahi, USA (KUM, 19.516°N, 154.811°W, 3 m asl)
#Niwot Ridge, USA (NWR, 40.1°N, 105.5°W, 3475 m asl)
#Trinidad Head, USA (THD, 41.0°N, 124.1°W, 120 m asl)
#Wisconsin, USA (LEF, 45.6°N, 90.27°W, 868 m asl; inlet is 396 m above ground)
#Harvard Forest, USA (HFM, 42.5°N, 72.2°W, 340 m asl; inlet is 29 m above ground)
#Mace Head, Ireland (MHD, 53.3°N, 9.9°W, 42 m asl)
#Barrow, USA (BRW, 71.3°N, 156.6°W, 8 m asl)
#Alert, Canada (ALT, 82.5°N, 62.3°W, 210 m asl)
#Summit, Greenland (SUM, 72.6°N, 38.4°W, 3200 m asl)
##   alt  sum  brw  mhd  thd  nwr  kum  mlo  smo  cgo  ush  psa spo
## alert summit barrow mace_head trinidad_head niwot_ridge kumukahi mauna_loa AmSamoa cape_grim Ushuaia PalmerSta SPole

##latKey = "brw:71.3;sum:72.5;nwr:40.04;mlo:19.5;smo:-14.3;spo:-67;alt:82.45;cgo:-40.68;kum:19.52;mhd:53.33;psa:-64.92;thd:41.05;lef:45.95;ush:-54.87;hfm:42.54;gmi:13.386;mid:28.210;asc:-7.967;eic:-27.160"
