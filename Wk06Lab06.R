options(repos ="http://cran.us.r-project.org")  # required to get latest libs
# Installing the packages we will play with today
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster)
# From previous weeks
pacman::p_load(EcoHydRology,rnoaa,curl,httr)
rm(list=objects())
setwd("~/Wk06Lab06/")
###Source TMWB model
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")
##source CNmodel function
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/CNmodel")

# Download a soils dataset for your basin based on the WebSoilSurvey method 
# and replace this url with your own
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/eei1fn32t45stkpd5nuhbnnh/wss_aoi_2022-03-03_14-26-07.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")
# use the function to get data from USGS 0205551460 
#LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
myflowgage_id="0205551460"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2022-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
mysoil=readOGR("wss_aoi_2022-03-03_14-26-07/spatial/soilmu_a_aoi.shp")    
# Explore the mysoil dataset which is returned
mybbox=c(mysoil@bbox)
# First associate mukey with cokey from component
mysoil$mukey=mysoil$MUKEY  # or rename the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
co2ch = SDA_query(q_co2ch)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r
mu2ch=merge(mu2co,co2ch)
mu2chmax=aggregate(mu2ch,list(mu2ch$mukey),max)
proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem=get_elev_raster(locations=mysoil, 
                      z = 11, prj =proj4string(mysoil) ,
                      src ="aws",clip="bbox",expand = 0.001)

stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
# We are looking for stations with elements that have PRCP, TMAX and TMIN 
# and current data (i.e. Year 2021). 
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
#summary(WXData)  #

# Create an aligned modeldata data frame to build our model in
modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm
# View(modeldata)  
# Compare your precipitation to the flow out of your basin
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=
  modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=
  modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MinTemp)/2.0
modeldata[is.na(modeldata)]=0 # A Quick BUT sloppy removal of NAs
TMWB=modeldata

#####################################
# Homework 1
#####################################
################################
# Calibrate CN model parameters
###############################
CNmodeldf = modeldata

f <- function (x) {
  CNopt=x[1]
  IaOpt=x[2]
  DAWCopt = x[3]
  zOpt = x[4]
  fcresOpt = x[5]
  CNmodelnew = CNmodel(CNmodeldf = CNmodeldf, CNavg = CNopt, IaFrac = IaOpt, 
                     func_DAWC = DAWCopt, func_z = zOpt, fnc_fcres = fcresOpt)
  return(1-NSE(CNmodelnew$Qmm,CNmodelnew$Qpred))
}

awcmin=min(mu2chmax$awc_r, na.rm  = T); awcmax=max(mu2chmax$awc_r, na.rm  = T); 
zmin = min(mu2chmax$hzdepb_r, na.rm  = T)*10; 
zmax = max(mu2chmax$hzdepb_r, na.rm  = T)*10
lower <- c(35, .01, awcmin, zmin, .1)
upper <- c(99, .25, awcmax, zmax, .5)
## run DEoptim and set a seed first for replicability
set.seed(1234)
DEoptim(f, lower, upper, control = DEoptim.control(itermax = 50))

detach(CNmodeldf)
#$optim$bestmem
#CNopt         IaOpt         DAWCopt         zOpt         fcresOpt 
#9.604356e+01 2.262338e-02 2.157368e-01 1.875752e+03 3.747216e-01 


CNmodelnew <- CNmodel(CNmodeldf, CNavg = 96.04, IaFrac = 0.022, 
                     func_DAWC = 0.215, func_z = 1876, fnc_fcres = 0.375 )
detach(CNmodeldf)
################################
# Calibrate TMWBmodel parameters
###############################
TMWB = modeldata

f2 <- function (x) {
  fcresOpt=x[1]
  SFTmpOpt=x[2]
  TlagOPt = x[3]
  AWCopt = x[4]
  
  TMWBnew = TMWBmodel(TMWB=TMWB, fcres=fcresOpt, SFTmp=SFTmpOpt, bmlt6=2.5,
                      bmlt12=1, Tlag=TlagOPt, AWCval=AWCopt, Slope=0)
  return(1-NSE(TMWBnew$Qmm, TMWBnew$Qpred))
}

awcmin=min(mu2chmax$awc_r, na.rm  = T); awcmax=max(mu2chmax$awc_r, na.rm  = T); 
zmin = min(mu2chmax$hzdepb_r, na.rm  = T)*10;
awcmin = awcmin*zmin
awcmax = awcmax*zmax
lower <- c(0.1, -5, 0.1, awcmin)
upper <- c(0.5, 5, 1, awcmax)


## run DEoptim and set a seed first for replicability
set.seed(1234)
DEoptim(f2, lower, upper, control = DEoptim.control(itermax = 50))

detach(TMWB)
#$optim$bestmem
#fcres        SFTmp         Tlag         AWC
#0.2856447  4.4838149  0.8145922 80.7296812 

TMWBnew = TMWBmodel(TMWB=TMWB, fcres=0.285, SFTmp=4.48, bmlt6=2.5,
                    bmlt12=1, Tlag=0.814, AWCval=80.73, Slope=0)
detach(TMWB)

##########Homework 1 Plotting ################

NSETMWB <- NSE(TMWBnew$Qmm,TMWBnew$Qpred)
NSECN <- NSE(CNmodelnew$Qmm,CNmodelnew$Qpred)
pacman::p_load(ggplot2)
Phw1 <- ggplot() +
  geom_line(data=TMWBnew,aes(x=date, y = Qmm,colour="Qmm")) +
  geom_line(data=TMWBnew,aes(x=date, y = Qpred,colour="Qpred_TMWB,NSE=0.346")) +
  geom_line(data=CNmodelnew,aes(x=date, y = Qpred,colour="Qpred_CN,NSE=0.567")) +
  labs(x = 'Date', y = 'Flow (mm)') +
  scale_colour_manual("", 
                      breaks = c("Qmm", "Qpred_TMWB,NSE=0.346", "Qpred_CN,NSE=0.567"),
                      values = c("black", "blue","red")) +
  theme(text = element_text(size = 10)) +
  ggtitle("Discharge Comparison between CN model 
          and TMWB against Observed data")
pdf("Phw1.pdf")
Phw1
dev.off()
NSECN
NSETMWB
#####################################
# Homework 2
#####################################
pacman::p_load(lubridate, data.table)
BasinTMWB_JO=TMWBnew[(month(TMWBnew$date) >= 5    #June to October
                      & month(TMWBnew$date) < 11),]
S_low <- (1000/95-10)*25.4
S_high <- (1000/37-10)*25.4
attach(BasinTMWB_JO)
f <- function (x) {
  Sest=x
  NSE(Qmm,dP^2/(dP+Sest))
}

Sest=optimize(f, c(S_low,S_high), tol = 0.0001,maximum = TRUE)$maximum
plot(dP,Qmm)
points(dP,dP^2/(dP+Sest),col="red") 
NSE(Qmm,dP^2/(dP+Sest))
detach(BasinTMWB_JO)

############### Topographic Index Classes #################
nTIclass=5
VSAsol=data.table(WetClass=seq(from=nTIclass,to=1),
                  As=seq(1:nTIclass)*(1/nTIclass),Wetfrac=(1/nTIclass))
VSAsol[,sSratio:=2*(sqrt(1-shift(As))-sqrt(1-As))/Wetfrac-1]
#
# Inspect what the previous command gives us, note it is just a fancy way of 
# shifting the index of a vector in the VSAsol data frame 
# using the data.table::shift() function.
#
VSAsol 
#
# Now fill in the missing value
#
VSAsol$sSratio[1]=2*(sqrt(1-0)-sqrt(1-VSAsol$As[1]))/VSAsol$Wetfrac[1]-1
#
VSAsol[,sigma:=Sest*sSratio]
VSAsol[,CN:=25400/(sigma+254)]
VSAsol
plot(VSAsol$As,VSAsol$sigma)
lines(VSAsol$As,VSAsol$sigma)
plot(VSAsol$As,VSAsol$CN)
lines(VSAsol$As,VSAsol$CN)

############## TIC Model #####################################
# Initialize the TI Class objects from top to bottom of slope
TIC05=modeldata
TIC04=modeldata
TIC03=modeldata
TIC02=modeldata
TIC01=modeldata
# For TIC05 CNavg=VSAsol$CN[1]
TIC01 = CNmodel(CNmodeldf = TIC01, CNavg=VSAsol$CN[5],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)

TIC02$P=TIC01$Excess+TIC02$P
TIC02 = CNmodel(CNmodeldf = TIC02, CNavg=VSAsol$CN[4],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)
TIC03$P=TIC02$Excess+TIC03$P
TIC03 = CNmodel(CNmodeldf = TIC03, CNavg=VSAsol$CN[3],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)
TIC04$P=TIC03$Excess+TIC04$P
TIC04 = CNmodel(CNmodeldf = TIC04, CNavg=VSAsol$CN[2],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)
TIC05$P=TIC04$Excess+TIC05$P
TIC05 = CNmodel(CNmodeldf = TIC05, CNavg=VSAsol$CN[1],
                func_DAWC=.3,IaFrac=0.05,
                func_z=1000,fnc_fcres=.3)

# Scale reservoir coefficient between the .2-.5 given in class
# Same as lab 5, just call the CNmodel function with each of the
# TIClass objects and route Qpred to ExcessIn below
#
# First, let's look at the differences in flow for different CNs
Phw2 <- ggplot() +
  geom_line(data=TIC05,aes(x=date, y = Qmm, linetype = "dashed")) +
  geom_line(data=TIC05,aes(x=date, y = Qpred*0.2,colour="Qpred_TIC05")) +
  geom_line(data=TIC04,aes(x=date, y = Qpred*0.2,colour="Qpred_TIC04")) +
  geom_line(data=TIC03,aes(x=date, y = Qpred*0.2,colour="Qpred_TIC03")) +
  geom_line(data=TIC02,aes(x=date, y = Qpred*0.2,colour="Qpred_TIC02")) +
  geom_line(data=TIC01,aes(x=date, y = Qpred*0.2,colour="Qpred_TIC01")) +
  labs(x = 'Date', y = 'Flow (mm)', linetype = "Observed flow", colour = "Estimated")+
  ggtitle("simple CN-VSA model")

pdf("Phw2.pdf")
Phw2
dev.off()
Average_runoff <- data.frame(TIC5 = mean(TIC05$Qpred*0.2), TIC4 = mean(TIC04$Qpred*0.2), 
                             TIC3 = mean(TIC03$Qpred*0.2), TIC2 = mean(TIC02$Qpred*0.2),
                             TIC1 = mean(TIC01$Qpred*0.2))
################################
### Homework 03 ##############
Phw3 <- ggplot() +
  
  geom_line(data=TIC05,aes(x=date, y = AW,colour="AW_TIC05")) +
  geom_line(data=TIC04,aes(x=date, y = AW,colour="AW_TIC04")) +
  geom_line(data=TIC03,aes(x=date, y = AW,colour="AW_TIC03")) +
  geom_line(data=TIC02,aes(x=date, y = AW,colour="AW_TIC02")) +
  geom_line(data=TIC01,aes(x=date, y = AW,colour="AW_TIC01")) +
  labs(x = 'Date', y = 'AW (mm)', colour = "Available Water")+
  ggtitle("simple CN-VSA model") +
  theme_classic()

pdf("Phw3.pdf")
Phw3
dev.off()
Average_AW <- data.frame(TIC5meanAW = mean(TIC05$AW), TIC4meanAW = mean(TIC04$AW), 
                         TIC3meanAW = mean(TIC03$AW), TIC2meanAW = mean(TIC02$AW),
                         TIC1meanAW = mean(TIC01$AW), 
                         TIC5minAW = min(TIC05$AW), TIC4minAW = min(TIC04$AW), 
                         TIC3minAW = min(TIC03$AW), TIC2minAW = min(TIC02$AW),
                         TIC1minAW = min(TIC01$AW),
                         TIC5maxAW = max(TIC05$AW), TIC4maxAW = max(TIC04$AW), 
                         TIC3maxAW = max(TIC03$AW), TIC2maxAW = max(TIC02$AW),
                         TIC1maxAW = max(TIC01$AW))

#############################
#####Graduate HW #######
Phw4 <- ggplot() +
    geom_line(data=TIC01,aes(x=date, y = ET,colour="ET_TIC01")) +
  geom_line(data=TIC02,aes(x=date, y = ET, colour="ET_TIC02")) +
  geom_line(data=TIC03,aes(x=date, y = ET, colour="ET_TIC03")) +
  geom_line(data=TIC04,aes(x=date, y = ET, colour="ET_TIC04")) +
  geom_line(data=TIC05,aes(x=date, y = ET, colour="ET_TIC05")) +
  labs(x = 'Date', y = 'ET (mm)', colour = "Estimated")+
  ggtitle("simple CN-VSA model") +
  theme_classic()

pdf("Phw4.pdf")
Phw4
dev.off()
Average_ET <- data.frame(TIC5 = mean(TIC05$ET), TIC4 = mean(TIC04$ET), 
                             TIC3 = mean(TIC03$ET), TIC2 = mean(TIC02$ET),
                             TIC1 = mean(TIC01$ET))
