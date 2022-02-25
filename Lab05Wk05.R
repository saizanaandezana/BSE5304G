if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology, rnoaa, elevatr,soilDB,rgdal,raster, rgeos, tidyverse, patchwork, aqp)
setwd(dir = "C:/Users/binyam/OneDrive - Virginia Tech/Advanced Watershed Modelling/Week 5/Wk05Lab05")
myflowgage_id = "0205551460"
myflowgage = get_usgs_gage(myflowgage_id, begin_date = "2015-01-01", end_date = "2022-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

stns = meteo_distance(
  station_data = ghcnd_stations(),
  lat = myflowgage$declat,
  long = myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
# We are looking for stations with elements that have PRCP, TMAX and TMIN 
# and current data (i.e. Year 2021). 
WXStn = stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData = meteo_pull_monitors(
  monitors = WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)

# Create an aligned modeldata data frame to build our model in
modeldata = merge(WXData, myflowgage$flowdata, by.x = "date", by.y = "mdate")
summary(modeldata)  #
modeldata$MaxTemp = modeldata$tmax/10 # Converting to C
modeldata$MinTemp = modeldata$tmin/10 # Converting to C
modeldata$P = modeldata$prcp/10 # Converting to mm
View(modeldata)  
# Compare your precipitation to the flow out of your basin
modeldata$P[is.na(modeldata$P)] = 0
modeldata$MinTemp[is.na(modeldata$MinTemp)] = 0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)] =
  modeldata$MinTemp[is.na(modeldata$MaxTemp)] + 1
modeldata$MaxTemp[modeldata$MaxTemp <= modeldata$MinTemp] =
  modeldata$MinTemp[modeldata$MaxTemp <= modeldata$MinTemp] + 1
modeldata$AvgTemp = (modeldata$MaxTemp + modeldata$MinTemp)/2.0

summary(modeldata)
modeldata[is.na(modeldata)] = 0 # A Quick BUT sloppy removal of NAs

########### CN amd TMWB functions##########################################
CNmodel <- function(CNmodeldf, CNavg = 75, IaFrac = 0.05, fnc_slope = 0, 
                    fnc_aspect=0, func_DAWC = .3, func_z = 1000, fnc_fcres = .3) {
    # Energy Balance based Snow Accumulation 
  # and Melt model from the EcoHydRology package.
  attach(CNmodeldf)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1, 
                      windHt = 2, groundAlbedo = 0.25,SurfEmissiv = 0.95, windSp = 2, 
                      forest = 0, startingSnowDepth_m = 0,startingSnowDensity_kg_m3=450)
  # We will update the -3 in the above to be a lapse rate adjustment
  detach(CNmodeldf)
  CNmodeldf$SNO=SNO_Energy$SnowWaterEq_mm
  CNmodeldf$SNOmlt=SNO_Energy$SnowMelt_mm
  CNmodeldf$SnowfallWatEq_mm=SNO_Energy$SnowfallWatEq_mm
  CNmodeldf$SnowMelt_mm=SNO_Energy$SnowMelt_mm
  attach(CNmodeldf)
  CNmodeldf$Albedo=.23
  CNmodeldf$Albedo[CNmodeldf$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),
                   Tmax_C = MaxTemp,Tmin_C = MinTemp,
                   lat_radians = myflowgage$declat*pi/180) * 1000
  CNmodeldf$PET=PET
  detach(CNmodeldf)
  rm(list="PET")
  
  CNmodeldf$AWC=func_DAWC*func_z
  # Oh, this we want to vary some of these around our watershed!
  CNmodeldf$dP = 0 # Initializing Net Precipitation
  CNmodeldf$ET = 0 # Initializing ET
  CNmodeldf$AW = 0 # Initializing AW
  CNmodeldf$Excess = 0 # Initializing Excess
  CNmodeldf$S =0 # Initializing S
  CNmodeldf$Qpred=0 # Initializing Qpred
  attach(CNmodeldf)
  SSCNavg=(1000/CNavg-10)*25.4
  SSCN=SoilStorage(S_avg=SSCNavg, field_capacity=func_DAWC*.9,
                   soil_water_content=0.1*func_DAWC, porosity=func_DAWC)
  Ia_init=IaFrac*SSCN   
  CNmodeldf$CNavg = CNavg
  CNmodeldf$SSCNavg = SSCNavg
  CNmodeldf$SSCN = SSCN
  detach(CNmodeldf)
  rm(list=c("CNavg", "SSCN", "SSCNavg"))
  CNmodeldf$Ia = Ia_init
  attach(CNmodeldf)
  # Those processes that are dependant on prior days conditions, we run as a 
  # loop through each of the days.
  for (t in 2:length(AW)){
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    # Calculating Net Precipitation which adds in slope above's Excess
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] + 
      SNO_Energy$SnowMelt_mm[t]    # CN Solution
    # Is the soil saturated, and thus can't take more dP? 
    if (AW[t-1] + dP[t]>=AWC[t]){
      Excess[t]=AW[t-1] + dP[t] -AWC[t]
      AW[t]=AWC[t]
      # Otherwise, if dP is less than the initial abstraction? 
      # https://en.wikipedia.org/wiki/Runoff_curve_number#Definition
    } else if (dP[t]<=Ia[t]) {
      Excess[t]=0.0
      AW[t]=AW[t-1] + dP[t]
    } else {
      Excess[t]=(dP[t]-Ia[t])^2/(dP[t]-Ia[t]+SSCN[t])
      AW[t]=AW[t-1] + dP[t] -Excess[t]
    }
    S[t]=S[t-1]+Excess[t]
    Qpred[t]=fnc_fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  CNmodeldf$ET=ET
  CNmodeldf$dP=dP
  CNmodeldf$AW=AW
  CNmodeldf$Excess=Excess
  CNmodeldf$S=S
  CNmodeldf$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  rm(list=c("AW", "dP", "ET", "Excess", "Qpred", "S"))
  detach(CNmodeldf)
  return(CNmodeldf)
}

source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")

####################################Soil and Topo data#####################
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/bky5q4i2wkdjmvm3t0d1gniq/wss_aoi_2022-02-11_09-31-24.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")
mysoil = readOGR("wss_aoi_2022-02-11_09-31-24/spatial/soilmu_a_aoi.shp")
# Explore the mysoil dataset which is returned

mybbox = c(mysoil@bbox)
# First associate mukey with cokey from component
mysoil$mukey = mysoil$MUKEY  # or rename the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey, cokey FROM component WHERE mukey IN ", mukey_statement, sep = "")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r

mu2ch = merge(mu2co,co2ch)
summary(mu2ch)
mu2chmax = aggregate(mu2ch,list(mu2ch$mukey), max)

proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem = get_elev_raster(locations = mysoil, 
                        z = 11, prj = proj4string(mysoil) ,
                        src ="aws",clip="bbox",expand = 0.001)

summary(terrain(mydem, opt='TRI',unit = "radians"))
plot(terrain(mydem, opt='TRI',unit = "radians"))

####################################CNModel#################################

TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata
##################Gather soil and topo parameters###########

summary(mu2chmax)
##topslope has the lowest Z, lowest slope as botslope, midslope has the #highest slope, we assume that awc in % is 
#equal for all three models
awcpercent=mean(mu2chmax$awc_r,na.rm=TRUE)
### for TopSlope model lowest depth
TopslopeZ=min(mu2chmax$hzdepb_r)*10 #in mm
MidslopeZ=mean(mu2chmax$hzdepb_r)*10
BotslopeZ=max(mu2chmax$hzdepb_r)*10
#####################AWCVAL is mm of awc= depth*awc(%)
AWCvalTop=TopslopeZ*awcpercent
AWCvalMid=MidslopeZ*awcpercent
AWCvalBot=BotslopeZ*awcpercent
####################################calculation of slope from terrain #function,note that the unit is in degree
slop_stat = summary(terrain(mydem, opt='Slope', unit = "radians"))
SlopeTop = slop_stat[[2]] #radians
SlopeBot = slop_stat[[2]]
SlopeMid = slop_stat[[5]]
# Running the model with the three HRU dataframes
# Low slope but highest ksat
# These 3 function calls are what you will vary for the Lab 04 homework 

# Call the new CNmodel() function with Top,Mid,BotSlope HRU objects,
# passing the Qpred into the lower HRUs HillslopeAboveExcess (as area scaled flow)
#TopSlopeCN = CNmodel(TopSlopeCN, CNavg = 60)
TopSlopeCN = CNmodel(CNmodeldf = TopSlopeCN, IaFrac = 0.05, CNavg = 60, fnc_slope = SlopeTop,
                     fnc_aspect = 0, func_DAWC=awcpercent,
                     func_z = TopslopeZ, fnc_fcres = .3)
MidSlopeCN$P = TopSlopeCN$Excess + MidSlopeCN$P
# Higher slope, medium ksat, fcres=0.5 
MidSlopeCN = CNmodel(CNmodeldf = MidSlopeCN, IaFrac = 0.05, CNavg = 60,fnc_slope=SlopeMid, 
                     fnc_aspect=0,func_DAWC=awcpercent,
                     func_z=MidslopeZ,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlopeCN$P = MidSlopeCN$Excess + BotSlopeCN$P
BotSlopeCN = CNmodel(CNmodeldf = BotSlopeCN, IaFrac = 0.05, CNavg = 60,fnc_slope=SlopeBot, 
                     fnc_aspect=0,func_DAWC=awcpercent,
                     func_z=BotslopeZ,fnc_fcres=.2)

TopSlopeCN = mutate(TopSlopeCN, TopAW = TopSlopeCN$AW)
MidSlopeCN = mutate(MidSlopeCN, MidAW = MidSlopeCN$AW)
BotSlopeCN = mutate(BotSlopeCN, BotAW = BotSlopeCN$AW)
plotAWCN = left_join(TopSlopeCN[, c("date","TopAW")], 
                   MidSlopeCN[, c("date","MidAW")], by = "date") %>% 
  left_join(BotSlopeCN[, c("date","BotAW")]) %>% pivot_longer(-date)
p1 <- ggplot(plotAWCN, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Available Water") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

#FOR Excess
TopSlopeCN = mutate(TopSlopeCN, TopEx = TopSlopeCN$Excess)
MidSlopeCN = mutate(MidSlopeCN, MidEx = MidSlopeCN$Excess)
BotSlopeCN = mutate(BotSlopeCN, BotEx = BotSlopeCN$Excess)
plotExCN = left_join(TopSlopeCN[, c("date","TopEx")], 
                   MidSlopeCN[, c("date","MidEx")], by = "date") %>% 
  left_join(BotSlopeCN[, c("date","BotEx")]) %>% pivot_longer(-date)
p2 <- ggplot(plotExCN, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Excess") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

#pdf("awex2.pdf")
p1/p2
#dev.off()


# CN Model Performance 
plot(BotSlopeCN$date,BotSlopeCN$Qpred,type="l")
lines(BotSlopeCN$date,BotSlopeCN$Qmm,type="l", col = "red")
plot(BotSlopeCN$date,cumsum(BotSlopeCN$Excess),type="l", col = "red")
NSeff(BotSlopeCN$Qmm,BotSlopeCN$Qpred)

# finish building all the hillslope HRUs….

####################################TMWBM#################################
TopSlope = modeldata
MidSlope = modeldata
BotSlope = modeldata

TopSlope = TMWBmodel(TMWB = TopSlope,SFTmp = 1, 
                     AWCval = AWCvalTop,
                     Tlag = .5,fcres=.3,Slope = SlopeTop)
MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5 
MidSlope = TMWBmodel(TMWB = MidSlope,SFTmp = 1, 
                       AWCval = AWCvalMid,
                       Tlag = .5,fcres=0.5,Slope = SlopeMid)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlope = TMWBmodel(TMWB = BotSlope,SFTmp = 1, 
                       AWCval = AWCvalBot,
                       Tlag = .5,fcres=0.2,Slope = SlopeBot)
##############AW Plots
TopSlope = mutate(TopSlope, TopAW = TopSlope$AW)
MidSlope = mutate(MidSlope, MidAW = MidSlope$AW)
BotSlope = mutate(BotSlope, BotAW = BotSlope$AW)
plotAW = left_join(TopSlope[, c("date","TopAW")], 
                   MidSlope[, c("date","MidAW")], by = "date") %>% 
  left_join(BotSlope[, c("date","BotAW")]) %>% pivot_longer(-date)
p3 <- ggplot(plotAW, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Available Water") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

#FOR Excess
TopSlope = mutate(TopSlope, TopEx = TopSlope$Excess)
MidSlope = mutate(MidSlope, MidEx = MidSlope$Excess)
BotSlope = mutate(BotSlope, BotEx = BotSlope$Excess)
plotEx = left_join(TopSlope[, c("date","TopEx")], 
                   MidSlope[, c("date","MidEx")], by = "date") %>% 
  left_join(BotSlope[, c("date","BotEx")]) %>% pivot_longer(-date)
p4 <- ggplot(plotEx, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Excess") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

#pdf("awex2.pdf")
p3/p4
#dev.off()


# TMWB Model Performance 
plot(BotSlope$date,BotSlope$Qpred,type="l")
lines(BotSlopeCN$date, BotSlopeCN$Qpred, type="l", col = "red")
NSeff(BotSlope$Qmm,BotSlope$Qpred)

CompareQpred = data.frame(modeldata$date, modeldata$Qmm, BotSlope$Qpred, BotSlopeCN$Qpred)
colnames(CompareQpred) <- c("date", "Qmm", "TMWBQpred", "CNQpred")
CompareQpred = CompareQpred %>% pivot_longer(-date)
p5 <- ggplot(CompareQpred, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "CN & TMWB models") + theme_classic() + xlab(NULL) +
  scale_y_continuous(name = "Depth (mm)") + ggtitle("Lick Run")

#################################################################
#################### HW3 ######################################
####################################################################
myflowgage_id = "12115700"
myflowgage = get_usgs_gage(myflowgage_id, begin_date = "2015-01-01", end_date = "2022-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

stns = meteo_distance(
  station_data = ghcnd_stations(),
  lat = myflowgage$declat,
  long = myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
# We are looking for stations with elements that have PRCP, TMAX and TMIN 
# and current data (i.e. Year 2021). 
WXStn = stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData = meteo_pull_monitors(
  monitors = WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)

# Create an aligned modeldata data frame to build our model in
modeldata = merge(WXData, myflowgage$flowdata, by.x = "date", by.y = "mdate")
summary(modeldata)  #
modeldata$MaxTemp = modeldata$tmax/10 # Converting to C
modeldata$MinTemp = modeldata$tmin/10 # Converting to C
modeldata$P = modeldata$prcp/10 # Converting to mm
View(modeldata)  
# Compare your precipitation to the flow out of your basin
modeldata$P[is.na(modeldata$P)] = 0
modeldata$MinTemp[is.na(modeldata$MinTemp)] = 0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)] =
  modeldata$MinTemp[is.na(modeldata$MaxTemp)] + 1
modeldata$MaxTemp[modeldata$MaxTemp <= modeldata$MinTemp] =
  modeldata$MinTemp[modeldata$MaxTemp <= modeldata$MinTemp] + 1
modeldata$AvgTemp = (modeldata$MaxTemp + modeldata$MinTemp)/2.0

summary(modeldata)
modeldata[is.na(modeldata)] = 0 # A Quick BUT sloppy removal of NAs


####################################Soil and Topo data#####################
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/tosambsrx5wdif00sp153f22/wss_aoi_2022-02-17_13-53-03.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")
mysoil = readOGR("wss_aoi_2022-02-17_13-53-03/spatial/soilmu_a_aoi.shp")
# Explore the mysoil dataset which is returned

mybbox = c(mysoil@bbox)
# First associate mukey with cokey from component
mysoil$mukey = mysoil$MUKEY  # or rename the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey, cokey FROM component WHERE mukey IN ", mukey_statement, sep = "")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r

mu2ch = merge(mu2co,co2ch)
summary(mu2ch)
mu2chmax = aggregate(mu2ch,list(mu2ch$mukey), max)

proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem = get_elev_raster(locations = mysoil, 
                        z = 11, prj = proj4string(mysoil) ,
                        src ="aws",clip="bbox",expand = 0.001)

summary(terrain(mydem, opt='TRI',unit = "radians"))
plot(terrain(mydem, opt='TRI',unit = "radians"))

####################################CNModel#################################

TopSlopeCN2=modeldata
MidSlopeCN2=modeldata
BotSlopeCN2=modeldata
##################Gather soil and topo parameters###########

summary(mu2chmax)
##topslope has the lowest Z, lowest slope as botslope, midslope has the #highest slope, we assume that awc in % is 
#equal for all three models
awcpercent=mean(mu2chmax$awc_r,na.rm=TRUE)
### for TopSlope model lowest depth
TopslopeZ=min(mu2chmax$hzdepb_r)*10 #in mm
MidslopeZ=mean(mu2chmax$hzdepb_r)*10
BotslopeZ=max(mu2chmax$hzdepb_r)*10
#####################AWCVAL is mm of awc= depth*awc(%)
AWCvalTop=TopslopeZ*awcpercent
AWCvalMid=MidslopeZ*awcpercent
AWCvalBot=BotslopeZ*awcpercent
####################################calculation of slope from terrain #function,note that the unit is in degree
slop_stat = summary(terrain(mydem, opt='Slope', unit = "radians"))
SlopeTop = slop_stat[[2]] #radians
SlopeBot = slop_stat[[2]]
SlopeMid = slop_stat[[5]]
# Running the model with the three HRU dataframes
# Low slope but highest ksat
# These 3 function calls are what you will vary for the Lab 04 homework 

# Call the new CNmodel() function with Top,Mid,BotSlope HRU objects,
# passing the Qpred into the lower HRUs HillslopeAboveExcess (as area scaled flow)
#TopSlopeCN = CNmodel(TopSlopeCN, CNavg = 60)
TopSlopeCN2 = CNmodel(CNmodeldf = TopSlopeCN2, IaFrac = 0.05, CNavg = 60, fnc_slope = SlopeTop,
                      fnc_aspect = 0, func_DAWC=awcpercent,
                      func_z = TopslopeZ, fnc_fcres = .3)
MidSlopeCN2$P = TopSlopeCN2$Excess + MidSlopeCN2$P
# Higher slope, medium ksat, fcres=0.5 
MidSlopeCN2 = CNmodel(CNmodeldf = MidSlopeCN2, IaFrac = 0.05, CNavg = 60,fnc_slope=SlopeMid, 
                      fnc_aspect=0,func_DAWC=awcpercent,
                      func_z=MidslopeZ,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlopeCN2$P = MidSlopeCN2$Excess + BotSlopeCN2$P
BotSlopeCN2 = CNmodel(CNmodeldf = BotSlopeCN2, IaFrac = 0.05, CNavg = 60,fnc_slope=SlopeBot, 
                      fnc_aspect=0,func_DAWC=awcpercent,
                      func_z=BotslopeZ,fnc_fcres=.2)

TopSlopeCN2 = mutate(TopSlopeCN2, TopAW = TopSlopeCN2$AW)
MidSlopeCN2 = mutate(MidSlopeCN2, MidAW = MidSlopeCN2$AW)
BotSlopeCN2 = mutate(BotSlopeCN2, BotAW = BotSlopeCN2$AW)
plotAWCN = left_join(TopSlopeCN2[, c("date","TopAW")], 
                     MidSlopeCN2[, c("date","MidAW")], by = "date") %>% 
  left_join(BotSlopeCN2[, c("date","BotAW")]) %>% pivot_longer(-date)
p7 <- ggplot(plotAWCN, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Available Water") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

#FOR Excess
TopSlopeCN2 = mutate(TopSlopeCN2, TopEx = TopSlopeCN2$Excess)
MidSlopeCN2 = mutate(MidSlopeCN2, MidEx = MidSlopeCN2$Excess)
BotSlopeCN2 = mutate(BotSlopeCN2, BotEx = BotSlopeCN2$Excess)
plotExCN = left_join(TopSlopeCN2[, c("date","TopEx")], 
                     MidSlopeCN2[, c("date","MidEx")], by = "date") %>% 
  left_join(BotSlopeCN2[, c("date","BotEx")]) %>% pivot_longer(-date)
p8 <- ggplot(plotExCN, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Excess") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

#pdf("awex2.pdf")
p7/p8
#dev.off()


# CN Model Performance 
plot(BotSlopeCN2$date,BotSlopeCN2$Qpred,type="l")
lines(BotSlopeCN2$date,BotSlopeCN2$Qmm,type="l", col = "red")
plot(BotSlopeCN2$date,cumsum(BotSlopeCN2$Excess),type="l", col = "red")
NSeff(BotSlopeCN2$Qmm,BotSlopeCN2$Qpred)

# finish building all the hillslope HRUs….

####################################TMWBM#################################
TopSlope2 = modeldata
MidSlope2 = modeldata
BotSlope2 = modeldata

TopSlope2 = TMWBmodel(TMWB = TopSlope2,SFTmp = 1, 
                      AWCval = AWCvalTop,
                      Tlag = .5,fcres=.3,Slope = SlopeTop)
MidSlope2$P=TopSlope2$Excess+MidSlope2$P
# Higher slope, medium ksat, fcres=0.5 
MidSlope2 = TMWBmodel(TMWB = MidSlope2,SFTmp = 1, 
                      AWCval = AWCvalMid,
                      Tlag = .5,fcres=0.5,Slope = SlopeMid)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope2$P=MidSlope2$Excess+BotSlope2$P
BotSlope2 = TMWBmodel(TMWB = BotSlope2,SFTmp = 1, 
                      AWCval = AWCvalBot,
                      Tlag = .5,fcres=0.2,Slope = SlopeBot)
##############AW Plots
TopSlope2 = mutate(TopSlope2, TopAW = TopSlope2$AW)
MidSlope2 = mutate(MidSlope2, MidAW = MidSlope2$AW)
BotSlope2 = mutate(BotSlope2, BotAW = BotSlope2$AW)
plotAW = left_join(TopSlope2[, c("date","TopAW")], 
                   MidSlope2[, c("date","MidAW")], by = "date") %>% 
  left_join(BotSlope2[, c("date","BotAW")]) %>% pivot_longer(-date)
p9 <- ggplot(plotAW, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Available Water") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

#FOR Excess
TopSlope2 = mutate(TopSlope2, TopEx = TopSlope2$Excess)
MidSlope2 = mutate(MidSlope2, MidEx = MidSlope2$Excess)
BotSlope2 = mutate(BotSlope2, BotEx = BotSlope2$Excess)
plotEx = left_join(TopSlope2[, c("date","TopEx")], 
                   MidSlope2[, c("date","MidEx")], by = "date") %>% 
  left_join(BotSlope2[, c("date","BotEx")]) %>% pivot_longer(-date)
p10 <- ggplot(plotEx, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Excess") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

p9/p10

# TMWB Model Performance 
plot(BotSlope2$date,BotSlope2$Qpred,type="l")
lines(BotSlope2CN2$date, BotSlope2CN2$Qpred, type="l", col = "red")
NSeff(BotSlope2$Qmm,BotSlope2$Qpred)

# Plot Qpred, Qmm for CN and TMWB
CompareQpred = data.frame(modeldata$date, modeldata$Qmm, BotSlope2$Qpred, BotSlopeCN2$Qpred)
colnames(CompareQpred) <- c("date", "Qmm", "TMWBQpred", "CNQpred")
CompareQpred = CompareQpred %>% pivot_longer(-date)
p11 <- ggplot(CompareQpred, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "CN & TMWB models") + theme_classic() + xlab(NULL) +
  scale_y_continuous(name = "Depth (mm)") 

pdf("HW3TMWB.pdf")
p11
dev.off()

CompareQpred3 = data.frame(BotSlopeCN2$date, BotSlopeCN2$Qmm, BotSlope2$Qpred, BotSlopeCN2$Qpred)
colnames(CompareQpred3) <- c("date", "Qmm", "TMWBQpred","CNQpred")
CompareQpred3 = CompareQpred3 %>% pivot_longer(-date)
CompareQpred4 = data.frame(BotSlopeCN2$date, BotSlopeCN2$SnowfallWatEq_mm, BotSlopeCN2$SnowMelt_mm)
colnames(CompareQpred4) <- c("date", "CNSnow Water Equ.", "CNSnow Melt")
CompareQpred4 = CompareQpred4 %>% pivot_longer(-date)
p12 <- ggplot(CompareQpred3, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "CN & TMWB models") + theme_classic() + xlab(NULL) +
  scale_y_continuous(name = "Depth (mm)") + ggtitle("Boulder Creek")
p13 <- ggplot(CompareQpred4, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "CN & TMWB models") + theme_classic() + xlab(NULL) +
  scale_y_continuous(name = "Depth (mm)") 
pdf("HW3 CNmodel.pdf")
p12/p13
dev.off()

pdf("HW3 2models.pdf")
p5/p12 
dev.off()

###############################################################################

#####################               HW 4               #######################

#################################################################################
myflowgage_id = "01421618"
myflowgage = get_usgs_gage(myflowgage_id, begin_date = "2015-01-01", end_date = "2022-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

stns = meteo_distance(
  station_data = ghcnd_stations(),
  lat = myflowgage$declat,
  long = myflowgage$declon,
  units = "deg",
  radius = 20,
  limit = NULL
)
# We are looking for stations with elements that have PRCP, TMAX and TMIN 
# and current data (i.e. Year 2021). 
WXStn = stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData = meteo_pull_monitors(
  monitors = WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)

# Create an aligned modeldata data frame to build our model in
modeldata = merge(WXData, myflowgage$flowdata, by.x = "date", by.y = "mdate")
summary(modeldata)  #
plot(modeldata$date, modeldata$Qmm)
lines(modeldata$date, modeldata$prcp/10)
modeldata$MaxTemp = modeldata$tmax/10 # Converting to C
modeldata$MinTemp = modeldata$tmin/10 # Converting to C
modeldata$P = modeldata$prcp/10 # Converting to mm
View(modeldata)  
# Compare your precipitation to the flow out of your basin
modeldata$P[is.na(modeldata$P)] = 0
modeldata$MinTemp[is.na(modeldata$MinTemp)] = 0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)] =
  modeldata$MinTemp[is.na(modeldata$MaxTemp)] + 1
modeldata$MaxTemp[modeldata$MaxTemp <= modeldata$MinTemp] =
  modeldata$MinTemp[modeldata$MaxTemp <= modeldata$MinTemp] + 1
modeldata$AvgTemp = (modeldata$MaxTemp + modeldata$MinTemp)/2.0

summary(modeldata)
modeldata[is.na(modeldata)] = 0 # A Quick BUT sloppy removal of NAs
modeldata <- modeldata %>% filter(date <= as.Date("2020-02-29") & date >= as.Date("2016-09-30"))

####################################Soil and Topo data#####################
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/mc215itubhtxlvmeoietapfh/wss_aoi_2022-02-23_19-07-26.zip"
download.file(url,"mysoil.zip")
unzip("mysoil.zip")
mysoil = readOGR("wss_aoi_2022-02-23_19-07-26/spatial/soilmu_a_aoi.shp")
# Explore the mysoil dataset which is returned

mybbox = c(mysoil@bbox)
# First associate mukey with cokey from component
mysoil$mukey = mysoil$MUKEY  # or rename the column
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
print(mukey_statement)
q_mu2co = paste("SELECT mukey, cokey FROM component WHERE mukey IN ", mukey_statement, sep = "")
print(q_mu2co)
mu2co = SDA_query(q_mu2co)
# Second associate cokey with ksat_r,awc_r,hzdepb_r from chorizon
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))
q_co2ch = paste("SELECT cokey,ksat_r,awc_r,hzdepb_r  FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
print(q_co2ch)
co2ch = SDA_query(q_co2ch)
# Last, bring them back together, and aggregate based on max values
# of ksat_r,awc_r, and hzdepb_r

mu2ch = merge(mu2co,co2ch)
summary(mu2ch)
mu2chmax = aggregate(mu2ch,list(mu2ch$mukey), max)

proj4_ll = "+proj=longlat"
proj4string(mysoil) = proj4_ll
mydem = get_elev_raster(locations = mysoil, 
                        z = 11, prj = proj4string(mysoil) ,
                        src ="aws",clip="bbox",expand = 0.001)

summary(terrain(mydem, opt='TRI',unit = "radians"))
plot(terrain(mydem, opt='TRI',unit = "radians"))

####################################CNModel#################################

TownBrookCN = modeldata

##################Gather soil and topo parameters###########

##topslope has the lowest Z, lowest slope as botslope, midslope has the #highest slope, we assume that awc in % is 

awcpercent=mean(mu2chmax$awc_r,na.rm=TRUE)
TownBrookZ=mean(mu2chmax$hzdepb_r)*10 #in mm
AWCvalTop=TownBrookZ*awcpercent
slop_stat = summary(terrain(mydem, opt='Slope', unit = "radians"))
slop_stat
SlopeMedian = slop_stat #radians
CNavg = c(50, 55, 60, 65, 70)
fnc_fcres= c(0.1, 0.2, 0.3, 0.4, 0.5)

# Call the new CNmodel() function with 1 HRU objects,
TownBrookCN = CNmodel(CNmodeldf = TownBrookCN, IaFrac = 0.05, CNavg = CNavg[1], fnc_slope = SlopeMedian[[1]],
                      fnc_aspect = 0, func_DAWC=awcpercent,
                      func_z = TownBrookZ, fnc_fcres = fnc_fcres[1])



plotAWCN = data.frame(TownBrookCN$date, TownBrookCN$Qmm, TownBrookCN$Qpred)
colnames(plotAWCN) <- c("date", "CNQmm", "CNQpred")
plotAWCN = plotAWCN  %>% pivot_longer(-date)
p15 <- ggplot(plotAWCN, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Available Water") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 


NSeff(TownBrookCN$Qmm,TownBrookCN$Qpred)












#FOR Excess
TownBrookCN = mutate(TownBrookCN, TopEx = TownBrookCN$Excess)
MidSlopeCN2 = mutate(MidSlopeCN2, MidEx = MidSlopeCN2$Excess)
BotSlopeCN2 = mutate(BotSlopeCN2, BotEx = BotSlopeCN2$Excess)
plotExCN = left_join(TownBrookCN[, c("date","TopEx")], 
                     MidSlopeCN2[, c("date","MidEx")], by = "date") %>% 
  left_join(BotSlopeCN2[, c("date","BotEx")]) %>% pivot_longer(-date)
p8 <- ggplot(plotExCN, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Excess") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

#pdf("awex2.pdf")
p7/p8
#dev.off()


# CN Model Performance 

plot(BotSlopeCN2$date,cumsum(BotSlopeCN2$Excess),type="l", col = "red")


# finish building all the hillslope HRUs….

####################################TMWBM#################################
TopSlope2 = modeldata


TopSlope2 = TMWBmodel(TMWB = TopSlope2,SFTmp = 1, 
                      AWCval = AWCvalTop,
                      Tlag = .5,fcres=.3,Slope = SlopeTop)

##############AW Plots
TopSlope2 = mutate(TopSlope2, TopAW = TopSlope2$AW)

plotAW = left_join(TopSlope2[, c("date","TopAW")], 
                   MidSlope2[, c("date","MidAW")], by = "date") %>% 
  left_join(BotSlope2[, c("date","BotAW")]) %>% pivot_longer(-date)
p9 <- ggplot(plotAW, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Available Water") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

#FOR Excess
TopSlope2 = mutate(TopSlope2, TopEx = TopSlope2$Excess)
MidSlope2 = mutate(MidSlope2, MidEx = MidSlope2$Excess)
BotSlope2 = mutate(BotSlope2, BotEx = BotSlope2$Excess)
plotEx = left_join(TopSlope2[, c("date","TopEx")], 
                   MidSlope2[, c("date","MidEx")], by = "date") %>% 
  left_join(BotSlope2[, c("date","BotEx")]) %>% pivot_longer(-date)
p10 <- ggplot(plotEx, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Excess") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 

p9/p10

# TMWB Model Performance 
plot(BotSlope2$date,BotSlope2$Qpred,type="l")
lines(BotSlope2CN2$date, BotSlope2CN2$Qpred, type="l", col = "red")
NSeff(BotSlope2$Qmm,BotSlope2$Qpred)

# Plot Qpred, Qmm for CN and TMWB
CompareQpred = data.frame(modeldata$date, modeldata$Qmm, BotSlope2$Qpred, BotSlopeCN2$Qpred)
colnames(CompareQpred) <- c("date", "Qmm", "TMWBQpred", "CNQpred")
CompareQpred = CompareQpred %>% pivot_longer(-date)
p11 <- ggplot(CompareQpred, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "CN & TMWB models") + theme_classic() + xlab(NULL) +
  scale_y_continuous(name = "Depth (mm)") 

pdf("HW3TMWB.pdf")
p11
dev.off()

CompareQpred3 = data.frame(BotSlopeCN2$date, BotSlopeCN2$Qmm, BotSlope2$Qpred, BotSlopeCN2$Qpred)
colnames(CompareQpred3) <- c("date", "Qmm", "TMWBQpred","CNQpred")
CompareQpred3 = CompareQpred3 %>% pivot_longer(-date)
CompareQpred4 = data.frame(BotSlopeCN2$date, BotSlopeCN2$SnowfallWatEq_mm, BotSlopeCN2$SnowMelt_mm)
colnames(CompareQpred4) <- c("date", "CNSnow Water Equ.", "CNSnow Melt")
CompareQpred4 = CompareQpred4 %>% pivot_longer(-date)
p12 <- ggplot(CompareQpred3, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "CN & TMWB models") + theme_classic() + xlab(NULL) +
  scale_y_continuous(name = "Depth (mm)") + ggtitle("Boulder Creek")
p13 <- ggplot(CompareQpred4, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "CN & TMWB models") + theme_classic() + xlab(NULL) +
  scale_y_continuous(name = "Depth (mm)") 
pdf("HW3 CNmodel.pdf")
p12/p13
dev.off()

pdf("HW3 2models.pdf")
p5/p12 
dev.off()
