if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology, rnoaa, elevatr,soilDB,rgdal,raster, rgeos, tidyverse, patchwork, aqp)

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
TMWB = modeldata

#CNModel
#
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
#
# Like before, initializing the 3 hillslope classes
#
TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata
# Call the new CNmodel() function with Top,Mid,BotSlope HRU objects,
# passing the Qpred into the lower HRUs HillslopeAboveExcess (as area scaled flow)
TopSlopeCN=CNmodel(TopSlopeCN, CNavg = 60)
TopSlopeCN = CNmodel(CNmodeldf = TopSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=500,fnc_fcres=.3)
MidSlope$P=TopSlope$Excess+MidSlope$P
# Higher slope, medium ksat, fcres=0.5 
MidSlopeCN = CNmodel(CNmodeldf = MidSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$P=MidSlope$Excess+BotSlope$P
BotSlopeCN = CNmodel(CNmodeldf = BotSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=1000,fnc_fcres=.2)


#AW Plots HW1
plot(BotSlopeCN$date,BotSlopeCN$AW,type="l",col=1,xlab="Date",ylab="AW (mm)")
lines(MidSlopeCN$date,MidSlopeCN$AW,type="l",col=2)
lines(TopSlopeCN$date,TopSlopeCN$AW,type="l",col=3)
# Excess Plots HW1
plot(BotSlopeCN$date,BotSlopeCN$Excess,type="l",col=1,xlab="Date",ylab="Excess (mm)")
lines(MidSlopeCN$date,MidSlopeCN$Excess,type="l",col=2)
lines(TopSlopeCN$date,TopSlopeCN$Excess,type="l",col=3)

# PET and ET HW2
plot(BotSlopeCN$date,BotSlopeCN$PET,type="l",col=1,xlab="Date",ylab="(P)ET (mm)")
lines(BotSlopeCN$date,BotSlopeCN$ET,type="l",col=2)
lines(MidSlopeCN$date,MidSlopeCN$ET,type="l",col=3)
lines(TopSlopeCN$date,TopSlopeCN$ET,type="l",col=4)
# or as cumulative summations
plot(TopSlopeCN$date,cumsum(BotSlopeCN$PET),type="l",
     xlab="Date",ylab="(P)ET")
lines(TopSlopeCN$date,cumsum(TopSlopeCN$ET),col="red")
lines(MidSlopeCN$date,cumsum(MidSlopeCN$ET),col="green")
lines(BotSlopeCN$date,cumsum(BotSlopeCN$ET),col="blue")


# Cumulative Summary of QPred is very informative
plot(BotSlopeCN$date,cumsum(BotSlopeCN$Qpred),type="l",
     xlab="Date",ylab="Flow Q Cumulative Summary (mm)")
lines(MidSlopeCN$date,cumsum(MidSlopeCN$Qpred),col="red")
lines(TopSlopeCN$date,cumsum(TopSlopeCN$Qpred),col="green")

# Model Performance 
plot(BotSlopeCN$date,BotSlopeCN$Qpred,type="l")
NSeff(BotSlopeCN$Qmm,BotSlopeCN$Qpred)

# finish building all the hillslope HRUs….
