options(repos ="http://cran.us.r-project.org")  # required to get latest libs
# Installing the packages we will play with today
if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,soilDB,rgdal,raster, rgeos, tidyverse, patchwork, aqp)
# From previous weeks
pacman::p_load(EcoHydRology,rnoaa,curl,httr)
rm(list=objects())
setwd("C:/Users/binyam/OneDrive - Virginia Tech/Advanced Watershed Modelling/Week 4/Homework/Lab04wk04")
#dir.create("~/Lab04")
#setwd("~/Lab04/")
# 3 Functions to calculate SWE and excess when soil is drying, 
#   wetting, and wetting above capacity

browseURL("https://github.com/vtdrfuka/BSE5304_2022/tree/main/functions")
browseURL("https://github.com/vtdrfuka/BSE5304_2022/blob/main/functions/TMWBmodel.R")
browseURL("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/NSE.R")
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")

# Download a soils dataset for your basin based on the WebSoilSurvey method 
# and replace this url with your own
url="https://websoilsurvey.sc.egov.usda.gov/DSD/Download/AOI/tosambsrx5wdif00sp153f22/wss_aoi_2022-02-17_13-53-03.zip"
download.file(url, "mysoil.zip")
unzip("mysoil.zip")
# https://cran.r-project.org/web/packages/elevatr/elevatr.pdf
# https://cran.r-project.org/web/packages/raster/raster.pdf
# https://cran.r-project.org/web/packages/soilDB/soilDB.pdf
# https://cran.r-project.org/web/packages/rgdal/rgdal.pdf
# use the function to get data from USGS 12115700
# BOULDER CREEK NEAR CEDAR FALLS, WA
myflowgage_id = "12115700"
myflowgage = get_usgs_gage(myflowgage_id, begin_date = "2015-01-01", end_date = "2022-03-01")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

# Read the soil data using the soilDB package
# What is this returning? Why do we care?
######### It reads in an OGR data and vector layer from USDA-NCSS soil databases
#into a suitable spatial vector in this case a shape file is returned in a 
#dataframe of spatial polygons

# This needs to be completed based on your download

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

summary(terrain(mydem, opt = 'slope', unit = "degrees"))
# What is this 'slope'? Use the man page for the terrain() function to answer
# terrain slope in degrees computed from DEM raster elevation values

plot(terrain(mydem, opt='TPI', unit = "degrees"))
# What is this 'TPI'? 
#Topographic Position Index, a method of terrain classification which 
#evaluate a point/raster altitude against its surrounding. + values for higher position
# - values for lower position

summary(terrain(mydem, opt='TRI',unit = "degrees"))
plot(terrain(mydem, opt='TRI',unit = "degrees"))
# What is this 'TRI'? 
# It is Terrain Ruggedness Index, it provides the mean difference between a 
# pixel/raster cell and its surroundings

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
WXStn=stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[1]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
summary(WXData)  #

# Create an aligned modeldata data frame to build our model in
modeldata = merge(WXData, myflowgage$flowdata, by.x = "date", by.y = "mdate")
summary(modeldata)  #
modeldata$MaxTemp = modeldata$tmax/10 # Converting to C
modeldata$MinTemp = modeldata$tmin/10 # Converting to C
modeldata$P = modeldata$prcp/10 # Converting to mm
# View(modeldata)  
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
# Last weeks homework example
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/TMWBmodel.R")
# Calibrating the parameters one at a time
for (fcres in seq(.1,.5,.1)){
  TMWBnew = TMWBmodel(TMWB = TMWB, fcres = fcres)
  print(paste(fcres, NSE(TMWBnew$Qmm, TMWBnew$Qpred)))
}
# fcres = 0.1 is the highest NSE = 0.187
for (SFTmp in seq(-5,20)){
  TMWBnew = TMWBmodel(TMWB = TMWB, fcres = .1, SFTmp = SFTmp)
  print(paste(SFTmp, NSE(TMWBnew$Qmm, TMWBnew$Qpred)))
}
# SFTmp = 3 is the highest NSE = 0.268
for(AWCval in seq(50,350,50)){   #updated AWC value ranges based on soil data set
  TMWBnew = TMWBmodel(TMWB = TMWB, fcres = .1, SFTmp = 3, Tlag = .5, AWCval = AWCval)
  print(paste(AWCval, NSE(TMWBnew$Qmm, TMWBnew$Qpred)))
}

# Best result for "BOULDER CREEK NEAR CEDAR FALLS, WA" NSE = 0.277. If processed as Lab03
TMWBnew = TMWBmodel(TMWB = TMWB, fcres = .1, SFTmp = 3, Tlag = .5, AWCval = 100)
print(paste(AWCval, NSE(TMWBnew$Qmm, TMWBnew$Qpred))) # Error in code AWCval 
#gives wrong value in the current scope (outside the function)

# For a simple spatial model, use TMWB, to initialize 3 
# slope components for Top, Mid, and Bottom of the hillside. 
# First compute Excess for Top area using estimated + from calibration above
#area specific change for Slope, fcres, AWC fn(Soil Depth & %AWC)

TopSlope = TMWB
AWCtop = mean(mu2chmax$hzdepb_r, na.rm = T) * mean(mu2chmax$awc_r, na.rm = T)*10 #convert to mm
slop_stat = summary(terrain(mydem, opt='Slope', unit = "radians"))
Slopetop = slop_stat[3]   # index value for Median slope in slop_stat
fcres =  0.1            #Average is 0.1 for basin, best value above
TopSlope = TMWBmodel(TMWB=TopSlope,fcres=fcres,SFTmp=3,bmlt6=2.5,bmlt12=1,Tlag=.5,
                     AWCval=AWCtop,Slope=Slopetop)
print(NSE(TopSlope$Qmm, TopSlope$Qpred)) # NSE comparison is not valid here as 
#the areas are supposed to generate different excess and Qpred is computed after BotSlope

#Compute MidSlope with updated TMWB Data
MidSlope = TMWB
mean(MidSlope$P)
MidSlope$P = MidSlope$P + TopSlope$Excess
mean(MidSlope$P)   # Quick check for change P

AWCmid = min(mu2chmax$hzdepb_r, na.rm = T) * min(mu2chmax$awc_r, na.rm = T)*10 #convert to mm
AWCmid
slop_stat = summary(terrain(mydem, opt='Slope', unit = "radians"))
slop_stat
Slopemid = slop_stat[4]   # index value for 3rd Quant slope in slop_stat
Slopemid
fcres =  0.1            #Average is 0.1 for basin, best value above
MidSlope = TMWBmodel(TMWB=MidSlope,fcres=fcres,SFTmp=3,bmlt6=2.5,bmlt12=1,Tlag=.5,
                     AWCval=AWCmid,Slope=Slopemid)
print(NSE(MidSlope$Qmm, MidSlope$Qpred)) # NSE comparison is not valid here as 
#the areas are supposed to generate different excess and Qpred is computed after BotSlope

#Compute BotSlope with updated TMWB Data
BotSlope = TMWB
mean(BotSlope$P)
BotSlope$P = BotSlope$P + MidSlope$Excess
mean(BotSlope$P)   # Quick check for change P

AWCbot = max(mu2chmax$hzdepb_r, na.rm = T) * max(mu2chmax$awc_r, na.rm = T)*10 #convert to mm
AWCbot
slop_stat = summary(terrain(mydem, opt='Slope', unit = "radians"))
slop_stat
Slopebot = slop_stat[2]   # index value for 1st Quant slope in slop_stat
Slopebot
fcres =  0.01            #Average is 0.1 for basin, best value above
BotSlope = TMWBmodel(TMWB=BotSlope,fcres=fcres,SFTmp=3,bmlt6=2.5,bmlt12=1,Tlag=.5,
                     AWCval=AWCbot,Slope=Slopebot)
print(NSE(BotSlope$Qmm, BotSlope$Qpred)) # NSE comparison is not valid here as 
#the areas are supposed to generate different excess and Qpred is computed after BotSlope
# Homework 1
# For AW 
TopSlope = mutate(TopSlope, TopAW = TopSlope$AW)
MidSlope = mutate(MidSlope, MidAW = MidSlope$AW)
BotSlope = mutate(BotSlope, BotAW = BotSlope$AW)
plotAW = left_join(TopSlope[, c("date","TopAW")], 
                   MidSlope[, c("date","MidAW")], by = "date") %>% 
  left_join(BotSlope[, c("date","BotAW")]) %>% pivot_longer(-date)
p1 <- ggplot(plotAW, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Available Water") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 
  
#FOR Excess
TopSlope = mutate(TopSlope, TopEx = TopSlope$Excess)
MidSlope = mutate(MidSlope, MidEx = MidSlope$Excess)
BotSlope = mutate(BotSlope, BotEx = BotSlope$Excess)
plotEx = left_join(TopSlope[, c("date","TopEx")], 
                   MidSlope[, c("date","MidEx")], by = "date") %>% 
  left_join(BotSlope[, c("date","BotEx")]) %>% pivot_longer(-date)
p2 <- ggplot(plotEx, aes(x = date, y = value, color = name)) +
  geom_line() + labs(color = "Excess") + theme_classic() +
  scale_y_continuous(name = "Depth (mm)") 
  
pdf("awex2.pdf")
p1/p2
dev.off()

#Homework 2
mybbox = c(mysoil@bbox)
BoulderCrk_SP <- fetchKSSL(bbox = c(-122,47,-121,48)) #I have widened the search 
#box as there were no pedons in the creek watershed

clay <- groupedProfilePlot(BoulderCrk_SP, groups='pedon_key', group.name.cex=0.65, color='clay', 
                   name='hzn_desgn', id.style='side', label='pedon_id', max.depth=200, group.line.lwd =1)
mtext('pedon_key', side=3, line=-1, at = c(0,0), adj = 0)

silt <- groupedProfilePlot(BoulderCrk_SP, groups='pedon_key', group.name.cex=0.65, color='silt', 
                           name='hzn_desgn', id.style='side', label='pedon_id', max.depth=200, group.line.lwd =1)
mtext('pedon_key', side=3, line=-1, at = c(0,0), adj = 0)

sand <- groupedProfilePlot(BoulderCrk_SP, groups='pedon_key', group.name.cex=0.65, color='sand', 
                           name='hzn_desgn', id.style='side', label='pedon_id', max.depth=200, group.line.lwd =1)
mtext('pedon_key', side=3, line=-1, at = c(0,0), adj = 0)

clay + silt + sand
dev.off()





