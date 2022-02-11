# Binyam Asfaw
# Solution to Lab 03 HW 03
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rnoaa,EcoHydRology,lattice, dataRetrieval, tidyverse, lubridate, patchwork)

# For a station near Blacksburg we look at
# https://maps.waterdata.usgs.gov/mapper/index.html
# or https://waterdata.usgs.gov/nwis/rt
# and like gage:

#find sites in a bounding box coords of bottom left, top right
SawNF <- c(-115.61, 43.39, -114.56, 44.36)

#get sites in this bounding box that have daily discharge
SawNF_sites <- whatNWISsites(bBox = SawNF, parameterCD = c("00060", "00010"), hasDataTypeCd = "dv", drainAreaMax = "200", Area = drainAreaMax)

start <- "2016-01-01"
end <- "2022-01-31"
params <- c("00060")

SawNF_dat <- readNWISdv(siteNumbers = SawNF_sites$site_no, parameterCd = params, startDate = start, endDate = end) %>%
  renameNWISColumns()

usgsid <- SawNF_dat$site_no[1]


# Finding weather station with names:waterdata.usgs.gov
myflowgage_id = usgsid
myflowgage=get_usgs_gage(myflowgage_id,
                         begin_date="2016-01-01",end_date="2022-01-30")

# Now use the functions meteo_distance and ghcnd_stations to start looking for weather stations nearby
 
stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)
# We are looking for stations with elements that have PRCP, TMAX and TMIN and current data (i.e. Year 2021). 
WXStn = stns[stns$element=="TMAX"&stns$last_year>=2021,]$id[4]
WXData = meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP")
)
summary(WXData)  #
plot(WXData$date, WXData$tmax, col = "black")
points(WXData$date, WXData$tmin, col = "red")
points(WXData$date, WXData$prcp, col = "blue")

# Creat an aligned modeldata data frame to build our model in
modeldata = merge(WXData, myflowgage$flowdata, by.x="date", by.y="mdate")
summary(modeldata)

# Convert Units and Area normalize flow to match (depth)
# flow(m^3/day) / (area(km^2) * 1000^2m/km^2) * 1000mm/m = flow(mm/day)

modeldata$Qmm = modeldata$flow/myflowgage$area/10^3

# It is good practice to use similar object names to the values in the documentation of the model 
#(P, Q, MaxTemp, MinTemp) 

modeldata$MaxTemp = modeldata$tmax/10 # Converting to C
modeldata$MinTemp = modeldata$tmin/10 # Converting to C
modeldata$P = modeldata$prcp/10 # Converting to mm

View(modeldata)  
# Compare your precipitation to the flow out of your basin 

mean(modeldata$Qmm)
mean(modeldata$P, na.rm = T)

modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]= modeldata$MinTemp[is.na(modeldata$MaxTemp)] + 1

summary(modeldata)

TMWB = modeldata

summary(TMWB)

# Building our Soil Wetting and Drying Functions

soilwetting <- function(AWprev,dP_func,AWC_func){
  AW_func <- AWprev+dP_func
  excess_func <- 0.0
  c(AW_func,excess_func)
} 

soildrying <- function(AWprev,dP_func,AWC_func){
  AW_func = AWprev*exp(dP_func/AWC_func)
  excess_func <- 0.0
  c(AW_func,excess_func)
}
# soil_wetting_above_capacity function
soil_wetting_above_capacity<-function(AWprev,dP_func,AWC_func){
  AW_func<-AWC_func
  excess_func<-AWprev+dP_func-AWC_func
  c(AW_func,excess_func)
}

# Add some some soil parameters to be associated with the area 
# above the flow gage.

myflowgage$FldCap = 0.45
myflowgage$WiltPt = 0.15
myflowgage$Z = 1000
TMWB$AWC = (myflowgage$FldCap - myflowgage$WiltPt)*myflowgage$Z # 


TMWB$PET = mean(TMWB$P,na.rm=T)-mean(TMWB$Qmm,na.rm=T)  # in mm/day
TMWB$ET = TMWB$PET # in mm/day
TMWB$dP = TMWB$P - TMWB$PET
TMWB$AW = NA  #Assigns all values in column with “NA” (Not available)
TMWB$AW[1] = 250
TMWB$Excess = NA
TMWB$Excess[1] = 0
head(TMWB)
attach(TMWB)
for (t in 2:length(date)){
  if (dP[t]< 0) {  
    values<-soildrying(AW[t-1],dP[t],AWC[t])
  } else if (AW[t-1]+dP[t]>AWC[t]) {
    values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
  } else {
    values<-soilwetting (AW[t-1],dP[t],AWC[t])
  }
  AW[t]<-values[1]
  Excess[t]<-values[2]
  
}
detach(TMWB)
TMWB$AW <-AW
TMWB$Excess<-Excess
rm(list=c("AW","Excess"))
TMWB$Qpred=NA
TMWB$Qpred[1]=0
TMWB$S=NA
TMWB$S[1]=0

attach(TMWB)
fcres=.3   # reservoir coefficient
for (t in 2:length(date)){
  S[t]=S[t-1]+Excess[t]     
  Qpred[t]=fcres*S[t]
  S[t]=S[t]-Qpred[t]
}
detach(TMWB) # IMPORTANT TO DETACH
TMWB$S=S
TMWB$Qpred=Qpred # UPDATE vector BEFORE DETACHING
rm(list=c("S","Qpred"))
summary(TMWB)
View(TMWB)
dev.off()
####################################
###let's start to plot
#####################################
#pcp and Qmm and Qpred
maxRange <- 1.1*(max(TMWB$P,na.rm = T) + max(TMWB$Qmm,na.rm = T))
p1<- ggplot() +
  # Use geom_tile to create the inverted hyetograph. geom_tile has a bug that displays a warning message for height and width, you can ignore it.
  geom_tile(data = TMWB, aes(x=date,y = -1*(P/2-maxRange), # y = the center point of each bar
                             height = P,
                             width = 1),
            fill = "black",
            color = "black") +
  # Plot your discharge data
  geom_line(data=TMWB,aes(x=date, y = Qmm, colour ="Qmm"), size=1) +
  geom_line(data=TMWB,aes(x=date, y = Qpred, colour= "Qpred"), size=1) +
  scale_colour_manual("", 
                      breaks = c("Qmm", "Qpred"),
                      values = c("red", "blue")) +
  # Create a second axis with sec_axis() and format the labels to display the original precipitation units.
  scale_y_continuous(name = "Discharge (mm/day)",
                     sec.axis = sec_axis(trans = ~-1*(.-maxRange),
                                         name = "Precipitation (mm/day)"))+
  scale_x_continuous(name = NULL,labels = NULL)+
  ggtitle(myflowgage$gagename)

#ET and Excess
p2 <- ggplot(TMWB, aes(x=date)) +
  geom_line(aes(y=Excess, colour="Excess"), size=1) +
  geom_line(aes(y=ET, colour="ET"), size=1) +
  scale_colour_manual("", 
                      breaks = c("ET", "Excess"),
                      values = c("red", "blue")) +
  scale_y_continuous(name = "Depth (mm/day)",) +
  scale_x_continuous(name = NULL,labels = NULL) 
#FOR AW
p3 <- ggplot(TMWB, aes(x=date)) +
  geom_line(aes(y=AW,colour="AW"), size=1) +
  scale_colour_manual("", 
                      breaks = c("AW"),
                      values = c("black")) +
  scale_y_continuous(
    # Features of the first axis
    name = "AW (mm)",
    
  )

p1 + p2 + p3 + plot_layout(ncol = 1, widths = c(2,2,1))

####################################
#Lab03wk03

#notice that there is an Energy Balance based Snow Accumulation 
#and Melt model in the EcoHydRology package.
############################################
### Problem 1
############################################
attach(TMWB)
slope = c(0, atan(10/100), atan(10/100), atan(45/100), atan(45/100))
aspect = c(0, 0, 180*(pi/180), (360-45)*(pi/180), (360-135)*(pi/180))
Output <- list()
for (i in 1:length(slope)) {
  Output[[i]] <- SnowMelt(date, P, MaxTemp, MinTemp, myflowgage$declat, 
                          slope = slope[i],
                          aspect = aspect[i], tempHt = 1, windHt = 2, groundAlbedo = 0.25,
                          SurfEmissiv = 0.95, windSp = 2, forest = 0, startingSnowDepth_m = 0,
                          startingSnowDensity_kg_m3=450)
  
}
detach(TMWB)
SNO_Energy_1 <- Output[[1]]
SNO_Energy_2 <- Output[[2]]
SNO_Energy_3 <- Output[[3]]
SNO_Energy_4 <- Output[[4]]
SNO_Energy_5 <- Output[[5]]
SNO_Energy_1_5 <- SNO_Energy_1[, c("Date", "SnowWaterEq_mm")]
SNO_Energy_1_5 <- SNO_Energy_1_5 %>% mutate("SnowWaterEq_mm_2" = SNO_Energy_2$SnowWaterEq_mm, 
                                            "SnowWaterEq_mm_3" = SNO_Energy_3$SnowWaterEq_mm,
                                            "SnowWaterEq_mm_4" = SNO_Energy_4$SnowWaterEq_mm,
                                            "SnowWaterEq_mm_5" = SNO_Energy_5$SnowWaterEq_mm)
summary(SNO_Energy_1_5)
#Plot Energy Balance Based SWE
pdf("SnoWaEq.pdf")
SNO_P_1 <- ggplot(SNO_Energy_1_5, aes(x = Date)) +
  geom_line(aes(y = SnowWaterEq_mm, color = "a")) +
  geom_line(aes(y = SnowWaterEq_mm_2, color = "b")) +
  geom_line(aes(y = SnowWaterEq_mm_3, color = "c")) +
  geom_line(aes(y = SnowWaterEq_mm_4, color = "d")) +
  geom_line(aes(y = SnowWaterEq_mm_5, color = "e")) +
  theme_classic() +
  labs(color = "Slope, aspect combination") +
  ggtitle("Effect of slope and aspect on snow accumulation")
SNO_P_2 <- SNO_P_1 + scale_x_date(name = "2017_Jan-Jun", limits = as.Date(c("2017-01-09", "2017-06-30")))
SNO_P_1 + SNO_P_2 + plot_layout(ncol = 1, nrow = 2) 
dev.off()
#############################################
#Problem 2
####################################

SFTmp = -5    # referred to as SFTMP in SWAT input (Table 1) snow fall temp.
bmlt6 = 7   # referred to as SMFMX in SWAT input (Table 1) melt factor for snow June21
bmlt12 = 7  # referred to as SMFMN in SWAT input adjusted for season. Minimum melt rate for snow during years 
Tmlt = SFTmp  # assumed to be same as Snow Fall Temperature
Tlag = 0.2    # referred to as TIMP in SWAT input (Table 1)
TMWB$AvgTemp <- (TMWB$MaxTemp + TMWB$MinTemp)/2
TMWB$bmlt = (bmlt6 + bmlt12)/2 + (bmlt6 - bmlt12)/2 *  sin(2*pi/365*(julian(TMWB$date,origin = as.Date("2000-01-01"))-81))

# Initialize SNO, Tsno as well as the first values of each

TMWB$SNO = 0  # Snow Depth (mm)
TMWB$Tsno = 0  # Snow Temp (C)
TMWB$SNOmlt = 0  # Snow Melt (mm)

attach(TMWB)
for (t in 2:length(date)){
  Tsno[t]= Tsno[t-1] * (1.0-Tlag) +  AvgTemp[t] * Tlag
  if(AvgTemp[t] < SFTmp){
    SNO[t]= SNO[t-1] + P[t]
  }  else {
    SNOmlt[t]= bmlt[t] * SNO[t-1] * ((Tsno[t] + MaxTemp[t])/2 - Tmlt) 
    SNOmlt[t]= min(SNOmlt[t],SNO[t-1])
    SNO[t]= SNO[t-1] - SNOmlt[t]
  }
  print(t)
}

detach(TMWB)

TMWB$Tsno = Tsno
TMWB$SNO = SNO
TMWB$SNOmlt = SNOmlt
rm(list = c("SNO", "SNOmlt", "Tsno"))


############################################
TMWB$Albedo=.23
TMWB$Albedo[TMWB$SNO>0]=.95

attach(TMWB)
PET = PET_fromTemp(Jday = (1+as.POSIXlt(date)$yday), Tmax_C = MaxTemp, Tmin_C = MinTemp,
                   albedo=Albedo,lat_radians = myflowgage$declat*pi/180) * 1000
TMWB$PET=PET
plot(date,PET)
detach(TMWB)
rm(list=c("PET"))

#############################################

myflowgage$FldCap=.45
myflowgage$WiltPt=.15
myflowgage$Z=1000
TMWB$AWC = 350 #(myflowgage$FldCap-myflowgage$WiltPt)*myflowgage$Z # 
TMWB$dP = 0 # Initializing Net Precipitation
TMWB$ET = 0 # Initializing ET
TMWB$AW = 0 # Initializing AW
TMWB$Excess = 0 # Initializing Excess

# Loop to calculate AW and Excess

attach(TMWB)
for (t in 2:length(AW)){
  # This is where Net Precipitation is now calculated

  ET[t] = min (AW[t-1], PET[t])
  ET[t] = (AW[t-1]/AWC[t-1])*PET[t] # New Model
  if(AvgTemp[t] >= SFTmp){
    dP[t] = P[t] - ET[t] + SNOmlt[t] 
  }  else {
    dP[t] = ET[t]
  }
  # compute soil water balance
  if (dP[t]<=0) {
    values <- soildrying(AW[t-1], dP[t], AWC[t])
  } else if((dP[t]>0) & (AW[t-1] + dP[t]) <= AWC[t]) {
    values <- soilwetting(AW[t-1], dP[t], AWC[t])
  } else {
    values <- soil_wetting_above_capacity(AW[t-1], dP[t], AWC[t])
  }
  AW[t]<-values[1]
  Excess[t]<-values[2]
  print(t)
}
TMWB$AW = AW
TMWB$Excess = Excess
TMWB$dP = dP
rm(list=c("AW","dP","ET", "Excess"))
detach(TMWB) # IMPORTANT TO DETACH

#####################################
TMWB$Qpred = NA
TMWB$Qpred[1] = 0
TMWB$S = NA
TMWB$S[1] = 0
attach(TMWB)
fcres = 0.18
for (t in 2:length(date)){
                        S[t] = S[t-1] + Excess[t]     
                        Qpred[t] = fcres*S[t]
                        S[t] = S[t] - Qpred[t]
}
TMWB$S = S
TMWB$Qpred = Qpred # UPDATE vector BEFORE DETACHING
plot(date,P,col="black")
lines(date,Qmm,type = "l",col="black")
lines(date,Qpred,col="blue")

detach(TMWB) # IMPORTANT TO DETACH
rm(list=c("Qpred","S"))

####################################
#collect temporary Qpreds
####################################
Qpred_df <- data.frame(date = TMWB$date)
Qpred_df <- mutate(Qpred_df, cal_6 = TMWB$Qpred)

#############################
#NSE for reservoir coefficients
NSE = function(Yobs,Ysim){
  return(1-sum((Yobs-Ysim)^2, na.rm=TRUE)/sum((Yobs-mean(Yobs, na.rm=TRUE))^2, na.rm=TRUE))
}
NSE_AWC <- c()
#convert_data = dataframe_name[[‘column_name’]]
Yobs <- TMWB[["Qmm"]]
Ysimu <- c("cal_1","cal_2","cal_3","cal_4","cal_5", "cal_6")
for (i in 1:length(Ysimu)) {
  Ysim <- Qpred_df[[Ysimu[i]]]
  value <- NSE(Yobs, Ysim)
  NSE_AWC <- append(NSE_AWC, value)
}
NSE_SFTmp
NSE_fcres
NSE_bmlt6
NSE_bmlt12
NSE_tlag
NSE_AWC
# Qmm and Qpred Plot
pdf("Simulated flow.pdf")
ggplot(TMWB, aes(x=date)) +
  geom_line(aes(y=Qmm, colour="Qmm"), size=1) +
  geom_line(aes(y=Qpred, colour="Qpred"), size=1) +
  scale_y_continuous(name = "Depth (mm/day)",) +
  ggtitle("Observed and simulated Discharge in mm")
dev.off()
