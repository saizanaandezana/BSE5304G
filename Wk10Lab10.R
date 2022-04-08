objects()
rm(list=objects())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,elevatr,raster,rgdal,
                 data.table,foreign,maptools,dataRetrieval,gdistance)
setwd("~/Wk10Lab10/")
#
# Note we have a new library to access USGS Waterdata
# https://owi.usgs.gov/R/dataRetrieval.html
# https://owi.usgs.gov/R/training-curriculum/usgs-packages/dataRetrieval-readNWIS/
#
sites <- whatNWISsites(bBox=c(-113.581, 35.743,-111.680,36.528),  # search for USGS sites using Bbox
                       parameterCd=c("00065","00060"),
                       hasDataTypeCd="uv")
#
# Need to figure out which data to download. 
# https://nwis.waterdata.usgs.gov/nwis/pmcodes?radio_pm_search=param_group&pm_group=All+--+include+all+parameter+groups&pm_search=&casrn_search=&srsname_search=&format=html_table&show=parameter_group_nm&show=parameter_nm&show=casrn&show=srsname&show=parameter_units
# 
##############################################
# 0205551460 LICK RUN ABOVE PATTON AVENUE AT ROANOKE, VA
##############################################
make_usgs_gage_list=function(siteNo = "0205551460",
parameterCd = c("00060","00065"),
start.date = "2017-05-01",  # Not frozen to not frozen
end.date = "2017-11-01"){   # to still not frozen

USGSlist=list()   # Organize the data in a nice list as in previous labs
USGSlist[["flowdata"]]<- readNWISuv(siteNumbers = siteNo,parameterCd = parameterCd,startDate = start.date,endDate = end.date)
head(USGSlist$flowdata)  # Note that we have 00060 and 00065...
#  agency_cd	site_no        	dateTime X_00060_00000 X_00060_00000_cd
#1  	USGS 0205551460 2017-05-01 04:00:00      	6.38            	A
#2  	USGS 0205551460 2017-05-01 04:05:00      	6.38            	A
#  X_00065_00000 X_00065_00000_cd tz_cd
#1      	2.74            	A   UTC
#2      	2.74            	A   UTC
#
# And of course we want to work in SI units so:
USGSlist$flowdata$depth_m=USGSlist$flowdata$X_00065_00000*0.3048
# m/ft depth
USGSlist$flowdata$cms=USGSlist$flowdata$X_00060_00000*.02832
# m3/ft3 flow
#
# Let's add in the USGS gage site information to the list and inspect
USGSlist[["site"]]=readNWISsite(siteNo)
head(USGSlist$site)
class(USGSlist$site$dec_lat_va)
#
# Set the Manning Coefficient in the USGS Gage's Site Table
#
USGSlist$site$man_n=.035/1.49
#
# Create a SpatialPointsDataFrame out of the site dataframe in the USGS list
coordinates(USGSlist$site)=~dec_long_va+dec_lat_va
#
return(USGSlist)
}

USGS02056000=make_usgs_gage_list(siteNo = "02056000")
USGS0205551460=make_usgs_gage_list(siteNo ="0205551460" )
USGS02055100=make_usgs_gage_list(siteNo ="02055100" )
USGS02055000=make_usgs_gage_list(siteNo ="02055000" )
USGS02054530=make_usgs_gage_list(siteNo ="02054530" )

ab_ll=rbind(USGS02056000$site,
              USGS0205551460$site,
              USGS02055100$site,
              USGS02055000$site,
              USGS02054530$site)
class(ab_ll)
ab_ll@proj4string
proj4_utm = paste0("+proj=utm +zone=",
                     trunc((180+coordinates(USGS02055000$site)[1])/6+1), 
                     " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
proj4string(ab_ll)=proj4_ll
ab_utm=spTransform(ab_ll,crs_utm)
ab_utm@coords
mydem=get_aws_terrain(locations=ab_utm@coords, 
                        z = 12, prj = proj4_utm,expand=1)
#
# Lets plot the DEM and the gage locations so we can guess 
# what gages connect with what gages
#
dev.off()
plot(mydem)
plot(ab_utm,add=T)
text(ab_utm, labels=ab_utm@data$site_no, cex=0.6, font=2,pos=1)
# From Lab02, I know I can get an overview of streams with the 
# USGS H
url="https://prd-tnm.s3.amazonaws.com/StagedProducts/Hydrography/NHD/HU8/HighResolution/Shape/NHD_H_03010101_HU8_Shape.zip"
curl_download(url,"NHD_H_03010101_HU8_Shape.zip")
unzip("NHD_H_03010101_HU8_Shape.zip",exdir="03010101")
streams=readOGR("03010101/Shape/NHDFlowline.dbf")
streams_utm=spTransform(streams,crs_utm)
plot(streams_utm,col="blue",add=T)
#zoom(mydem)

# A quick readthrough of the Example 1: Hiking around Maunga Whau
# in the package vignette. 
 vignette("Overview", package = "gdistance")
# Set the starting and ending locations
# determine the river reach length and slope using the gdistance package.
#
A=SpatialPoints(USGS0205551460$site)# Up gradient site Lick Run
B=SpatialPoints(USGS02056000$site) # Down gradient site ROA River at Niagara
proj4string(A)=proj4_ll
proj4string(B)=proj4_ll
A_utm=spTransform(A,crs_utm)
B_utm=spTransform(B,crs_utm)
# Cut the DEM down to a more manageable size
cropmydem=crop(mydem,extend(extent(ab_utm),600))
cropmydem=trim(cropmydem)
cropmydem=cropmydem*1000.0
plot(cropmydem)
plot(ab_utm,add=T)
plot(B_utm, add =T, col = "red")
# Set up the weighting functions
altDiff <- function(x){x[2] - x[1]}
hd <- transition(cropmydem, altDiff, 8, symm=FALSE)
slope <- geoCorrection(hd)
adj <- adjacent(cropmydem, cells=1:ncell(cropmydem), pairs=TRUE, directions=8)
speed <- slope
speed[adj] <- 6 * exp(-3.5 * abs(slope[adj] + 0.05))
Conductance <- geoCorrection(speed)
# Find and plot the flow path


route <- function(A_utm, Gauge, n, B_utm) {
AtoB <- shortestPath(Conductance, A_utm, B_utm, output="SpatialLines")
#dev.off()
#plot(AtoB,add=T)
#plot(streams_utm,col="blue",add=T)
#plot(AtoB,add=T)
SpatialLinesLengths(AtoB)
L=SpatialLinesLengths(AtoB) # km to m
L # reach length in m
#zoom()
#
# Getting slope, we will extract the slope for points A and B from the DEM and # divide the difference by the length in m, this gives us a much better 
# estimate of slope than taking the point slopes at the gage site
#
slope=(extract(mydem,A_utm)-
                               extract(mydem,B_utm))/L
slope

attach(Gauge)

B=(Gauge$site$man_n*
     Gauge$flowdata$cms)/(Gauge$flowdata$depth_m^(5/3)*
                                                               sqrt(slope))
#head(Gauge$flowdata)
#  agency_cd	site_no        	dateTime X_00060_00000 X_00060_00000_cd
#1  	USGS 05267000 2017-05-01 04:00:00      	6.38            	A
#2  	USGS 05267000 2017-05-01 04:05:00      	6.38            	A
#  X_00065_00000 X_00065_00000_cd tz_cd   	cms  depth_m    	B
#1      	2.74            	A   UTC 0.1806816 0.835152 0.103032
#2      	2.74            	A   UTC 0.1806816 0.835152 0.103032
#
# Lets look at how B changes with flow.    
#plot(Gauge$flowdata$dateTime,Gauge$flowdata$B, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")
# Does this seem reasonable (...like order of magnitude reasonable)? You can 
# perform a quick and dirty check using google earth and measuring the channel 
# width in a few places.
#
#plot(Gauge$flowdata$cms,Gauge$flowdata$depth_m, main="LICK RUN TO ROANOKE RIVER AT NIAGARA, VA")


# ck
#Gauge$flowdata$ck = ???
  # ANS
  ck =
  5/3*sqrt(slope)/Gauge$site$man_n*
  (Gauge$flowdata$depth_m^(2/3))
#
#Gauge$flowdata$dt = ???
  dt =
  L/ck

#plot(Gauge$flowdata$dateTime,Gauge$flowdata$dt)
outTime=Gauge$flowdata$dateTime+
  dt

# Find beginning of  Waves
newwave=
  Gauge$flowdata$cms *1.1 <
  data.table::shift(Gauge$flowdata$cms, n = n)

# Add plot of the point found
len=length(newwave)
newwave[is.na(newwave)]=F              # Investigate this
# Removes repeated finds by going through loop backwords
for (i in seq(len,2)){
  print(i)
  if(newwave[i]==T &
     newwave[i-1]==T){
    newwave[i]=F
  }
}
Gauge$site$L = L
Gauge$site$slope = slope
Gauge$flowdata$B = B
Gauge$flowdata$ck = ck
Gauge$flowdata$dt = dt
Gauge$flowdata$outTime = outTime
Gauge$flowdata$newwave = newwave
detach(Gauge)
rm(L, slope, B, ck, dt, outTime, newwave)
return(Gauge)
}


#route flow from upstream gauge A to downstream gauge B and calculate travel time dt
A1=SpatialPoints(USGS02055100$site)# Up gradient site Lick Run
A2=SpatialPoints(USGS02054530$site)# Up gradient site Lick Run
A3=SpatialPoints(USGS02055000$site)# Up gradient site Lick Run
proj4string(A1)=proj4_ll
proj4string(A2)=proj4_ll
proj4string(A3)=proj4_ll
A1_utm=spTransform(A1,crs_utm)
A2_utm=spTransform(A2,crs_utm)
A3_utm=spTransform(A3,crs_utm)

plot(cropmydem)
points(A1_utm, col = "red")
points(A2_utm, col = "black") 
points(A3_utm, col = "blue")
points(B_utm, col = "purple")
lines(streams_utm, add = T)

USGS02055100 = route(A_utm = A1_utm, Gauge =  USGS02055100, n=5, B_utm = B_utm)

plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms,type="l")
points(USGS02055100$flowdata$dateTime[USGS02055100$flowdata$newwave],
       USGS02055100$flowdata$cms[USGS02055100$flowdata$newwave],col=2)
which(USGS02055100$flowdata$newwave == TRUE)

##########
USGS02054530 = route(A_utm = A2_utm, Gauge =  USGS02054530, n=5, B_utm = B_utm)

plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$cms,type="l")
points(USGS02054530$flowdata$dateTime[USGS02054530$flowdata$newwave],
       USGS02054530$flowdata$cms[USGS02054530$flowdata$newwave],col=2)
which(USGS02054530$flowdata$newwave == TRUE)

########################
USGS02055000 = route(A_utm = A3_utm, Gauge =  USGS02055000, n=5, B_utm = B_utm)

plot(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$cms,type="l")
points(USGS02054530$flowdata$dateTime[USGS02054530$flowdata$newwave],
       USGS02054530$flowdata$cms[USGS02054530$flowdata$newwave],col=2)
which(USGS02054530$flowdata$newwave == TRUE)
# Find the time locations where waves begin
dev.off()
plot(USGS02055100$flowdata$dateTime,USGS02055100$flowdata$cms,
       type="l",xlim=c(USGS02055100$flowdata$dateTime[414],
                       USGS02055100$flowdata$dateTime[414+200]), ylim=c(0,300))
lines(USGS02055100$flowdata$outTime,USGS02055100$flowdata$cms,col=2)
lines(USGS02056000$flowdata$dateTime,USGS02056000$flowdata$cms,col=5, ylim=c(0,200))
lines(USGS02054530$flowdata$dateTime,USGS02054530$flowdata$cms,
     type="l",xlim=c(USGS02054530$flowdata$dateTime[414],
                     USGS02054530$flowdata$dateTime[414+200]))
lines(USGS02054530$flowdata$outTime,USGS02054530$flowdata$cms,col=7)
lines(USGS02055000$flowdata$dateTime,USGS02055000$flowdata$cms,
      type="l",xlim=c(USGS02055000$flowdata$dateTime[414],
                      USGS02055000$flowdata$dateTime[414+200]))
lines(USGS02055000$flowdata$outTime,USGS02055000$flowdata$cms,col=12)


pdf("HW2.pdf")

ggplot() +
  geom_line(data=USGS02056000$flowdata, aes(x= dateTime, y = cms, color = "Observed")) + 
  geom_line(data=USGS02055100$flowdata, aes(x=dateTime,y=cms, color = "Observed")) +
  geom_line(data=USGS02054530$flowdata, aes(x=dateTime,y=cms, color = "Observed")) +
  geom_line(data=USGS02055000$flowdata, aes(x=dateTime,y=cms, color = "Observed")) +
  geom_line(data=USGS02055100$flowdata, aes(x=outTime,y=cms, color = "USGS02055100_cms")) +
  geom_line(data=USGS02054530$flowdata, aes(x=outTime,y=cms, color = "USGS02054530_cms")) +
  geom_line(data=USGS02055000$flowdata, aes(x=outTime,y=cms, color = "USGS02055000_cms")) +
#  scale_color_manual("Observed" = "black")
#  scale_x_date() +
  coord_cartesian(xlim = c(USGS02055100$flowdata$dateTime[400],
                           USGS02055100$flowdata$dateTime[410+200])) +
  xlab(NULL)+
  ylab("Discharge, cms") +
  theme_classic()

dev.off() 
 
x1 <- USGS02055100$site$L
x2 <- USGS02054530$site$L
x3 <- USGS02055000$site$L
t1 <- USGS02055100$flowdata$dt[414]/3600
t2 <- USGS02054530$flowdata$dt[449]/3600
t3 <- USGS02055000$flowdata$dt[449]/3600

lag <- data.frame(gauge = c("USGS02055100","USGS02054530","USGS02055000"), dist=c(x1,x2,x3), time=c(t1,t2,t3))
lag
