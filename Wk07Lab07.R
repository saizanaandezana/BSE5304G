objects()
rm(list=objects())
#
# Build a working directory for this weeks lab and change working dir
# Note you might have to specify the path explicitly  as some 
# computers in the lab were not working correctly, to do this go to
# Session>Set Working Directory
dir.create("~/Week07Lab07/HW")
setwd("~/Week07Lab07/HW")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(EcoHydRology,curl,httr,rnoaa,raster,shapefiles,rgdal,elevatr,soilDB, tidyverse, patchwork, rgeos, FedData)

# Get our gold standard flow data from USGS 12115700 BOULDER CREEK NEAR CEDAR FALLS, WA `
myflowgage_id="12115700"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",
                         end_date = "2021-01-01")

# We want Q in mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3

## We can determine our UTM zone and build our proj4 projection string
trunc((180+myflowgage$declon)/6+1)
proj4_utm = paste0("+proj=utm +zone=", trunc((180+myflowgage$declon)/6+1), " +datum=WGS84 +units=m +no_defs")
print(proj4_utm)
# Lat/Lon (_ll) is much easier!
proj4_ll = "+proj=longlat"

# Now we will build our proj4strings which define our “Coordinate 
# Reference Systems” or CRS in future geographic manipulations. 
crs_ll=CRS(proj4_ll)
crs_utm=CRS(proj4_utm)
print(crs_ll)
print(crs_utm)

myflowgage$area   # area in km2
# If the watershed was square, which it is not, the size would be the 
# square root of the area. Also, the gage/pour point is not in the center
# so we will search around the gage.
# Build sp point for USGS gage location, in both _ll and _utm
latlon <- cbind(myflowgage$declon,myflowgage$declat)
myflowgage$gagepoint_ll <- SpatialPoints(latlon)
proj4string(myflowgage$gagepoint_ll)=proj4_ll
myflowgage$gagepoint_utm=spTransform(myflowgage$gagepoint_ll,crs_utm)
# Open up maps.google.com to guesstimate area/lengths
url=paste0("https://www.google.com/maps/@",
           myflowgage$declat,",",myflowgage$declon,",18z")
browseURL(url)

# For our search we are going to multiply the area by 6 and
# to get the distance
searchlength=sqrt(myflowgage$area*3)*1000
pourpoint=SpatialPoints(myflowgage$gagepoint_utm@coords,proj4string = crs_utm)
bboxpts=myflowgage$gagepoint_utm@coords
bboxpts=rbind(bboxpts,bboxpts+searchlength)
bboxpts=rbind(bboxpts,bboxpts-searchlength)
bboxpts
bboxpts=rbind(bboxpts,c(min(bboxpts[,1]),max(bboxpts[,2])))
bboxpts=rbind(bboxpts,c(max(bboxpts[,1]),min(bboxpts[,2])))
bboxpts
bboxpts=SpatialPoints(bboxpts,proj4string = crs_utm)
# From Lab04, get your DEM
mydem=get_aws_terrain(locations=bboxpts@coords, 
                      z = 12, prj = proj4_utm,src ="aws",expand=0)
res(mydem)
plot(mydem)
plot(bboxpts,add=T)
plot(pourpoint,add=T,col="red")

# Write our raster to a geotiff file that can be used with
# OS level hydrological models 
writeRaster(mydem,filename = "mydem.tif",overwrite=T)

old_path <- Sys.getenv("PATH")
old_path
Sys.setenv(PATH = paste(old_path,
                        paste0(Sys.getenv("HOME"),"/TauDEM/bin"), 
                        sep = ":"))
system("mpirun aread8")
z=raster("mydem.tif")
plot(z)

# Pitremove
system("mpiexec -n 8 pitremove -z mydem.tif -fel mydemfel.tif")
fel=raster("mydemfel.tif")
plot(fel)

plot(z-fel)

# D8 flow directions
system("mpiexec -n 8 d8flowdir -p mydemp.tif -sd8 mydemsd8.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
p=raster("mydemp.tif")
plot(p)
sd8=raster("mydemsd8.tif")
plot(sd8)

# Contributing area
system("mpiexec -n 8 aread8 -p mydemp.tif -ad8 mydemad8.tif")
ad8=raster("mydemad8.tif")
plot(log(ad8))
# zoom(log(ad8))


# Grid Network 
system("mpiexec -n 8 gridnet -p mydemp.tif -gord mydemgord.tif -plen mydemplen.tif -tlen mydemtlen.tif")
gord=raster("mydemgord.tif")
plot(gord)
#zoom(gord)

# DInf flow directions
system("mpiexec -n 8 dinfflowdir -ang mydemang.tif -slp mydemslp.tif -fel mydemfel.tif",show.output.on.console=F,invisible=F)
ang=raster("mydemang.tif")
plot(ang)
slp=raster("mydemslp.tif")
plot(slp)


# Dinf contributing area
system("mpiexec -n 8 areadinf -ang mydemang.tif -sca mydemsca.tif")
sca=raster("mydemsca.tif")
plot(log(sca))
# zoom(log(sca))

# Threshold
system("mpiexec -n 8 threshold -ssa mydemad8.tif -src mydemsrc.tif -thresh 2000")
src=raster("mydemsrc.tif")
plot(src)
plot(pourpoint, add=T)
# zoom(src)

# a quick R function to write a shapefile
makeshape.r=function(sname="shape",n=1)
{
  xy=locator(n=n)
  points(xy)
  
  #Point
  dd <- data.frame(Id=1:n,X=xy$x,Y=xy$y)
  ddTable <- data.frame(Id=c(1),Name=paste("Outlet",1:n,sep=""))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 1)
  write.shapefile(ddShapefile, sname, arcgis=T)
}
makeshape.r("ApproxOutlets")
# Move Outlets
system("mpiexec -n 8 moveoutletstostrm -p mydemp.tif -src mydemsrc.tif -o ApproxOutlets.shp -om Outlet.shp")
outpt=read.shp("Outlet.shp")
approxpt=read.shp("ApproxOutlets.shp")

plot(src)
points(outpt$shp[2],outpt$shp[3],pch=19,col=2)
points(approxpt$shp[2],approxpt$shp[3],pch=19,col=4)

# Contributing area upstream of outlet
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif")
ssa=raster("mydemssa.tif")
plot(ssa) 

# Threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc1.tif -thresh 8000")
src1=raster("mydemsrc1.tif")
plot(src1)


# Stream Reach and Watershed
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc1.tif -o Outlet.shp -ord mydemord.tif -tree mydemtree.txt -coord mydemcoord.txt -net mydemnet.shp -w mydemw.tif")
plot(raster("mydemord.tif"), add = T)
# zoom(raster("mydemord.tif"))
plot(raster("mydemw.tif"))

# Plot streams using stream order as width
snet=read.shapefile("mydemnet")
ns=length(snet$shp$shp)
for(i in 1:ns){
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$strmOrder[i])
}

# Peuker Douglas stream definition
system("mpiexec -n 8 peukerdouglas -fel mydemfel.tif -ss mydemss.tif")
ss=raster("mydemss.tif")
plot(ss)
#zoom(ss)

#  Accumulating candidate stream source cells
system("mpiexec -n 8 aread8 -p mydemp.tif -o Outlet.shp -ad8 mydemssa.tif -wg mydemss.tif")
ssa=raster("mydemssa.tif")
plot(ssa)
#zoom(ssa)
#  Drop Analysis
system("mpiexec -n 8 dropanalysis -p mydemp.tif -fel mydemfel.tif -ad8 mydemad8.tif -ssa mydemssa.tif -drp mydemdrp.txt -o Outlet.shp -par 5 500 10 0")

# Stream raster by threshold
system("mpiexec -n 8 threshold -ssa mydemssa.tif -src mydemsrc2.tif -thresh 300")
raster("mydemsrc2.tif")

# Stream network
system("mpiexec -n 8 streamnet -fel mydemfel.tif -p mydemp.tif -ad8 mydemad8.tif -src mydemsrc2.tif -ord mydemord2.tif -tree mydemtree2.dat -coord mydemcoord2.dat -net mydemnet2.shp -w mydemw2.tif -o Outlet.shp",show.output.on.console=F,invisible=F)

plot(raster("mydemw2.tif"))
snet=read.shapefile("mydemnet2")
ns=length(snet$shp$shp)
for(i in 1:ns){
  lines(snet$shp$shp[[i]]$points,lwd=snet$dbf$dbf$Order[i])
}

# Wetness Index
system("mpiexec -n 8 slopearearatio -slp mydemslp.tif -sca mydemsca.tif -sar mydemsar.tif", show.output.on.console=F, invisible=F)
sar=raster("mydemsar.tif")
wi=sar
wi[,]=-log(sar[,])
plot(wi)

# Distance Down
system("mpiexec -n 8 dinfdistdown -ang mydemang.tif -fel mydemfel.tif -src mydemsrc2.tif -m ave v -dd mydemdd.tif",show.output.on.console=F,invisible=F)
plot(raster("mydemdd.tif"))


mydemw=raster("mydemw.tif")
mybasinmask=trim(mydemw,padding=2)
mybasindem=crop(mydem,mybasinmask)
mybasindem=mask(mybasindem,mybasinmask)
plot(mybasindem)

# Wetness Index
mybasinslp=crop(slp,mybasinmask)
mybasinslp=mask(mybasinslp,mybasinmask)
plot(mybasinslp)

mybasinsca=crop(sca,mybasinmask)
mybasinsca=mask(mybasinsca,mybasinmask)
plot(mybasinsca)

TI = log( (mybasinsca+1)/(mybasinslp+0.00001) )
plot(TI)
#zoom(TI)

pacman::p_load(classInt)
nTIclass=5 #number of TI classes, currently equal area, can adjust method various ways e.g., classIntervals(v, n = nTIclass, style = "jenks")
v=values(TI)
v=v[!is.na(v)]
brks.qt = classIntervals(v, n = nTIclass, style = "quantile")$brks #length nTIclass+1 of just the numeric breakpoints
TIC = cut(TI, breaks=brks.qt, include.lowest = T, right=T)
#
# A series of plots to show all of the components
#
par(mfrow = c(2, 3))
plot(TI)
plot(TIC)
plot(mybasinsca)
plot(mybasinslp)
z = crop(z, mybasinmask)
z = mask(z, mybasinmask)
fel = crop(fel, mybasinmask)
fel = mask(fel, mybasinmask)
dif = z-fel
plot(dif)
writeRaster(dif,filename = "dif.tif",overwrite=T)
dev.off()
plot(TIC,col=rainbow(5))


ddf <- rasterToPoints(dif); ddf <- data.frame(ddf)
colnames(ddf) <- c("X","Y","DEM")
b.dem <- seq(min(ddf$DEM),max(ddf$DEM),length.out=2)
ramp <- colorRamp(c("gray", "brown"))
p1 <- ggplot()+
  geom_tile(data=ddf,mapping=aes(X,Y,fill=DEM))+
  scale_fill_gradientn(name="Difference",colours = rgb(ramp(seq(0, 1, length = 3)), max = 255),breaks=b.dem)+
#  scale_x_continuous(name=expression(paste("Longitude (",degree,")")),limits=c(-4,2),expand=c(0,0))+
#  scale_y_continuous(name=expression(paste("Latitude (",degree,")")),limits=c(4,12),expand=c(0,0))+
  coord_equal()
p1

Filled <- rasterToPoints(z); Filled <- data.frame(Filled)
colnames(Filled) <- c("X","Y","DEM")
b.dem <- seq(min(Filled$DEM),max(Filled$DEM),length.out=11)

p2 <- ggplot()+
  geom_tile(data=Filled,mapping=aes(X,Y,fill=DEM))+
  scale_fill_gradientn(name="Elevation",colours = rainbow(seq(0,1 - 1/12,length.out = 11)),breaks=b.dem)+
  #  scale_x_continuous(name=expression(paste("Longitude (",degree,")")),limits=c(-4,2),expand=c(0,0))+
  #  scale_y_continuous(name=expression(paste("Latitude (",degree,")")),limits=c(4,12),expand=c(0,0))+
  coord_equal()
p2

TIdf <- rasterToPoints(TI); TIdf <- data.frame(TIdf)
colnames(TIdf) <- c("X","Y","TI")
b.TI <- seq(min(TIdf$TI),max(TIdf$TI),length.out=5)
ramp <- colorRamp(c("gray", "blue"))
p3 <- ggplot()+
  geom_tile(data=TIdf,mapping=aes(X,Y,fill=TI))+
  scale_fill_gradientn(name="TI",colours = rgb(ramp(seq(0, 1, length = 3)), max = 255),breaks=b.TI)+
  #  scale_x_continuous(name=expression(paste("Longitude (",degree,")")),limits=c(-4,2),expand=c(0,0))+
  #  scale_y_continuous(name=expression(paste("Latitude (",degree,")")),limits=c(4,12),expand=c(0,0))+
  coord_equal()

TICdf <- rasterToPoints(TIC); TICdf <- data.frame(TICdf)
colnames(TICdf) <- c("X","Y","TIC")
b.TIC <- seq(min(TICdf$TIC),max(TICdf$TIC),length.out=5)
ramp <- colorRamp(c("gray", "blue"))
p4 <- ggplot()+
  geom_tile(data=TICdf,mapping=aes(X,Y,fill=TIC))+
  scale_fill_gradientn(name="TIC",colours = rgb(ramp(seq(0, 1, length = 5)), max = 255),breaks=b.TIC)+
  #  scale_x_continuous(name=expression(paste("Longitude (",degree,")")),limits=c(-4,2),expand=c(0,0))+
  #  scale_y_continuous(name=expression(paste("Latitude (",degree,")")),limits=c(4,12),expand=c(0,0))+
  coord_equal()

pdf("BOULDER CREEK_1.pdf")
p1+p2+p3+p4+plot_layout(ncol = 2, nrow = 2)
dev.off()

# create polygon using raster subwatershed boundaries
mydemw_poly=rasterToPolygons(mydemw,dissolve = T,na.rm = T) # takes to long to compute
plot(mydemw_poly,add=T,border="red")
mydemw_poly
# go and grab soils data using the polygon from websoilsurvey
x <- mapunit_geom_by_ll_bbox(bbox(spTransform(mydemw_poly,crs_ll)))
x=spTransform(x,crs_utm) # convert LL to UTM
out <- gIntersection(x, mydemw_poly, byid=TRUE)
plot(out,col=rainbow(5, 0.5))

outraster <- rasterize(out, mydemw, fun='first')
plot(outraster, ext = mydemw_poly)


library(FedData)
polygon1 <- polygon_from_extent(raster::extent(592723.1, 604662.2, 5240720, 5252660), 
                                  proj4string='+proj=utm +datum=WGS84 +zone=10')

LU <- get_nlcd(polygon1, label = "creek", year = 2011, dataset = "landcover") # this funcition is under maintenance and returns error


soildf <- rasterToPoints(outraster); soildf <- data.frame(soildf)
colnames(soildf) <- c("X","Y","soil")
b.soil <- seq(min(soildf$soil),max(soildf$soil),length.out=10)
ramp <- colorRamp(c("gray", "brown"))
p5 <- ggplot()+
  geom_tile(data=soildf,mapping=aes(X,Y,fill=soil))+
  scale_fill_gradient2(name="soil", low = "gray", mid = "yellow", high = "brown")+
  #  scale_x_continuous(name=expression(paste("Longitude (",degree,")")),limits=c(-4,2),expand=c(0,0))+
  #  scale_y_continuous(name=expression(paste("Latitude (",degree,")")),limits=c(4,12),expand=c(0,0))+
  coord_equal() 
#  theme_classic()
p5

pdf("BOULDER CREEK_2.pdf")
p2 + p4 + p5 + plot_layout(ncol = 1, nrow = 3)
dev.off()
