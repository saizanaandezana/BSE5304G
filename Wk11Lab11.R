objects()
rm(list=objects())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,EcoHydRology,curl,data.table,multisensi, sensitivity)


##################### PET_fromTemp() ##################

T <- seq(from = 5, to = 365, by = 5)
#Solution for PET_fromTemp
# Trick is you have to notice that "lat_radians" has replaced "lat" and
# there is no "units" variable... and... notice that the function has to
# be fixed to allow Jday to be a vector of different size than Tmax and Tmin


PET_fromTemp <- function (Jday, Tmax_C, Tmin_C, lat_rad, AvgT = (Tmax_C + Tmin_C)/2, 
                          albedo = 0.18, TerrestEmiss = 0.97, aspect = 0, slope = 0, 
                          forest = 0, PTconstant=1.26, AEparams=list(vp=NULL, opt="linear"))
{
  cloudiness <- EstCloudiness(Tmax_C, Tmin_C)
  DailyRad <- NetRad(lat_rad, Jday, Tmax_C, Tmin_C, albedo, forest, slope, 
                     aspect, AvgT, cloudiness, TerrestEmiss, AvgT, AEparams=AEparams)
  potentialET <- PTpet(DailyRad, AvgT, PTconstant)
  potentialET[which(potentialET < 0)] <- 0
  potentialET[which(Tmax_C == -999 | Tmin_C == -999)] <- (-999)
  return(potentialET)
}



PET_fromTemp2 <- function(X, t = T) {
  out <- matrix(nrow = nrow(X), ncol = length(t), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- PET_fromTemp(lat_rad=X$lat_rad[i], 
                             Tmax_C=X$Tmax_C[i], 
                             Tmin_C=(X$Tmax_C[i]-X$Trange[i]),
                             albedo=X$albedo[i],
                             slope=X$slope[i],
                             aspect=X$aspect[i],
                             Jday=t)
  }
  out <- as.data.frame(out)
  names(out) <- paste("t", t, sep = "")
  return(out)
}


PET.seq.fast <- multisensi(design = fast99, model = PET_fromTemp2,
                             center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                             design.args=list(factors=c("lat_rad",
                                                         "Tmax_C",
                                                         "Trange",
                                                         "albedo",
                                                         "slope",
                                                         "aspect"),
                                               n = 1000, 
                                               q = "qunif",
                                               q.arg = list(
                                                 list(min = 0, max = pi/3), 
                                                 list(min = 1, max = 40), 
                                                 list(min = 1, max = 10), 
                                                 list(min = 0, max = 1), 
                                                 list(min = 0, max = 0.2), 
                                                 list(min = 0, max = pi*2))),
                             analysis.args=list(keep.outputs=FALSE))


pdf("PET.pdf")
plot(PET.seq.fast, normalized = T, color = terrain.colors, gsi.plot = F)
title(xlab = "Time in half-decades")
plot(PET.seq.fast, normalized = F, color = terrain.colors, gsi.plot = F)
title(xlab = "Time in half-decades")
dev.off()

##################### NetRad() ##################

NetRad(lat, Jday, Tx, Tn, albedo = 0.18, forest = 0, slope = 0, 
       aspect = 0, airtemp = (Tn+Tx)/2, cloudiness = "Estimate", 
       surfemissivity = 0.97, surftemp = (Tn+Tx)/2, units = "kJm2d", 
       AEparams=list(vp=NULL, opt="linear"))

NetRad2 <- function(X, t = T) {
  out <- matrix(nrow = nrow(X), ncol = length(t), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- NetRad(lat=X$lat[i], 
                      Tx=X$Tx[i], 
                      Tn=(X$Tx[i]-X$Trange[i]),
                      albedo=X$albedo[i],
                      slope=X$slope[i],
                      aspect=X$aspect[i],
                      Jday=t, units = "Wm2")
    
  }
  out <- as.data.frame(out)
  names(out) <- paste("t", t, sep = "")
  return(out)
}


NetRad.seq.fast <- multisensi(design = fast99, model = NetRad2,
                             center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                             design.args=list( factors=c("lat",
                                                         "Tx",
                                                         "Trange",
                                                         "albedo",
                                                         "slope",
                                                         "aspect"),
                                               n=1000, q = "qunif",
                                               q.arg = list(
                                                 
                                                 list(min = 0, max = pi/3), 
                                                 list(min = 1, max = 40), 
                                                 list(min = 1, max = 10), 
                                                 list(min = 0, max = 1), 
                                                 list(min = 0, max = 0.2), 
                                                 list(min = 0, max = pi*2))),
                             analysis.args=list(keep.outputs=FALSE))
pdf("NetRd.pdf")
plot(NetRad.seq.fast, normalized = T, color = terrain.colors, gsi.plot = F)
title(xlab = "Time in half-decades")

plot(NetRad.seq.fast, normalized = F, color = terrain.colors, gsi.plot = F)
title(xlab = "Time in half-decades")
dev.off()



############## SoilStorage(S_avg, field_capacity, soil_water_content, porosity) #######
CN <- seq(95, 40, -5)
S <- 25400/CN - 254

porosity <- c(0.2,0.4)
swc <- c(min(porosity)*0.25, max(porosity)*0.95)

SoilStorage2 <- function(X, s=S) { # S for CN(40, 95)
  out <- matrix(nrow = nrow(X), ncol = length(s), NA)
  for (i in 1:nrow(X)) {
    out[i, ] <- SoilStorage(S_avg=s, 
                       field_capacity=X$FC[i], 
                       soil_water_content=X$swc[i],
                       porosity=X$porosity[i]
                       )
  }
  out <- as.data.frame(out)
  names(out) <- paste("S", round(S,0), sep = "")
  return(out)
}


SoilStorage.seq.fast <- multisensi(design = fast99, model = SoilStorage2,
                              center = FALSE, reduction = NULL, analysis = analysis.sensitivity,
                              design.args=list( factors=c(
                                                          "FC",
                                                          "swc",
                                                          "porosity"),
                                                n=1000, q = "qunif",
                                                q.arg = list(
                                                  
                                                    # CN 40 to 95
                                                  list(min = 0.09, max = 0.38),       # sand to clay
                                                  list(min = min(swc), max = max(swc)), 
                                                  list(min = min(porosity), max = max(porosity)))),
                              analysis.args=list(keep.outputs=FALSE))

pdf("SoilStorage.pdf")
plot(SoilStorage.seq.fast, normalized = F, color = terrain.colors, gsi.plot = F)
title(xlab = "Storage for given watershed property")
plot(SoilStorage.seq.fast, normalized = T, color = terrain.colors, gsi.plot = F)
title(xlab = "Storage for given watershed property")
dev.off()
