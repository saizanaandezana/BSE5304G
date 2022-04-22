rm(list=objects())
setwd("~")
dir.create("Wk12Lab12")
setwd("~/Wk12Lab12/")
#install.packages("EcoHydRology", repos="http://R-Forge.R-project.org")
library(EcoHydRology)
pacman::p_load(rnoaa,ggplot2,fastmatch,forcats,dplyr)
#--------------source CN model function from previous lab----------
source("https://raw.githubusercontent.com/vtdrfuka/BSE5304_2022/main/functions/CNmodel")
#LITTLE OTTER CREEK AT FERRISBURG, VT.
#---------Getting streamflow data from USGS-------------------
myflowgage_id="12115700"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2010-01-01",end_date = "2022-04-20")
# Note that flow returned is in m3/day, but we want mm/day for the basin
myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
#--------Getting weather data from FillMissWX (noaa) function
WXData=FillMissWX(declat=myflowgage$declat, declon=myflowgage$declon,60,    ### How does this function consider data homogeneity (how many station get considered?)
                    date_min="2010-01-01",
                    date_max="2022-04-20")
#---building modeldata dataframe by merging streamflow and weather data-------

WXData[length(WXData$date),] = WXData[length(WXData$date)-1,]


modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")



#---Running CNmodel function to check the model-------------------------------
CN_DF = CNmodel(CNmodeldf = modeldata,CNavg = 70,
                IaFrac = 0.05,fnc_slope =  0.05,
                fnc_aspect =0.0,func_DAWC = 0.35,
                func_z=1000,fnc_fcres=0.20,
                declat=myflowgage$declat,
                declon=myflowgage$declon)
NSeff(CN_DF$Qmm, CN_DF$Qpred)
#-----------getting ready for optimization using DEoptim----------------------
f <- function (x) {
  CNopt=x[1]
  IaOpt=x[2]
  fnc_slopeOpt=x[3]
  fnc_aspectOpt=x[4]
  func_DAWCOpt=x[5]
  func_zOpt=x[6]
  fnc_fcresOpt=x[7]
  
  CNmodelnew = CNmodel(CNmodeldf =modeldata,
                     CNavg = CNopt,IaFrac = IaOpt,
                     fnc_slope=fnc_slopeOpt,fnc_aspect=fnc_aspectOpt,
                     func_DAWC=func_DAWCOpt,func_z=func_zOpt,
                     fnc_fcres=fnc_fcresOpt,declat=myflowgage$declat,
                     declon=myflowgage$declon)
  mNSE = 1 - NSeff(CNmodelnew$Qmm,CNmodelnew$Qpred)
  
  return(mNSE)
}
#            Educated guesses for ranges… or not.        
#            CN,  Ia, Slp, Asp, AWC,Depth,ResCoef
lower <- c(30,0.01, 0.0, 0.0, 0.2,  500,   0.1)
upper <- c(99,0.20, 0.2, 4, 0.4, 2000,   0.5)

# run DEoptim and set a seed first for replicability
set.seed(123)
#
# The cluster will show you 128 cores, but the management software
# will limit you to running on however many cores you requested.
#
cl <- parallel::makeCluster(16)
outDEoptim = DEoptim(f, lower, upper,
                     DEoptim.control(cluster = cl, strategy = 6,
                                     NP = 16, itermax = 100, parallelType = 1,
                                     packages = c("EcoHydRology"),
                                     parVar = c("CNmodel","SnowMelt","PET_fromTemp",
                                              "SoilStorage","NSeff","modeldata","myflowgage")))

outDEoptim$optim$bestmem
# NSE = 1- bestval
1-outDEoptim$optim$bestval
x=outDEoptim$optim$bestmem
f(x)
names(outDEoptim$member$lower) = c("CNavg","IaFrac","fnc_slope", 
                                   "fnc_aspect","func_DAWC","func_z",
                                   "fnc_fcres")
plot(outDEoptim)
plot(outDEoptim,plot.type="bestvalit")
#plot(outDEoptim,plot.type="bestmemit")
dev.off()

##-------plot NSE vs iteration------------
NSE = data.frame(1-outDEoptim$member$bestvalit)
colnames(NSE) = "NSE"
plot(NSE$NSE,type = "l", ylab = "NSE", xlab = "Iteration number")
#----another way to plot all parameters of differnt iteration----
BestMemit = data.frame(outDEoptim$member$bestmemit)
colnames(BestMemit) = c("CNavg","IaFrac","fnc_slope", 
                        "fnc_aspect","func_DAWC","func_z",
                        "fnc_fcres")
#plotting the progression of improvements similar to above
for(i in 1:length(BestMemit)){
  plot(BestMemit[,i], ylab= "",xlab="",ylim=c(lower[i],upper[i]),cex=0.3)
  mtext(side=1, line=2, "iteration", col="black",cex=1)
  mtext(side=2, line=3, colnames(BestMemit)[i], col="black", font=2, cex=1)
}

#---------merging NSE and BestMemit datafram ---------------------------
Sensitivity_df = cbind(BestMemit, NSE)
#---calculate changes between two succesive rows-----------------------------
deltaparams = (tail(Sensitivity_df,-1)-Sensitivity_df[1:(length(Sensitivity_df[,1])-1),])
# 
#----calculate the relative sensitivity of each parameter------
#       Remember this slide 13 from the sensitivity lecture?

shifted_sensitivity_df = Sensitivity_df[1:(length(Sensitivity_df[,1]))-1,]

for (i in colnames(deltaparams)){
  #i="CNavg"
  nam = paste0("junkpar_",i)
  assign(nam, deltaparams[!deltaparams[,fmatch(i,names(deltaparams))]==0,])
  nam2=paste0("initialval_",i)
  assign(nam2,shifted_sensitivity_df[!deltaparams[,fmatch(i,names(deltaparams))]==0,])
  nam3=paste0("Sr_",i)
  assign(nam3,abs(((get(paste0("junkpar_",i))$NSE)/
                     (get(paste0("junkpar_",i))[fmatch(i,names(deltaparams))]))*
                    (get(paste0("initialval_",i))[fmatch(i,names(deltaparams))]/
                       get(paste0("initialval_",i))$NSE)))
}

# 
# With sr for each params and you can get the mean of them to compare the sr
# Making a dataframe with bbox and violin plots
RelSensi=base::as.data.frame(stack(list(CNavg=Sr_CNavg$CNavg,
                                          IaFrac=Sr_IaFrac$IaFrac,
                                          fnc_slope=Sr_fnc_slope$fnc_slope, 
                                          fnc_aspect=Sr_fnc_aspect$fnc_aspect,
                                          func_DAWC=Sr_func_DAWC$func_DAWC,
                                          func_z=Sr_func_z$func_z,
                                          fnc_fcres=Sr_fnc_fcres$fnc_fcres
)))
names(RelSensi)[2]="Parameter"
# Basic violin plot
p <- ggplot(RelSensi, aes(x=Parameter, y=values, color=Parameter)) + 
  coord_cartesian(ylim = c(0,3))+
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme(legend.position = "None")

p

pdf("HW1.pdf")
p + stat_summary(fun.y=mean, geom="point", shape=8, size=4,col="black") +
  xlab("Parameter")+
  ylab("Relative Sensitivity")+
  theme(axis.text = element_text(size =10))+
  theme(axis.text = element_text(face="bold"))+
  theme(axis.title = element_text(size =10))+
  theme(axis.title = element_text(face="bold"))
dev.off()


###########################################
############ CNCNCNCN ####################
              ##########
              #########3
dat=data.frame(Qm3ps=CNmodelnew$Qmm,
               Jday=strptime(CNmodelnew$date,"%Y-%m-%d")$yday+1,Precip=CNmodelnew$P,
               TMX=CNmodelnew$MaxTemp,TMN=CNmodelnew$MinTemp,
                 jdate=julian(CNmodelnew$date))
# Need the previous days flow
# The NN for this example requires the previous days flow… remember this for 
# the homework questions. 
dat$Qm3ps_lag=dplyr::lag(dat$Qm3ps)
dat[is.na(dat)]=0
View(dat)
head(dat)
datmax <- apply(dat, 2, max)
datmin <- apply(dat, 2, min)
# now split the datat into train, validation and test sets.
scaled <- scale(dat, center = datmin, scale = datmax - datmin)

# sub sample to get the test and train set
index <- seq(1, nrow(dat))
train_index <- sample(index, 0.8*length(index))
train <- scaled[train_index,]
test <- scaled[-train_index,]
#===================================================
#             Now, setup and train the neural network!!!
#===================================================
# first set a random seed to a lucky number, but it does give results!
pacman::p_load(neuralnet,ggplot2,Metrics,reshape2)
set.seed(4)
# setup and train the neural network
net <- neuralnet(Qm3ps ~ Qm3ps_lag+Precip+jdate+Jday+TMX+TMN,
                   data = train,
                   hidden = c(6, 3), # can adjust here
                   threshold = 0.01, # if you network does not run, 
                   #                         try turning this up just a little
                   stepmax = 1e+05, # can adjust if you don't converge
                   rep = 1,
                   #startweights = NULL,
                   #learningrate.limit = NULL,
                   #learningrate.factor = list(minus = 0.5, plus = 1.2),
                   learningrate = 0.0027,
                   lifesign = "full",
                   #lifesign.step = 1000,
                   algorithm = "backprop",
                   err.fct = 'sse',
                   act.fct = "logistic",
                   linear.output = FALSE,
                   exclude = NULL,
                   constant.weights = NULL)
#===================================================
#                     Plot the Structure of the Neural Net
#===================================================
# View the structure of the neural net
plot(net)
#===================================================
#             Make Predictions to See how the neural network did
#===================================================
# look at predictions from the train and test sets
#  to see how the network did

# make the predictions on the train data
nn_prediction = predict(net, train[
  ,c("Qm3ps_lag","Precip","jdate","Jday","TMX","TMN")])

# un-scale the predictions
unscaled_train_prediction = (nn_prediction * (max(dat$Qm3ps) - min(dat$Qm3ps))) + min(dat$Qm3ps)

train_jdate = (train[,c("jdate")] * (max(dat$jdate) - min(dat$jdate))) + min(dat$jdate)
plot(train_jdate,unscaled_train_prediction)
# But, notice how predictions are created
lines(train_jdate,unscaled_train_prediction,type="l")

# make predictions on the test data
nn_prediction = predict(net, test[, c("Qm3ps_lag","Precip","jdate","Jday","TMX","TMN")])

# un-scale
unscaled_test_prediction = (nn_prediction * (max(dat$Qm3ps) - min(dat$Qm3ps))) + min(dat$Qm3ps)
test_jdate = (test[,c("jdate")] * (max(dat$jdate) - min(dat$jdate))) + min(dat$jdate)
#===================================================
#           Plot the predictions for both train and test
#===================================================
# plot the train and test and predictions
train_fit <- data.frame(NNjdate = train_jdate, NNFlow = unscaled_train_prediction, color="red",data_set = 'train')
test_fit = data.frame(NNjdate = test_jdate, NNFlow  = unscaled_test_prediction, color="green",data_set = 'test')
# Build a results dataframe for the NN
NNFlowModel=rbind(train_fit,test_fit)
NNFlowModel=NNFlowModel[order(NNFlowModel$NNjdate),]
#
# Rebuild a date from Julian Date
NNFlowModel$date=as.Date(NNFlowModel$NNjdate,origin = as.Date("1970-01-01"))

# Neatly combine results into a single dataframe.
NNFlowModel_test=subset(NNFlowModel,data_set=="test")
NNFlowModel_train=subset(NNFlowModel,data_set=="train")
plot(CNmodelnew$date, CNmodelnew$Qmm,type="l",col="black")
lines(NNFlowModel_train$date,NNFlowModel_train$NNFlow,col="red",type="l")
lines(NNFlowModel_test$date,NNFlowModel_test$NNFlow,col="green",type="l")
# Compare this NN Model to that of Lab12
lines(CNmodelnew$date, CNmodelnew$Qpred,type="l",col="blue")
legend("topleft",legend = c("Observed Flow","NN Train Flow","NN Test Flow", "CN_Model Flow"),col = c("black","red","green","blue"),lty = 1:2, cex = 0.6)

# Calculate NSEs for the CN_Model and NNFlowModel
comparedf=merge(CNmodelnew,NNFlowModel)
NN_NSE<- NSeff(comparedf$Qmm, comparedf$NNFlow) 
print(NN_NSE)
CNmodelnew_NSE<- NSeff(CNmodelnew$Qmm, CNmodelnew$Qpred)
print(CNmodelnew_NSE)


########################## assuming no record for last three years


index1 <- seq(1, nrow(dat)-1460)
train_index1 <- sample(index1, 0.8*length(index1))
train_1 <- scaled[train_index1,]
test_1 <- scaled[-train_index1,]

set.seed(6)
net_1 <- neuralnet(Qm3ps ~ Qm3ps_lag+Precip+jdate+Jday+TMX+TMN,
                 data = train_1,
                 hidden = c(6, 3), # can adjust here
                 threshold = 0.011, # if you network does not run, 
                 #                         try turning this up just a little
                 stepmax = 1e+05, # can adjust if you don't converge
                 rep = 1,
                 #startweights = NULL,
                 #learningrate.limit = NULL,
                 #learningrate.factor = list(minus = 0.5, plus = 1.2),
                 learningrate = 0.0027,
                 lifesign = "full",
                 #lifesign.step = 1000,
                 algorithm = "backprop",
                 err.fct = 'sse',
                 act.fct = "logistic",
                 linear.output = FALSE,
                 exclude = NULL,
                 constant.weights = NULL)
#===================================================
#                     Plot the Structure of the Neural Net
#===================================================
# View the structure of the neural net
plot(net)
#===================================================
#             Make Predictions to See how the neural network did
#===================================================
# look at predictions from the train and test sets
#  to see how the network did

# make the predictions on the train data
nn_prediction_1 = predict(net_1, train_1[
  ,c("Qm3ps_lag","Precip","jdate","Jday","TMX","TMN")])

# un-scale the predictions
unscaled_train_prediction_1 = (nn_prediction_1 * (max(dat$Qm3ps) - min(dat$Qm3ps))) + min(dat$Qm3ps)

train_jdate_1 = (train_1[,c("jdate")] * (max(dat$jdate) - min(dat$jdate))) + min(dat$jdate)
plot(train_jdate_1,unscaled_train_prediction_1)
# But, notice how predictions are created
lines(train_jdate_1,unscaled_train_prediction_1,type="l")

# make predictions on the test data
nn_prediction_1 = predict(net_1, test_1[, c("Qm3ps_lag","Precip","jdate","Jday","TMX","TMN")])

# un-scale
unscaled_test_prediction_1 = (nn_prediction_1 * (max(dat$Qm3ps) - min(dat$Qm3ps))) + min(dat$Qm3ps)
test_jdate_1 = (test_1[,c("jdate")] * (max(dat$jdate) - min(dat$jdate))) + min(dat$jdate)
#===================================================
#           Plot the predictions for both train and test
#===================================================
# plot the train and test and predictions
train_fit_1 <- data.frame(NNjdate = train_jdate_1, NNFlow = unscaled_train_prediction_1, color="red",data_set = 'train')
test_fit_1 = data.frame(NNjdate = test_jdate_1, NNFlow  = unscaled_test_prediction_1, color="green",data_set = 'test')
# Build a results dataframe for the NN
NNFlowModel_1=rbind(train_fit_1,test_fit_1)
NNFlowModel_1=NNFlowModel_1[order(NNFlowModel_1$NNjdate),]
#
# Rebuild a date from Julian Date
NNFlowModel_1$date=as.Date(NNFlowModel_1$NNjdate,origin = as.Date("1970-01-01"))

# Calculate NSEs for the CN_Model and NNFlowModel
comparedf_1=merge(CNmodelnew,NNFlowModel_1)
NN_NSE_1<- NSeff(comparedf$Qmm, comparedf_1$NNFlow) 
print(NN_NSE)
CNmodelnew_NSE<- NSeff(CNmodelnew$Qmm, CNmodelnew$Qpred)
print(CNmodelnew_NSE)
# Neatly combine results into a single dataframe.
NNFlowModel_test_1=subset(NNFlowModel_1,data_set=="test")
NNFlowModel_train_1=subset(NNFlowModel_1,data_set=="train")

pdf("HW3.pdf")
ggplot() +
  geom_line(data=NNFlowModel_test_1, aes(x = date, y = NNFlow, color = "NNFlow_test")) +
  geom_vline(xintercept = NNFlowModel_train_1$date[length(NNFlowModel_train_1$date)]) +
  geom_line(data=NNFlowModel_train_1, aes(x = date, y = NNFlow, color = "NNFlow_train")) +
  geom_point(data=CNmodelnew, aes(x = date, y = Qmm, color = "Qmm") , size = 0.5) +
  geom_point(data=CNmodelnew, aes(x = date, y = Qpred, color = "Qpred"), size = 0.5) +
  ylab("Discharge, mm/d") +
  ggtitle(paste0("NNModel trained to date 2018-02-28 NSE =", round(NN_NSE_1,2), "_&_", "CNmodel NSE=",round(CNmodelnew_NSE,2))) +
  theme_classic()
dev.off()


#################      Predict Next 3 days Flow ##################
 # Using the net_2 NN model
###############################################################

dat_2=data.frame(Qm3ps=CNmodelnew$Qmm,
               Jday=strptime(CNmodelnew$date,"%Y-%m-%d")$yday+1,Precip=CNmodelnew$P,
               TMX=CNmodelnew$MaxTemp,TMN=CNmodelnew$MinTemp,
               jdate=julian(CNmodelnew$date))
# Need the previous days flow
# The NN for this example requires the previous days flow… remember this for 
# the homework questions. 
dat_2$Qm3ps_lag=dplyr::lag(dat_2$Qm3ps)
dat[is.na(dat_2)]=0
View(dat_2)
head(dat_2)

# Vector of today, tomorro, frid, sat, sun

yes = dat_2[4491,]
thurs = c(2.14, 110, 0.61, 25, 21, 19102, 2.14)
fri =  c(2.14, 111, 0.21, 28, 20, 19103, 2.14)
sat =  c(2.14, 112, 1.06, 27, 20, 19104, 2.14) 
sund =  c(2.14, 113, 0.76, 29, 19, 19105, 2.14)

forcast_df = rbind(thurs, fri, sat, sund)

colnames(forcast_df) = colnames(dat_2)
rownames(forcast_df) = c(4492, 4493, 4494, 4495)

forcast_df = rbind(dat_2, forcast_df)

forcast_df$Qm3ps_lag[1] =0

datmax <- apply(forcast_df, 2, max)
datmin <- apply(forcast_df, 2, min)
# now split the datat into train, validation and test sets.
scaled <- scale(forcast_df, center = datmin, scale = datmax - datmin)

# sub sample to get the test and train set
index_3 <- seq(1, nrow(forcast_df)-4)
train_index_3 <- sample(index_3, 0.8*length(index_3))
train_3 <- scaled[train_index_3,]
test_3 <- scaled[-train_index_3,]


NNForecast = predict(net, test_3[900:903, c("Qm3ps_lag","Precip","jdate","Jday","TMX","TMN")])
unscaled_NNForecast = (NNForecast * (max(dat_2$Qm3ps) - min(dat_2$Qm3ps))) + min(dat_2$Qm3ps)
rownames(unscaled_NNForecast) = c("thursday", "friday", "saturday", "sunday")
colnames(unscaled_NNForecast) = c("Qpred")

forcast_df = forcast_df[4492:4495,]
forcast_df$Qforecast = unscaled_NNForecast[,1]
rownames(forcast_df) = c("thursday", "friday", "saturday", "sunday")
