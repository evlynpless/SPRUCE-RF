#R script for RF Model using mRates and training/testing data

library("randomForestSRC")
par(mfrow=c(1,1))

#Read in deme level values from all spatial variables
Altitude <- read.table("altitude_1KMmedian_MERIT_MAPS1.txt")
Slope <- read.table("slope_1KMmedian_MERIT_MAPS1.txt")
Aridity <- read.table("AI_annual_MAPS1.txt")
MeanTemp <- read.table("bio1_mean_MAPS1.txt")
MaxTemp <- read.table("bio5_mean_MAPS1.txt")
MinTemp <- read.table("bio6_mean_MAPS1.txt")
MeanPrec <- read.table("bio12_mean_MAPS1.txt")
PrecWet <- read.table("bio13_mean_MAPS1.txt")
PrecDry <- read.table("bio14_mean_MAPS1.txt")
GPP <- read.table("GPP_mean_MAPS1.txt")
Fusc <- read.table("fuscgroup_MAPS1.txt")
Mors <- read.table("morsgroup_MAPS1.txt")
Palp <- read.table("palpgroup_MAPS1.txt")
LangAA <- read.table("LanguageMap_AA_MAPS1.txt")
LangNC <- read.table("LanguageMap_NC_MAPS1.txt")
LangNS <- read.table("LanguageMap_NS_MAPS1.txt")
Evergreen <- read.table("consensus_full_class_2_MAPS1.txt")
Decid <- read.table("consensus_full_class_3_MAPS1.txt")
Tree <- read.table("consensus_full_class_4_MAPS1.txt")
Shrub <- read.table("consensus_full_class_5_MAPS1.txt")
Herb <- read.table("consensus_full_class_6_MAPS1.txt")
Crop <- read.table("consensus_full_class_7_MAPS1.txt")
Barren <- read.table("consensus_full_class_11_MAPS1.txt")
Water <- read.table("consensus_full_class_12_MAPS1.txt")
Kernel <- read.table("kernel.txt")
Access <- read.table("accessibility_to_cities_2015_v1.0_res_MAPS1.txt")
Building <- read.table("GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_MAPS1.txt")
Friction <- read.table("friction_surface_2015_v1.0_MAPS1.txt")

#Merge the environmental datasets into one dataframe: each deme has own row, and each variable has own column
Env.Table <- data.frame(Altitude, Slope, Aridity, MeanTemp, MaxTemp, MinTemp, MeanPrec, PrecWet, PrecDry, GPP, Fusc, Mors, Palp, LangAA, LangNC, LangNS, Evergreen, Decid, Tree, Shrub, Herb, Crop, Barren, Water, Kernel, Access, Building, Friction)
names(Env.Table) <- c("Altitude", "Slope", "Aridity", "MeanTemp", "MaxTemp", "MinTemp", "MeanPrec", "PrecWet", "PrecDry", "GPP", "Fusc", "Mors", "Palp", "LangAA", "LangNC", "LangNS", "Evergreen", "Decid", "Tree", "Shrub", "Herb", "Crop", "Barren", "Water", "Kernel", "Access", "Building", "Friction")

#Read in deme level values for migration surfaces (each column in the mRates.txt file is a deme, and the order of demes is specified by demes.txt) for three independent runs of MAPS
#(Make sure demes are in same order as the environmental data)
#Convert to migration rate and take log value to normalize response variable
mRates1 <- read.table("2-4cM_5million_Run1/mRates.txt", header = F)
mRates1 <- 10^(mRates1)
mRates1 <- log(mRates1)
mRates1 <- colMeans(mRates1)

mRates2 <- read.table("2-4cM_5million_Run2/mRates.txt", header = F)
mRates2 <- 10^(mRates2)
mRates2 <- log(mRates2)
mRates2 <- colMeans(mRates2)

mRates3 <- read.table("2-4cM_5million_Run3/mRates.txt", header = F)
mRates3 <- 10^(mRates3)
mRates3 <- log(mRates3)
mRates3 <- colMeans(mRates3)

mRates1_2_3 <- cbind(mRates1, mRates2, mRates3)
mRatesAvg <- rowMeans(mRates1_2_3 )
Full.Table = cbind(Env.Table, mRatesAvg)

#"Full Model" random forest with tuning
#FullModel_noTuning_option = rfsrc(mRatesAvg ~ Altitude + Slope + Aridity + MeanTemp + MaxTemp + MinTemp + MeanPrec + PrecWet + PrecDry + GPP + Fusc + Mors + Palp +  LangAA + LangNC + LangNS + Evergreen + Decid + Tree + Shrub + Herb + Crop + Barren + Water + Kernel, importance=TRUE, na.action=c("na.omit"), data=Full.Table)
FullModel_tune = tune(mRatesAvg ~ Altitude + Slope + Aridity + MeanTemp + MaxTemp + MinTemp + MeanPrec + PrecWet + PrecDry + GPP + Fusc + Mors + Palp +  LangAA + LangNC + LangNS + Evergreen + Decid + Tree + Shrub + Herb + Crop + Barren + Water + Kernel + Access + Building + Friction, importance=TRUE, na.action=c("na.omit"), data=Full.Table)
FullModel_tune$optimal[["mtry"]]
FullModel_tune$optimal[["nodesize"]]
FullModel_withTuning = rfsrc(mRatesAvg ~ Altitude + Slope + Aridity + MeanTemp + MaxTemp + MinTemp + MeanPrec + PrecWet + PrecDry + GPP + Fusc + Mors + Palp +  LangAA + LangNC + LangNS + Evergreen + Decid + Tree + Shrub + Herb + Crop + Barren + Water + Kernel + Access + Building + Friction, importance=TRUE, na.action=c("na.omit"), mtry = FullModel_tune$optimal[["mtry"]], nodesize =  FullModel_tune$optimal[["nodesize"]], data=Full.Table)

#Find correlation between predicted and actual migration rates, and RMSE of the model
R = cor(FullModel_withTuning$predicted.oob, Full.Table$mRatesAvg)
RMSE = sqrt(mean((FullModel_withTuning$predicted.oob - Full.Table$mRatesAvg)^2))

#Print results and create figures
FullModel_withTuning
plot(FullModel_withTuning)
paste("R", R)
paste("RMSE", RMSE)

pdf("mRates_2-4_varImp.pdf", 7, 7)
plot(FullModel_withTuning, m.target = NULL, plots.one.page = TRUE, sorted = TRUE, verbose = TRUE)
dev.off() 

pdf("mRates_2-4_FullScatter.pdf", 5, 5)
plot(FullModel_withTuning$predicted.oob, Full.Table$mRatesAvg, xlab="Predicted migration", ylab="Observed migration (MAPS)")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(R,3))), cex=0.7)
dev.off()


#Run random forest again with 10-fold cross-validation

#Randomly shuffle the data
yourData<-Full.Table[sample(nrow(Full.Table)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(yourData)),breaks=10,labels=FALSE)

#Perform 10-fold cross-validation
for(i in 1:10){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- yourData[testIndexes, ]
  assign(paste0("testData_", i), testData)
  trainData <- yourData[-testIndexes, ]
  assign(paste0("trainData_", i), trainData)
  #write.csv(trainData, paste0("/Users/eviepless/Documents/Postdoc/Postdoc_research/EasternAfrica/MAPS_FS_structure/Pagani_Scheinfeldt/ExcludeNorthernSudan/RF_Out/10_fold_CV/6_Inf_trainData_", i, ".csv"))
}

for(i in 1:10) {
  Full.Table.train = assign(paste0("trainData_",i), get(paste0("trainData_", i)))
  Full.Table.valid = assign(paste0("testData_",i), get(paste0("testData_", i)))
  TrainModel_tune = tune(mRatesAvg ~ Altitude + Slope + Aridity + MeanTemp + MaxTemp + MinTemp + MeanPrec + PrecWet + PrecDry + GPP + Fusc + Mors + Palp +  LangAA + LangNC + LangNS + Evergreen + Decid + Tree + Shrub + Herb + Crop + Barren + Water, importance=TRUE, na.action=c("na.omit"), data=Full.Table.train)
  TrainModel_tune$optimal[["mtry"]]
  TrainModel_tune$optimal[["nodesize"]]
  TrainModel = rfsrc(mRatesAvg ~ Altitude + Slope + Aridity + MeanTemp + MaxTemp + MinTemp + MeanPrec + PrecWet + PrecDry + GPP + Fusc + Mors + Palp +  LangAA + LangNC + LangNS + Evergreen + Decid + Tree + Shrub + Herb + Crop + Barren + Water, importance=TRUE, na.action=c("na.omit"), mtry = TrainModel_tune$optimal[["mtry"]], nodesize =  TrainModel_tune$optimal[["nodesize"]], data=Full.Table.train)
  
  #TrainModel
  #plot(TrainModel)
  #results_df <- data.frame("Rtrain", "Rtest", "RMSEtrain", "RMSEtest")
  Rtrain = cor(TrainModel$predicted.oob, Full.Table.train$mRatesAvg)
  Rtest = cor((predict.rfsrc(TrainModel, Full.Table.valid))$predicted, Full.Table.valid$mRatesAvg)
  RMSEtrain = sqrt(mean((TrainModel$predicted.oob - Full.Table.train$mRatesAvg)^2))
  RMSEtest = sqrt(mean((predict.rfsrc(TrainModel, Full.Table.valid)$predicted - Full.Table.valid$mRatesAvg)^2))
  resultsList <- c(Rtrain, Rtest, RMSEtrain, RMSEtest)
  print(resultsList)
  print(TrainModel)
}


