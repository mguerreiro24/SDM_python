library(dismo)
library(ncdf4)
library(dplyr)
library(leaflet)
library(htmlwidgets)


get_env_data <- function(DEPTH){
  temp <-raster(paste('new_temp',DEPTH,'.nc',sep=''))
  sal <-raster(paste('new_sal',DEPTH,'.nc',sep=''))
  o2 <-raster(paste('new_o2',DEPTH,'.nc',sep=''))
  pr <-raster(paste('new_pr',DEPTH,'.nc',sep=''))
  bd <-raster(paste('new_bd',DEPTH,'.nc',sep=''))
  ed <-raster(paste('new_ed',DEPTH,'.nc',sep=''))
  em <-raster(paste('new_em',DEPTH,'.nc',sep=''))
  #mud <-raster(paste('new_mud',DEPTH,'.nc',sep=''))
  musm <-raster(paste('new_musm',DEPTH,'.nc',sep=''))
  mumm <-raster(paste('new_mumm',DEPTH,'.nc',sep=''))
  #mld <-raster(paste('new_mld',DEPTH,'.nc',sep=''))
  mlsm <-raster(paste('new_mlsm',DEPTH,'.nc',sep=''))
  #mlmm <-raster(paste('new_mlmm',DEPTH,'.nc',sep=''))
  mlhmm <-raster(paste('new_mlhmm',DEPTH,'.nc',sep=''))
  #biomass <-raster(paste('new_biomass',DEPTH,'.nc',sep=''))
  
  env_data_to_map <- stack(temp,
                           sal,
                           o2,
                           pr,
                           bd,
                           ed,
                           em,
                           #mud,
                           musm,
                           mumm,
                           #mld,
                           mlsm,
                           #mlmm,
                           mlhmm
                           # biomass
  )
  return(env_data_to_map)
}

leafMap_prob <- function(MOD,name,add_title){
  points <- read.table("D:/PhD/articles/article SDM/data_occurrences.tsv",sep='\t',header = T)
  points <- na.omit(points)
  points <- points[points$Species==name,]
  palete <- colorNumeric(c('#00FF00','#FF0000'), values(MOD), na.color = 'transparent')#c(0,1)
  m <-  leaflet()%>% #leaflet(HR)
    addProviderTiles(provider = 'Esri.OceanBasemap') %>%
    addScaleBar()%>%
    addCircles(points$Longitude,points$Latitude, group='observations')%>%
    #addCircleMarkers(lng=~Longitude, lat=~Latitude,
    #                 clusterOptions = leaflet::markerClusterOptions())
    addRasterImage(colors = palete, MOD, opacity = .6,group=add_title)%>%
    addLegend(position = 'topright',pal = palete, values(MOD), title = paste(name,add_title,sep=' - '))%>%
    addLayersControl(overlayGroups=c(add_title,'observations'))#baseGroups=c(,,),
  saveWidget(m, file=paste(name,add_title,'.html',sep = ''))
  return(m)
}

leafMap <- function(MOD,name,add_title){
  points <- read.table("D:/PhD/articles/article SDM/data_occurrences.tsv",sep='\t',header = T)
  points <- na.omit(points)
  points <- points[points$Species==name,]
  #palete <- colorNumeric(c('#FF0000','#FFFFCC','#0000FF'), values(MOD), na.color = 'transparent')#values(MOD)
  palete <- colorNumeric(c('#0000FF','#FFFFCC','#FF0000'), values(MOD), na.color = 'transparent')#values(MOD)
  m <-  leaflet()%>% #leaflet(HR)
    addProviderTiles(provider = 'Esri.OceanBasemap') %>%
    addScaleBar()%>%
    addCircles(points$Longitude,points$Latitude)%>%
    #addCircleMarkers(lng=~Longitude, lat=~Latitude,
    #                 clusterOptions = leaflet::markerClusterOptions())
    addRasterImage(colors = palete, MOD, opacity = .6)%>%
    addLegend(position = 'topright',pal = palete, values(MOD), title = paste(name,add_title,sep=' - '))
  saveWidget(m, file=paste(name,add_title,'.html',sep = ''))
  return(m)
}

BuildingMaps <- function(MOD,name,deeper_sea,Weights,depths){
  MOD <- stack(paste(name,'.grd',sep = ''))#when data saved
  MODmax <- max(MOD, na.rm=T)
  MODmax <- mask(MODmax,deeper_sea)
  
  MODsd <- weighted.mean(MOD,Weights, na.rm=T)
  MODsd <- mask(MODsd,deeper_sea)
  
  MM <- which.max(MOD)
  MM <- subs(MM,data.frame(id=1:125, v=depths))
  MM <- mask(MM,deeper_sea)
  
  msd <- leafMap_prob(MODsd,name,'weighted Mean')
  mmax <- leafMap_prob(MODmax,name,'Max')
  mm <- leafMap(MM,name,'Depth')
}

into_the_DEEP <- function(depths,d,bio){#,maxent_model
for (i in 0:124){
  cat(i,'\n')
  
deeper_sea <-raster(paste('new_bd',i,'.nc',sep=''))
deeper_sea[deeper_sea<0] <- NA

env_data_to_map <- get_env_data(i)
env_data_to_map <- mask(env_data_to_map,deeper_sea)
cat('envData sorted, ')

predictionsb <- predict(env_data_to_map, bio, type="response")
#plot(predictionsb)
cat('bioclim done, ')

predictionsd <-predict(env_data_to_map, d, type="response")
#plot(predictionsd)
cat('dismo done, ')

predictionsm <-predict(env_data_to_map, maxent_model, type="response")
#plot(predictionsm)
cat('MaxEnt done.')


mods <- stack(predictionsb,predictionsd,predictionsm)
mods <- mean(mods, na.rm = T)
#plot(mods)
if (i==0){mod <- c(mods)}else {mod <- c(mod,mods)}
cat('\n')
}
MOD <- brick(mod)
names(MOD) <- depths
  return(MOD)
}

get_depth <- function(value,depths){
  #test[3] <- get_depth(test[,3],depths)
  f <-Biobase::matchpt(value,depths)[,1]
  return(f-1)
}

set_depth <- function(value,depths){
  new_value <- depths[value]
  return(new_value)
}

check_corr <- function(dat){
  data <- na.omit(dat)
  for (i in 1:(dim(data)[2]-1)){
    for (ii in (i+1):dim(data)[2]){
      cat(names(data)[i],' vs. ',names(data)[ii],': ',cor(data[i], data[ii], method = "pearson"))
      cat('\n')
    }
  }
      
}

####################################################################################################################################
setwd("D:/PhD/articles/article SDM/data/")

depths <- c(0.49402499198913574,0.5057600140571594,1.5413750410079956,1.5558550357818604,2.6456689834594727,2.667681932449341,3.8194949626922607,3.8562800884246826,5.078224182128906,5.1403608322143555,6.440614223480225,6.543034076690674,7.92956018447876,8.09251880645752,9.572997093200684,9.822750091552734,11.404999732971191,11.773679733276367,13.467140197753906,13.991040229797363,15.810070037841797,16.525320053100586,18.495559692382812,19.429800033569336,21.598819732666016,22.757619857788086,25.211410522460938,26.558300018310547,29.444730758666992,30.87455940246582,34.43415069580078,35.74020004272461,40.344051361083984,41.18001937866211,47.211891174316406,47.37369155883789,53.85063934326172,55.76428985595703,61.11283874511719,65.80726623535156,69.02168273925781,77.61116027832031,77.85385131835938,86.92942810058594,92.3260726928711,97.04131317138672,108.0302963256836,109.72930145263672,120.0,130.66600036621094,133.0758056640625,147.4062042236328,155.85069274902344,163.1645050048828,180.54989624023438,186.12559509277344,199.7899932861328,221.14120483398438,222.47520446777344,244.89059448242188,266.0403137207031,271.3564147949219,300.88751220703125,318.1274108886719,333.86279296875,370.6885070800781,380.2130126953125,411.7939147949219,453.9377136230469,457.6256103515625,508.639892578125,541.0889282226562,565.2922973632812,628.0260009765625,643.5667724609375,697.2587280273438,763.3331298828125,773.3682861328125,856.6790161132812,902.3392944335938,947.4478759765625,1045.85400390625,1062.43994140625,1151.990966796875,1245.291015625,1265.8609619140625,1387.376953125,1452.2509765625,1516.364013671875,1652.5679931640625,1684.2840576171875,1795.6710205078125,1941.8929443359375,1945.2960205078125,2101.027099609375,2225.077880859375,2262.422119140625,2429.02490234375,2533.3359375,2600.3798828125,2776.0390625,2865.702880859375,2955.570068359375,3138.56494140625,3220.820068359375,3324.64111328125,3513.446044921875,3597.031982421875,3704.656982421875,3897.98193359375,3992.48388671875,4093.158935546875,4289.953125,4405.22412109375,4488.15478515625,4687.5810546875,4833.291015625,4888.06982421875,5089.47900390625,5274.7841796875,5291.68310546875,5494.5751953125,5698.06103515625,5727.9169921875,5902.05810546875)
Weights <- depths
Weights[2:125] <- Weights[2:125]-Weights[1:124]

deeper_sea <-raster('new_bd0.nc')
deeper_sea[deeper_sea<200] <- NA

#var <-read.table("Ancistroteuthis_new.tsv",sep='\t',header = T)
#name <- 'Ancistroteuthis lichtensteinii'
#absences_pseudo<-read.table('background_points_Ancistroteuthis_lichtensteinii_new.tsv',sep='\t',header = T)
#var<-read.table("Tsagittatus_new.tsv",sep='\t',header = T)
#absences_pseudo<-read.table('background_points_Todarodes_sagittatus_new.tsv',sep='\t',header = T)
#name <- 'Todarodes sagittatus'
#var<-read.table("Averanyi_new.tsv",sep='\t',header = T)
#name <- 'Abralia veranyi'
#absences_pseudo<-read.table('background_points_Abralia_veranyi_new.tsv',sep='\t',header = T)
#var<-read.table("Cscabra_new.tsv",sep='\t',header = T)
#name<-'Cranchia scabra'
#absences_pseudo<-read.table('background_points_Cranchia_scabra_new.tsv',sep='\t',header = T)
#var<-read.table("Hbonnellii_new.tsv",sep='\t',header = T)
#name <-'Histioteuthis bonnellii'
#absences_pseudo<-read.table('background_points_Histioteuthis_bonnellii_new.tsv',sep='\t',header = T)
#var<-read.table("Hreversa_new.tsv",sep='\t',header = T)
#name <-'Histioteuthis reversa'
#absences_pseudo<-read.table('background_points_Histioteuthis_reversa_new.tsv',sep='\t',header = T)

var<-read.table("Obanksii_d.tsv",sep='\t',header = T)
name <-'Onychoteuthis banksii'
absences_pseudo<-read.table('background_points_Cranchia_scabra_new.tsv',sep='\t',header = T)


##ONly biomass
# var <- var[sample(1:nrow(var)), c(4:8,18)]#shuffle
# absences_pseudo <- absences_pseudo[,c(4:8,18)]

##all exceptbiomass
#var <- var[sample(1:nrow(var)), c(4:17)]#shuffle
#absences_pseudo <- absences_pseudo[,c(4:17)]

##all except biomass and depth of mesopelagics (correlated to epipelagic depth)
##and mesopelagic_l_migratory_mass(correlated with mesopelagic_u_static_mass)
var <- var[sample(1:nrow(var)), c(4:10,12,13,15,17)]#,18)]#shuffle
absences_pseudo <- absences_pseudo[,c(4:10,12,13,15,17)]


check_corr(var)


data1 <- sort(sample(nrow(var), nrow(var)*.7))
data2 <- sort(sample(nrow(absences_pseudo), nrow(absences_pseudo)*.7))
#test/train
absences_test <- absences_pseudo[-data2,]
presences_test <- var[-data1,]

absences_train <- absences_pseudo[data2,]
train <- var[data1,]#train
#var <- unique(var)
#var <- na.omit(var)

Presences <- train
Absences <- absences_train
Presences$present <- 1
Absences$present <- 0
presabs <- rbind(Presences,Absences)

maxent_model <- maxent(x=presabs[,1:11],
                       p=presabs$present)

d <- domain(train)
bio <- bioclim(train)



evaluation <- evaluate(presences_test,absences_test,d)
plot(evaluation,'ROC')
evaluation <- evaluate(presences_test,absences_test,bio)
plot(evaluation,'ROC')
evaluation <- evaluate(presences_test,absences_test,maxent_model)
plot(evaluation,'ROC')




MOD <- into_the_DEEP(depths,d,bio)#,maxent_model
writeRaster(MOD,filename=paste(name,'.grd',sep = ''), bandorder='BIL', overwrite=TRUE)


#MOD <- stack(paste(name,'.grd',sep = ''))#when data saved
MODmax <- max(MOD, na.rm=T)
maskMAX <- MODmax
maskMAX[maskMAX<0.25] <- NA

MODsd <- weighted.mean(MOD,Weights, na.rm=T)
maskWM <- MODsd
maskWM[maskWM<.2] <- NA

MM <- which.max(MOD)
MM <- subs(MM,data.frame(id=1:125, v=depths))
maskMM <- MM
maskMM[maskMM>3000] <- NA

MODsd <- mask(MODsd,deeper_sea)
MODsd <- mask(MODsd,maskMAX)
MODsd <- mask(MODsd,maskWM)
MODsd <- mask(MODsd,maskMM)
msd <- leafMap_prob(MODsd,name,'weighted Mean')
msd

MODmax <- mask(MODmax,deeper_sea)
MODmax <- mask(MODmax,maskMAX)
MODmax <- mask(MODsd,maskWM)
MODmax <- mask(MODsd,maskMM)
mmax <- leafMap_prob(MODmax,name,'Max')
mmax

MM <- mask(MM,deeper_sea)
MM <- mask(MM,maskMAX)
MM <- mask(MM,maskWM)
MM <- mask(MM,maskMM)
mm <- leafMap(MM,name,'Depth')
mm
#################################################################################

BuildingMaps(MOD,'Ancistroteuthis lichtensteinii',deeper_sea,Weights,depths)
BuildingMaps(MOD,'Todarodes sagittatus',deeper_sea,Weights,depths)
BuildingMaps(MOD,'Abralia veranyi',deeper_sea,Weights,depths)
BuildingMaps(MOD,'Cranchia scabra',deeper_sea,Weights,depths)
BuildingMaps(MOD,'Histioteuthis bonnellii',deeper_sea,Weights,depths)
BuildingMaps(MOD,'Histioteuthis reversa',deeper_sea,Weights,depths)

#################################################################################

var<-read.table("Obanksii_d.tsv",sep='\t',header = T)
var<-read.table("Cscabra_d.tsv",sep='\t',header = T)
absences_pseudo<-read.table('background_points_Cranchia_scabra_new.tsv',sep='\t',header = T)
name <-'Onychoteuthis banksii'
MOD <- stack(paste(name,'.grd',sep = ''))#when data saved
var <- var[, c(4:10,12,13,15,17,18)]
absences_pseudo <- absences_pseudo[,c(4:10,12,13,15,17)]
absences_pseudo$density <- 0
var <- na.omit(var)
absences_pseudo <- absences_pseudo[1:nrow(var),]
density_a <- rbind(var,absences_pseudo)


density_a$md <- predict(object=d,x=density_a[,1:11],type='response')
density_a$mb <- predict(object=bio,x=density_a[,1:11],type='response')
density_a$mM <- predict(object=maxent_model,x=density_a[,1:11],type='response')
density_a$m <- rowMeans(density_a[,13:15])
#l1 <- loess(density_a$density ~ density_a$m,span=.1)
#l4 <- loess(density_a$density ~ density_a$m,span=.25)
#l2 <- loess(density_a$density ~ density_a$m,span=.5)
LM <- lm(density_a$density[1:13] ~ density_a$m[1:13])
#LM <- lm(density_a$density ~ density_a$m)
#mk <- (MOD*3.037e-06)-1.455e-07
mk <- ((MOD*2.837e-06)-4.950e-08) * (111000*0.08333333)^2#km^2
mk <- mk*50*.5*.06#grams, 50% deposition, .06 carbon

MK <- which.max(mk)
MK <- subs(MK,data.frame(id=1:125, v=depths))
MK <- mask(MK,deeper_sea)
MKK <- MK
MKK[MKK>250] <-NA
MK <- mask(MK,MKK)
KK <- leafMap(MK,name,'Max density depth')
KK


mk <- weighted.mean(mk,Weights, na.rm=T)
deeper_sea <-raster('new_bd0.nc')
mk <- mk * deeper_sea
deeper_sea[deeper_sea<200] <- NA
#mk <- max(mk, na.rm=T)
mk <- mask(mk,deeper_sea)
mk <- mask(mk,MKK)
cellStats(mk,sum)
mkk <- mk
mkk[mkk<150000] <- NA
mk <- mask(mk,mkk)
kk <- leafMap(mk,name,'gCarbon')
kk





plot(mk[1])
hist(mk[1])
########################
library(gam)

setwd("D:/PhD/articles/article SDM/data/")

depths <- c(0.49402499198913574,0.5057600140571594,1.5413750410079956,1.5558550357818604,2.6456689834594727,2.667681932449341,3.8194949626922607,3.8562800884246826,5.078224182128906,5.1403608322143555,6.440614223480225,6.543034076690674,7.92956018447876,8.09251880645752,9.572997093200684,9.822750091552734,11.404999732971191,11.773679733276367,13.467140197753906,13.991040229797363,15.810070037841797,16.525320053100586,18.495559692382812,19.429800033569336,21.598819732666016,22.757619857788086,25.211410522460938,26.558300018310547,29.444730758666992,30.87455940246582,34.43415069580078,35.74020004272461,40.344051361083984,41.18001937866211,47.211891174316406,47.37369155883789,53.85063934326172,55.76428985595703,61.11283874511719,65.80726623535156,69.02168273925781,77.61116027832031,77.85385131835938,86.92942810058594,92.3260726928711,97.04131317138672,108.0302963256836,109.72930145263672,120.0,130.66600036621094,133.0758056640625,147.4062042236328,155.85069274902344,163.1645050048828,180.54989624023438,186.12559509277344,199.7899932861328,221.14120483398438,222.47520446777344,244.89059448242188,266.0403137207031,271.3564147949219,300.88751220703125,318.1274108886719,333.86279296875,370.6885070800781,380.2130126953125,411.7939147949219,453.9377136230469,457.6256103515625,508.639892578125,541.0889282226562,565.2922973632812,628.0260009765625,643.5667724609375,697.2587280273438,763.3331298828125,773.3682861328125,856.6790161132812,902.3392944335938,947.4478759765625,1045.85400390625,1062.43994140625,1151.990966796875,1245.291015625,1265.8609619140625,1387.376953125,1452.2509765625,1516.364013671875,1652.5679931640625,1684.2840576171875,1795.6710205078125,1941.8929443359375,1945.2960205078125,2101.027099609375,2225.077880859375,2262.422119140625,2429.02490234375,2533.3359375,2600.3798828125,2776.0390625,2865.702880859375,2955.570068359375,3138.56494140625,3220.820068359375,3324.64111328125,3513.446044921875,3597.031982421875,3704.656982421875,3897.98193359375,3992.48388671875,4093.158935546875,4289.953125,4405.22412109375,4488.15478515625,4687.5810546875,4833.291015625,4888.06982421875,5089.47900390625,5274.7841796875,5291.68310546875,5494.5751953125,5698.06103515625,5727.9169921875,5902.05810546875)
Weights <- depths
Weights[2:125] <- Weights[2:125]-Weights[1:124]

deeper_sea <-raster('new_bd0.nc')
deeper_sea[deeper_sea<200] <- NA

var<-read.table("Obanksii_d.tsv",sep='\t',header = T)
var <- var[, c(4:10,12,13,15,17,18)]#shuffle
name <-'Onychoteuthis banksii'
#name <-'Octopoteuthis sicula'
absences_pseudo<-read.table('background_points_Cranchia_scabra_new.tsv',sep='\t',header = T)
absences_pseudo <- absences_pseudo[1:12,c(4:10,12,13,15,17)]
absences_pseudo$density <- 0
data1 <- sort(sample(nrow(var), nrow(var)*.8))
data2 <- sort(sample(nrow(absences_pseudo), nrow(absences_pseudo)*.8))
#test/train
absences_test <- absences_pseudo[-data2,]
presences_test <- var[-data1,]

absences_train <- absences_pseudo[data2,]
train <- var[data1,]#train


presabs <- rbind(train,absences_train)


GAM <- gam(density ~ temperature + salinity + O2 + pressure + benthos_distance + epipelagic_depth + epipelagic_mass + mesopelagic_u_static_mass + mesopelagic_l_static_mass + mesopelagic_u_migratory_mass + mesopelagic_l_HighlyMigratory_mass,#s(temperature) + s(salinity) + s(O2) + s(pressure) + s(benthos_distance) + s(epipelagic_depth) + s(epipelagic_mass) + s(mesopelagic_u_static_mass) + s(mesopelagic_l_static_mass) + s(mesopelagic_u_migratory_mass) + s(mesopelagic_l_HighlyMigratory_mass),
           #method = 'REML',
           family = binomial(link = "logit"),
           data = presabs)
logistic_model <- glm(density ~ temperature + 
                        salinity + O2 + pressure + 
                        benthos_distance + 
                        epipelagic_depth + epipelagic_mass +
                        mesopelagic_u_static_mass + mesopelagic_u_migratory_mass +
                        mesopelagic_l_static_mass + mesopelagic_l_HighlyMigratory_mass, family = binomial(link = "logit"), data = presabs)


evaluation <- evaluate(presences_test,absences_test,logistic_model)
plot(evaluation,'ROC')

plot(evaluate(presences_test,absences_test,GAM),'ROC')


for (i in 0:124){
  cat(i,'\n')
  
  deeper_sea <-raster(paste('new_bd',i,'.nc',sep=''))
  deeper_sea[deeper_sea<0] <- NA
  
  env_data_to_map <- get_env_data(i)
  env_data_to_map <- mask(env_data_to_map,deeper_sea)
  cat('envData sorted, ')
  
  predictionsb <- predict(object=logistic_model,x=env_data_to_map, type="response")
  #plot(predictionsb)
  cat('glm done, ')
  
  if (i==0){mod <- c(predictionsb)}else {mod <- c(mod,predictionsb)}
  cat('\n')
}
MOD <- brick(mod)
names(MOD) <- depths









var <- na.omit(var)
logistic_model <- glm(present ~ temperature + salinity + 
                        O2 + pressure + benthos_distance +
                        Biomass,
                      family = binomial(link = "logit"),
                      data = presabs)

presences <- filter(var,present == 1)
absences <- filter(var,present == 0)
#background <-
#evaluation <- evaluate(presences,absences,bio)
evaluation <- evaluate(presences,absences,logistic_model)
plot(evaluation,'ROC')
predictions <- predict(env_data_to_map,
                       logistic_model,
                       type="response")
plot(predictions)
####################################




  #addCircles(lng = ,lat = )%>%
  #addRasterImage(colors = palete, f, opacity = .6)%>%
  #addLegend(position = 'topright',pal = palete, values(f), title = 'O2')%>%
  #addMarkers(lat = 48, lng = 16, label = 'Wien', group = 'name to keep in legend')
  #addLayersControl(baseGroups=c(,,),overlayGroups=c(,,))
  #
#m
hist(TE$Bathimetry,main='Todaropsis eblanae')
#l

setwd("D:/PhD/articles/article SDM/")
points <- read.table("D:/PhD/articles/article SDM/data_occurrences.tsv",sep='\t',header = T)
points <- na.omit(points)
AV <- points[points$Species=='Abralia veranyi',]
AL <- points[points$Species=='Ancistroteuthis lichtensteinii',]
HB <- points[points$Species=='Histioteuthis bonnellii',]
HR <- points[points$Species=='Histioteuthis reversa',]
TS <- points[points$Species=='Todarodes sagittatus',]
TE <- points[points$Species=='Todaropsis eblanae',]
GF <- points[points$Species=='Gonatus fabricii',]
#webshot::install_phantomjs()

