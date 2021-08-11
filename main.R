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
  # ed <-raster(paste('new_ed',DEPTH,'.nc',sep=''))
  # em <-raster(paste('new_em',DEPTH,'.nc',sep=''))
  # mud <-raster(paste('new_mud',DEPTH,'.nc',sep=''))
  # musm <-raster(paste('new_musm',DEPTH,'.nc',sep=''))
  # mumm <-raster(paste('new_mumm',DEPTH,'.nc',sep=''))
  # mld <-raster(paste('new_mld',DEPTH,'.nc',sep=''))
  # mlsm <-raster(paste('new_mlsm',DEPTH,'.nc',sep=''))
  # mlmm <-raster(paste('new_mlmm',DEPTH,'.nc',sep=''))
  # mlhmm <-raster(paste('new_mlhmm',DEPTH,'.nc',sep=''))
  biomass <-raster(paste('new_biomass',DEPTH,'.nc',sep=''))
  
  env_data_to_map <- stack(temp,
                           sal,
                           o2,
                           pr,
                           bd,
                           # ed,
                           # em,
                           # mud,
                           # musm,
                           # mumm,
                           # mld,
                           # mlsm,
                           # mlmm,
                           # mlhmm
                           biomass
  )
  return(env_data_to_map)
}

leafMap <- function(MOD,name){
  points <- read.table("D:/PhD/articles/article SDM/data_occurrences.tsv",sep='\t',header = T)
  points <- na.omit(points)
  points <- points[points$Species==name,]
  palete <- colorNumeric(c('#0000FF','#FFFFCC','#FF0000'), values(MOD), na.color = 'transparent')
  m <-  leaflet()%>% #leaflet(HR)
    addProviderTiles(provider = 'Esri.OceanBasemap') %>%
    addScaleBar()%>%
    addCircles(points$Longitude,points$Latitude)%>%
    #addCircleMarkers(lng=~Longitude, lat=~Latitude,
    #                 clusterOptions = leaflet::markerClusterOptions())
    addRasterImage(colors = palete, MOD, opacity = .6)%>%
    addLegend(position = 'topright',pal = palete, values(MOD), title = name)
  saveWidget(m, file=paste(name,'.html',sep = ''))
  return(m)
}

into_the_DEEP <- function(depths,d,bio,deeper_sea,test){
for (i in 0:124){
  cat(i,'\n')
  

env_data_to_map <- get_env_data(i)
env_data_to_map <- mask(env_data_to_map,deeper_sea)
cat('envData sorted, ')

predictionsb <- predict(env_data_to_map,
                       bio,
                       type="response")
#plot(predictionsb)
cat('bioclim done, ')
predictionsd <-predict(env_data_to_map,
                      d,
                      type="response")
#plot(predictionsd)
cat('dismo done.')


mods <- stack(predictionsb,predictionsd)
mods <- mean(mods, na.rm = T)
#plot(mods)
if (i==0){mod <- c(mods)}else {mod <- c(mod,mods)}
cat('\n')
}
MOD <- brick(mod)
names(MOD) <- depths
  return(list(MOD=MOD,evals_p=evals_p))
}

get_depth <- function(value,depths){
  #test[3] <- get_depth(test[,3],depths)
  f <-Biobase::matchpt(value,depths)[,1]
  return(f-1)
}

####################################################################################################################################
setwd("D:/PhD/articles/article SDM/data/")

depths <- c(0.49402499198913574,0.5057600140571594,1.5413750410079956,1.5558550357818604,2.6456689834594727,2.667681932449341,3.8194949626922607,3.8562800884246826,5.078224182128906,5.1403608322143555,6.440614223480225,6.543034076690674,7.92956018447876,8.09251880645752,9.572997093200684,9.822750091552734,11.404999732971191,11.773679733276367,13.467140197753906,13.991040229797363,15.810070037841797,16.525320053100586,18.495559692382812,19.429800033569336,21.598819732666016,22.757619857788086,25.211410522460938,26.558300018310547,29.444730758666992,30.87455940246582,34.43415069580078,35.74020004272461,40.344051361083984,41.18001937866211,47.211891174316406,47.37369155883789,53.85063934326172,55.76428985595703,61.11283874511719,65.80726623535156,69.02168273925781,77.61116027832031,77.85385131835938,86.92942810058594,92.3260726928711,97.04131317138672,108.0302963256836,109.72930145263672,120.0,130.66600036621094,133.0758056640625,147.4062042236328,155.85069274902344,163.1645050048828,180.54989624023438,186.12559509277344,199.7899932861328,221.14120483398438,222.47520446777344,244.89059448242188,266.0403137207031,271.3564147949219,300.88751220703125,318.1274108886719,333.86279296875,370.6885070800781,380.2130126953125,411.7939147949219,453.9377136230469,457.6256103515625,508.639892578125,541.0889282226562,565.2922973632812,628.0260009765625,643.5667724609375,697.2587280273438,763.3331298828125,773.3682861328125,856.6790161132812,902.3392944335938,947.4478759765625,1045.85400390625,1062.43994140625,1151.990966796875,1245.291015625,1265.8609619140625,1387.376953125,1452.2509765625,1516.364013671875,1652.5679931640625,1684.2840576171875,1795.6710205078125,1941.8929443359375,1945.2960205078125,2101.027099609375,2225.077880859375,2262.422119140625,2429.02490234375,2533.3359375,2600.3798828125,2776.0390625,2865.702880859375,2955.570068359375,3138.56494140625,3220.820068359375,3324.64111328125,3513.446044921875,3597.031982421875,3704.656982421875,3897.98193359375,3992.48388671875,4093.158935546875,4289.953125,4405.22412109375,4488.15478515625,4687.5810546875,4833.291015625,4888.06982421875,5089.47900390625,5274.7841796875,5291.68310546875,5494.5751953125,5698.06103515625,5727.9169921875,5902.05810546875)
Weights <- depths
Weights[2:125] <- Weights[2:125]-Weights[1:124]

#var <-read.table("Ancistroteuthis.tsv",sep='\t',header = T)
#name <- 'Ancistroteuthis lichtensteinii'
var<-read.table("Tsagittatus_new.tsv",sep='\t',header = T)
absences_pseudo<-read.table('background_points_Todarodes_sagittatus_new.tsv',sep='\t',header = T)
name <- 'Todarodes sagittatus'
#var<-read.table("Averanyi.tsv",sep='\t',header = T)
#name <- 'Abralia veranyi'
#var<-read.table("Cscabra_new.tsv",sep='\t',header = T)
#name<-'Cranchia scabra'
#var<-read.table("Hbonnellii.tsv",sep='\t',header = T)
#name <-'Histioteuthis bonnellii'
#var<-read.table("Hreversa.tsv",sep='\t',header = T)
#name <-'Histioteuthis reversa'


var <- var[sample(1:nrow(var)), c(4:8,18)]
absences_pseudo <- absences_pseudo[,c(4:8,18)]

data1 <- sort(sample(nrow(var), nrow(var)*.7))
data2 <- sort(sample(nrow(absences_pseudo), nrow(absences_pseudo)*.7))
#test/train
absences_test <- absences_pseudo[-data2,]
presences_test <- var[-data1,]

absences_train <- absences_pseudo[data2,]
var <- var[data1,]#train
var <- unique(var)
#var <- na.omit(var)

Presences <- var
Absences <- absences_train
Presences$present <- 1
Absences$present <- 0
presabs <- rbind(Presences,Absences)

logistic_model <- glm(present ~ temperature + salinity + 
                        O2 + pressure + benthos_distance +
                        Biomass,
                      family = binomial(link = "logit"),
                      data = presabs)

maxent_model <- maxent(x=presabs[,1:6],
                       p=presabs$present)

d <- domain(var)
bio <- bioclim(var)






presences <- test
evaluation <- evaluate(presences,absences,d)
plot(evaluation,'ROC')
evaluation <- evaluate(presences,absences,bio)
plot(evaluation,'ROC')




deeper_sea <-raster('new_bd0.nc')
deeper_sea[deeper_sea<200] <- NA

M <- into_the_DEEP(depths,d,bio,deeper_sea,test)
MOD <- M$MOD
evals_p <- M$evals_p
writeRaster(MOD,filename=paste(name,'.grd',sep = ''), bandorder='BIL', overwrite=TRUE)
MOD <- weighted.mean(MOD,Weights, na.rm=T)
#MOD <- max(MOD, na.rm=T)
MOD <- mask(MOD,deeper_sea)
#plot(MOD)
m <- leafMap(MOD,name)

m





########################
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
  #addMarkers(lat = 48, lng = 16, label = 'Wien')
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

