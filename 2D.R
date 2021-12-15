library(dismo)
library(ncdf4)
library(dplyr)
library(raster)
library(leaflet)
library(htmlwidgets)
library(sdmpredictors)
library(ENMTools)
library(nnet)#neural networks
library(gbm)#BRT
library(mgcv)#GAM GARP?



setwd("D:/PhD/articles/article SDM/")
#list_layers( datasets="Bio-ORACLE" )$name

spThin <- function(xy,dist,rep=1) {
  #dist -> distance between points in kilometers
  #rep -> replicates
  .r <- raster(resolution=0.008333333) # empty raster with 1km resolution
  .a <- raster(ext=extent(.r),resolution=res(.r)[1]*dist)
  .ac <- cellFromXY(.a,xy)
  .tbl <- data_frame(cell=.ac,x=xy[,1],y=xy[,2]) %>% group_by(cell)
  o <- list()
  for (i in 1:rep) {
    .s <- .tbl %>% sample_n(1) %>% ungroup()
    o[[i]] <- data.frame(.s %>% dplyr::select(x,y))
  }
  o
}

TSS <- function(eval){
  return(max(eval@TPR+eval@TNR-1))
}


TSSweight <- function(eval){
  if (TSS(eval)-0.2>0){
    (TSS(eval)-0.2)^2
  }else {
    0
  }
}

weights_tsss <- function(eval){
  if (eval-0.2>0){
    (eval-0.2)^2
  }else {
    0
  }
}

species_locations <- function(locations,species_name){
  cranch_locations <- select(locations[locations$species==species_name,],
                             longitude,
                             latitude)
  return(cranch_locations)
}

plotting <- function(rasteR,species,TSS,algo,points){
  #mods already went through mask
  mods <- rasteR
  palete = colorNumeric(c("#5E85B8","#EDF0C0","#C13127"), c(0,1), na.color = 'transparent')
  name <- species
  add_title <- algo
  m <-  leaflet()%>%
    addProviderTiles(provider = 'Esri.OceanBasemap') %>%
    addScaleBar()%>%
    addRasterImage(colors = palete, mods, opacity = .6)%>%
    addCircles(points$x,points$y)%>%
    addLegend(position = 'topright',pal = palete, values(mods), title = paste(name,' - ',add_title,', TSS: ',TSS))
  return(m)
}

main <- function(species_name,oegopsids,env,background_macaronesia_test,background_macaronesia_train.env,macaronesia.extent){
  cat(paste(species_name),'-',dim(species_locations(oegopsids,species_name))[1])
  presence_locations <- species_locations(oegopsids,species_name)
  j <- na.omit(extract(env,presence_locations))
  presence_locations <- presence_locations[-c(attr(x=j, which='na.action')),]
  presence_locations <- spThin(presence_locations,dist=500)[[1]]
  cat('-->',paste(dim(presence_locations)[1],'\n'))
  tssB <- c()
  predictionsB <- c()
  tssa <- c()
  predictionsa <- c()
  tssmx <- c()
  predictionsmX <- c()
  tssml <- c()
  predictionsml <- c()
  tssma <- c()
  predictionsma <- c()
  tssd <- c()
  predictionsd <- c()
  tssb <- c()
  predictionsb <- c()
  for (i in 1:10){
  cat(paste(i,'\n'))
  group <- kfold(presence_locations, 5)
  presence_locations.pres_train <- presence_locations[group != 1, ]
  presence_locations.pres_test <- presence_locations[group == 1, ]
  sp.env <- extract(env,presence_locations.pres_train)
  presence_absence_env <- as.data.frame(rbind(sp.env,background_macaronesia_train.env))
  presence_absence <- c(rep(1,dim(sp.env)[1]),rep(0,dim(background_macaronesia_train.env)[1]))
  ENV_PA <- presence_absence_env
  ENV_PA$presence <- presence_absence
  #BRT
  # cat('BRT\n')
  # BRT <-
  # eB <- evaluate(presence_locations.pres_test,background_macaronesia_test,GAM,env)
  # tB <- threshold(eB,stat='prevalence')
  # predictionsB <- c(predictionsB,predict(env, BRT, ext=macaronesia.extent,type="response"))
  # pB <- predictionsB>tB 
  # mods <- predictionsB
  # add_title <- 'BRT'
  # tssB <- c(tssB,eB)
  # m <- plotting(mods,species_name,tss,add_title,presence_locations)
  # saveWidget(m, file=paste(species_name,'_',add_title,'.html'))


  #cat('GAM model\n')
  GAM <- gam(presence ~ poly(BO_sstmean,3)+poly(BO_dissox,3)+poly(BO2_tempmean_bdmin,3)+poly(BO2_salinitymean_bdmin,3)+poly(BO2_salinitymean_ss,3)+poly(BO2_carbonphytomean_ss,3)+poly(BO21_carbonphytomean_bdmin,3),
              family=binomial(link='logit'),
              data=ENV_PA)
  ea <- evaluate(presence_locations.pres_test,background_macaronesia_test,GAM,env)
  #ta <- threshold(ea,stat='prevalence')
  predictionsa <- c(predictionsa,predict(env, GAM, ext=macaronesia.extent,type="response"))
  #pa <- predictionsa>ta
  tssa <- c(tssa,ea)
  

  #NN
  
  
  #cat('MaxEnt model\n')
  #maxEnt <- maxent(env,presence_locations.pres_train,)
  maxEnt <- maxent(x=presence_absence_env,presence_absence)
  #evaluation and threshold
  em <- evaluate(presence_locations.pres_test,background_macaronesia_test,maxEnt,env)
  #tm <- threshold(em,stat='prevalence')
  #prediction
  predictionsmX <- c(predictionsmX,predict(env, maxEnt, ext=macaronesia.extent,type="response"))
  #pm <- predictionsmX>tm
  tssmx <- c(tssmx,em)
  
  
  #cat('glm model\n')
  glmm <- glm(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
              family=binomial(link='logit'),
              data=ENV_PA)
  elm <- evaluate(presence_locations.pres_test,background_macaronesia_test,glmm,env)
  #tlm <- threshold(elm,stat='prevalence')
  predictionsml <- c(predictionsml,predict(env, glmm, ext=macaronesia.extent,type="response"))
  #plm <- predictionsml>tlm
  tssml <- c(tssml,elm)

  
  #cat('Mahalanobis model\n')
  maha <- mahal(sp.env)
  ema <- evaluate(presence_locations.pres_test,background_macaronesia_test,maha,env)
  #tma <-  threshold(ema,stat='prevalence')
  predictionsma <- c(predictionsma, predict(env, maha, ext=macaronesia.extent,type="response"))
  #predictionsma[predictionsma<0] <- 0
  #pma <- predictionsma>tma
  tssma <- c(tssma,ema)
  
  
  #cat('Domain model\n')
  #d <- domain(env,presence_locations.pres_train)
  d <- domain(sp.env)
  #evaluation and threshold
  ed <- evaluate(presence_locations.pres_test,background_macaronesia_test,d,env)
  #td <- threshold(ed,stat='prevalence')
  #prediction
  predictionsd <- c(predictionsd,predict(env, d, ext=macaronesia.extent,type="response"))
  #pd <- predictionsd>td
  tssd <- c(tssd,ed)
  
  
  #cat('bioclim model\n')
  # bio <- bioclim(env,presence_locations.pres_train)
  bio <- bioclim(sp.env)
  #evaluation and threshold
  eb <- evaluate(presence_locations.pres_test,background_macaronesia_test,bio,env)
  #tb <- threshold(eb,stat='prevalence')
  #prediction
  predictionsb <- c(predictionsb,predict(env, bio, ext=macaronesia.extent,type="response"))
  #pb <- predictionsb>tb
  tssb <- c(tssb,eb)
  }


  #pB <- stack()
  #pN <- stack()
  
  cat('GAM model\n')
  pa <- stack(predictionsa)
  wa <- sapply(tssa, function(x) TSSweight(x))
  tssa <- weighted.mean(sapply(tssa, function(x) TSS(x)), wa, na.rm=TRUE)
  pa <- weighted.mean(pa, wa, na.rm=TRUE)
  mods <- pa
  add_title <- 'GAM'
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))

  
  cat('GLM model\n')
  pml <- stack(predictionsml)
  wml <- sapply(tssml, function(x) TSSweight(x))
  tssml <- weighted.mean(sapply(tssml, function(x) TSS(x)), wml, na.rm=TRUE)
  pml <- weighted.mean(pml, wml, na.rm=TRUE)
  mods <- pml
  add_title <- 'GLM'
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))


  cat('Mahalanobis model\n')
  pma <- stack(predictionsma)
  wma <- sapply(tssma, function(x) TSSweight(x))
  #tssma <- sapply(tssma, function(x) TSS(x))
  tssma <- weighted.mean(sapply(tssma, function(x) TSS(x)), wma, na.rm=TRUE)
  pma <- weighted.mean(pma, wma, na.rm=TRUE)
  pma[pma<0] <- 0
  #pma[pma==NA] <- 0
  mods <- pma
  add_title <- 'Mahalanobis'
  m <-  plotting(mods,species_name,tssma,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))


  cat('Bioclim model\n')
  pb <- stack(predictionsb)
  wb <- sapply(tssb, function(x) TSSweight(x))
  #tssb <- sapply(tssb, function(x) TSS(x))
  tssb <- weighted.mean(sapply(tssb, function(x) TSS(x)), wb, na.rm=TRUE)
  pb <- weighted.mean(pb, wb, na.rm=TRUE)
  mods <- pb
  add_title <- 'Bioclim'
  m <-  plotting(mods,species_name,tssb,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))


  cat('Domain model\n')
  pd <- stack(predictionsd)
  wd <- sapply(tssd, function(x) TSSweight(x))
  #tssd <- sapply(tssd, function(x) TSS(x))
  tssd <- weighted.mean(sapply(tssd, function(x) TSS(x)), wd, na.rm=TRUE)
  pd <- weighted.mean(pd, wd, na.rm=TRUE)
  mods <- pd
  add_title <- 'Domain'
  m <-  plotting(mods,species_name,tssd,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))


  cat('MaxEnt model\n')
  pmX <- stack(predictionsmX)
  wmX <- sapply(tssmx, function(x) TSSweight(x))
  #tssmx <- sapply(tssmx, function(x) TSS(x))
  tssmx <- weighted.mean(sapply(tssmx, function(x) TSS(x)), wmX, na.rm=TRUE)
  pmX <- weighted.mean(pmX, wmX, na.rm=TRUE)
  mods <- pmX
  add_title <- 'MaxEnt'
  m <-  plotting(mods,species_name,tssmx,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))

  
  
  tsss <- c(tssa, tssml, tssma, tssb, tssd, tssmx)
  cat('Ensemble\n')
  models <- stack(pa, pml, pma, pb, pd, pmX)
  names(models) <- c("GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  #w <- sapply(list(ea, elm, ema, eb, ed, em), function(x) TSSweight(x))
  #Mean <- weighted.mean( models[[c("GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], w)
  ws <- weights_tsss(tsss)
  Mean <- weighted.mean( models[[c("GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  tsss <- weighted.mean(tsss,ws, na.rm=TRUE)
  mods <- Mean
  #mods <- mask(mods,mask_depth)
  add_title <- 'Ensemble'
  m <-  plotting(mods,species_name,tsss,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))

}




env <- load_layers(list_layers( datasets="Bio-ORACLE" )[c(21,12,158,212,310,298,335),]$layer_code)
#25,21,12,69,119,158,212,310,298,335
#bathi, SST,O2, "min deep O2",light bottom, temp bottom, sal bottom, sal surface,biomass,biomass deep
mask_depth <- load_layers(list_layers( datasets="Bio-ORACLE" )[25,]$layer_code)
mask_depth <- mask_depth[[1]]
#mask_depth <- crop(mask_depth,macaronesia.extent)
mask_depth[mask_depth>-200] <- NA
env <- mask(env,mask_depth)
macaronesia.extent <- extent(-35, -10, 12, 42)#lon,lon,lat,lat
macaronesia.env <- crop(env,macaronesia.extent)




oegopsids <- read.csv("occorrences.csv")
oegopsid_locations <- select(oegopsids,
                             longitude,
                             latitude)

squid_occorrences <-  leaflet()%>%
  addProviderTiles(provider = 'Esri.OceanBasemap') %>%
  addScaleBar()%>%
  addCircles(oegopsid_locations$longitude,oegopsid_locations$latitude)
#saveWidget(m, file=paste(name,add_title,'.html',sep = ''))
squid_occorrences


oegopsid_locations <- spThin(oegopsid_locations,dist=500)[[1]]
squid_occorrences <-  leaflet()%>%
  addProviderTiles(provider = 'Esri.OceanBasemap') %>%
  addScaleBar()%>%
  addCircles(oegopsid_locations$x,oegopsid_locations$y)
#saveWidget(m, file=paste(name,add_title,'.html',sep = ''))
squid_occorrences
#HB <- species_locations(oegopsids,'Histioteuthis bonnellii')
#squid_occorrences <-  leaflet()%>%
#  addProviderTiles(provider = 'Esri.OceanBasemap') %>%
#  addScaleBar()%>%
#  addCircles(HB$longitude,HB$latitude)
#saveWidget(m, file=paste(name,add_title,'.html',sep = ''))
#squid_occorrences



background_macaronesia <- randomPoints(env, n=1000, ext=macaronesia.extent, extf = 1.25)

colnames(background_macaronesia) = c('longitude', 'latitude')
group <- kfold(background_macaronesia, 5)
background_macaronesia_train <- background_macaronesia[group != 1, ]
background_macaronesia_train.env <- extract(env, background_macaronesia_train)
background_macaronesia_test <- background_macaronesia[group == 1, ]


for (species_name in c('Chiroteuthis veranii',#)){#unique(oegopsids$species)
                       'Cranchia scabra',
                       'Histioteuthis bonnellii',
                       'Histioteuthis reversa',
                       'Illex coindetii',
                       'Illex illecebrosus',
                       'Liocranchia reinhardtii',
                       'Lycoteuthis lorigera',
                       'Nototodarus gouldi',
                       'Nototodarus sloanii',
                       'Ommastrephes bartramii',
                       'Onychoteuthis banksii',
                       'Onykia robsoni',
                       'Pterygioteuthis gemmata',
                       'Pterygioteuthis giardi',#not working, perhaps no need to remove NAs?
                       'Pyroteuthis margaritifera',
                       'Sthenoteuthis oualaniensis',
                       'Teuthowenia megalops',
                       'Teuthowenia pellucida',
                       'Todarodes filippovae',
                       'Todarodes sagittatus',
                       'Todaropsis eblanae')){#unique(oegopsids$species)
  if (dim(species_locations(oegopsids,species_name))[1]>100){
    main(species_name,oegopsids,env,background_macaronesia_test,background_macaronesia_train.env,macaronesia.extent)
    }
}






main('Pterygioteuthis giardi',oegopsids,env,background_macaronesia_test,background_macaronesia_train.env,macaronesia.extent)


species_name <- "Cranchia scabra"
cranch_locations <- species_locations(oegopsids,species_name)
cat(dim(cranch_locations)[1],'\n')
cranch_locations <- spThin(cranch_locations,dist=500)[[1]]
group <- kfold(cranch_locations, 5)
cranch_locations.pres_train <- cranch_locations[group != 1, ]
cranch_locations.pres_test <- cranch_locations[group == 1, ]

cranch_env <- extract(env,cranch_locations.pres_train)
cranch_env <- na.omit(cranch_env)



#models

d <- domain(cranch_env)
bio <- bioclim(cranch_env)
maxEnt <- maxent(env,cranch_locations.pres_train)
#evaluation and threshold
ed <- evaluate(cranch_locations.pres_test,background_macaronesia_test,d,env)
td <- threshold(ed,stat='prevalence')
eb <- evaluate(cranch_locations.pres_test,background_macaronesia_test,bio,env)
tb <- threshold(eb,stat='prevalence')
em <- evaluate(cranch_locations.pres_test,background_macaronesia_test,maxEnt,env)
tm <- threshold(em,stat='prevalence')
#prediction

predictionsb <- predict(env, bio, ext=macaronesia.extent,type="response")
pb <- predictionsb>tb
predictionsd <- predict(env, d, ext=macaronesia.extent,type="response")
pd <- predictionsd>td
predictionsmX <- predict(env, maxEnt, ext=macaronesia.extent,type="response")
pm <- predictionsmX>tm
models <- stack(predictionsb,predictionsd,predictionsmX)
names(models) <- c("domain", "bioclim", "maxent")
auc <- sapply(list(ed, eb, em), function(x) x@auc)
w <- (auc-0.5)^2
Mean <- weighted.mean( models[[c("domain", "bioclim", "maxent")]], w)




#plotting
mods <- predictionsmX
palete = colorNumeric(c("#5E85B8","#EDF0C0","#C13127"), values(mods), na.color = 'transparent')
add_title <- 'm'
m <-  leaflet()%>%
	addProviderTiles(provider = 'Esri.OceanBasemap') %>%
	addScaleBar()%>%
	addRasterImage(colors = palete, mods, opacity = .6)%>%
	addLegend(position = 'topright',pal = palete, values(mods), title = paste(species_name,add_title,sep=' - '))
saveWidget(m, file=paste(species_name,'_',add_title,'.html'))

m



#oegopsids <- read.csv("occorrences.csv")
#ceph1 <- enmtools.species()
#ceph1$species.name <- "Cranchia scabra"
#ceph1$presence.points <- oegopsids[oegopsids$species==ceph1$species.name,]
#ceph1 <- check.species(ceph1)
#cepth1.mx <- enmtools.maxent(ceph1, env, test.prop = 0.3)
#cepth1.mx$test.evaluation
#cepth1.mx$training.evaluation
#response.plots
#ceph1.bc <- enmtools.bc(ceph1, env, test.prop = 0.3)
#macaronesia.batho <- crop(env[[1]],macaronesia.extent)
#mask_depth <- macaronesia.batho
#mask_depth[mask_depth>-200] <- NA
#macaronesia.batho <- mask(macaronesia.batho,mask_depth)
for (species_name in unique(oegopsids$species)){
  if (dim(species_locations(oegopsids,species_name))[1]>100){
    cat(paste(species_name),' - ',dim(species_locations(oegopsids,species_name))[1],'\n')
  }}
#####



#####
cat(paste(species_name),' - ',dim(species_locations(oegopsids,species_name))[1],'\n')
presence_locations <- species_locations(oegopsids,species_name)
group <- kfold(presence_locations, 5)
presence_locations.pres_train <- presence_locations[group != 1, ]
presence_locations.pres_test <- presence_locations[group == 1, ]
sp.env <- extract(env,presence_locations.pres_train)
sp.env <- na.omit(sp.env)

#models
d <- domain(env,presence_locations.pres_train)
bio <- bioclim(env,presence_locations.pres_train)
maxEnt <- maxent(env,presence_locations.pres_train)

#evaluation and threshold
ed <- evaluate(presence_locations.pres_test,background_macaronesia_test,d,env)
td <- threshold(ed,stat='prevalence')
eb <- evaluate(presence_locations.pres_test,background_macaronesia_test,bio,env)
tb <- threshold(eb,stat='prevalence')
em <- evaluate(presence_locations.pres_test,background_macaronesia_test,maxEnt,env)
tm <- threshold(em,stat='prevalence')

#prediction rasters
predictionsb <- predict(env, bio, ext=macaronesia.extent,type="response")
pb <- predictionsb>tb
predictionsd <- predict(env, d, ext=macaronesia.extent,type="response")
pd <- predictionsd>td
predictionsmX <- predict(env, maxEnt, ext=macaronesia.extent,type="response")
pm <- predictionsmX>tm
models <- stack(predictionsb,predictionsd,predictionsmX)
names(models) <- c("domain", "bioclim", "maxent")
auc <- sapply(list(ed, eb, em), function(x) x@auc)
w <- (auc-0.5)^2
Mean <- weighted.mean( models[[c("domain", "bioclim", "maxent")]], w)

#plotting
#ensemble
mods <- Mean
#mods <- mask(mods,mask_depth)
add_title <- 'Ensemble'
m <-  plotting(mods,species_name,'?',add_title)
saveWidget(m, file=paste(species_name,'_',add_title,'.html'))
#maxent
mods <- predictionsmX
#mods <- mask(mods,mask_depth)
add_title <- 'MaxEnt'
AUC <- em@auc
m <-  plotting(mods,species_name,AUC,add_title)
saveWidget(m, file=paste(species_name,'_',add_title,'.html'))
#bioclim
mods <- predictionsb
#mods <- mask(mods,mask_depth)
add_title <- 'Bioclim'
AUC <- eb@auc
m <-  plotting(mods,species_name,AUC,add_title)
saveWidget(m, file=paste(species_name,'_',add_title,'.html'))
#Domain
mods <- predictionsd
#mods <- mask(mods,mask_depth)
add_title <- 'Domain'
AUC <- ed@auc
m <-  plotting(mods,species_name,AUC,add_title)
saveWidget(m, file=paste(species_name,'_',add_title,'.html'))
