library(dismo)
library(ncdf4)
library(dplyr)
library(raster)
library(leaflet)
library(htmlwidgets)
library(sdmpredictors)
#library(ENMTools)
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
  if (TSS(eval)-0.7>0){
    (TSS(eval)-0.5)^2
  }else {
    0
  }
}

weights_tsss <- function(eval){
  if (eval-0.7>0){
    (eval-0.5)^2
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
    #addCircles(points$x,points$y)%>%
    addLegend(position = 'topright',pal = palete, values(mods), title = paste(name,' - ',add_title,', TSS: ',TSS))
  return(m)
}

plottingCompare <- function(rasteR,rasterG,species,TSS,algo){
  #mods already went through mask
  paleteR = colorNumeric(c("#000000", "#0000FF", "#FF0000"), c(-1,1), na.color = 'transparent')
  paleteG = colorNumeric(c("#000000", "#0000FF", "#00FF00"), c(-1,1), na.color = 'transparent')
  name <- species
  add_title <- algo
  m <-  leaflet()%>%
    addProviderTiles(provider = 'Esri.OceanBasemap') %>%
    addScaleBar()%>%
    addRasterImage(colors = paleteR, rasterR, opacity = .5)%>%
    addRasterImage(colors = paleteG, rasterG, opacity = .5)%>%
    # addLegend(position = 'topright',pal = paleteR, values(c(-1,1)), title = paste(name,' - 2100 - ',add_title,', TSS: ',TSS))%>%
    # addLegend(position = 'bottomright',pal = paleteG, values(c(-1,1)), title = paste(name,' - 2050'))
  return(m)
}

main180 <- function(species_name,
                    oegopsids,
                    env,
                    future_env1,
                    future_env2,
                    future_env3,
                    future_env4,
                    future_env5,
                    future_env6,
                    future_env7,
                    future_env8,
                    background_macaronesia_test,
                    background_macaronesia_train.env,
                    macaronesia.extent,
                    macaronesia.extent2){
  
  cat(paste(species_name),'-',dim(species_locations(oegopsids,species_name))[1])
  presence_locations <- species_locations(oegopsids,species_name)
  j <- na.omit(extract(env,presence_locations))
  presence_locations <- presence_locations[-c(attr(x=j, which='na.action')),]
  presence_locations <- spThin(presence_locations,dist=20)[[1]]
  cat('-->',paste(dim(presence_locations)[1],'\n'))
  
  envs <- list(
    env,
    future_env1,
    future_env2,
    future_env3,
    future_env4,
    future_env5,
    future_env6,
    future_env7,
    future_env8
    )
  env_names = list(
    '_Present_',
    '_2050.RCP26_',
    '_2050.RCP45_',
    '_2050.RCP60_',
    '_2050.RCP85_',
    '_2100.RCP26_',
    '_2100.RCP45_',
    '_2100.RCP60_',
    '_2100.RCP85_'
    )
  predictionsmX = list(
      predictions = list(c(),c()),
      predictions1 = list(c(),c()),
      predictions2 = list(c(),c()),
      predictions3 = list(c(),c()),
      predictions4 = list(c(),c()),
      predictions5 = list(c(),c()),
      predictions6 = list(c(),c()),
      predictions7 = list(c(),c()),
      predictions8 = list(c(),c()))
  tssmX <- c()#True Skill score

  threshesmX <- c()#thresholds
  for (i in 1:1){
    cat(paste(i,'\n'))
    group <- kfold(presence_locations, 4)
    presence_locations.pres_train <- presence_locations[group != 1, ]
    presence_locations.pres_test <- presence_locations[group == 1, ]
    sp.env <- extract(env,presence_locations.pres_train)
    presence_absence_env <- as.data.frame(rbind(sp.env,background_macaronesia_train.env))
    presence_absence <- c(rep(1,dim(sp.env)[1]),rep(0,dim(background_macaronesia_train.env)[1]))
    ENV_PA <- presence_absence_env
    ENV_PA$presence <- presence_absence
    

    maxEnt <- maxent(x=presence_absence_env,presence_absence,silent=TRUE)
    em <- evaluate(presence_locations.pres_test,background_macaronesia_test,maxEnt,env)
    tssmX <- c(tssmX,em)
    w <- TSSweight(em)
    if (i==1){
      for (ii in 1:length(envs)){
        predictionsmX[[ii]][[1]] <- predict(envs[[ii]], maxEnt, ext=macaronesia.extent,type="response")*w
        predictionsmX[[ii]][[2]] <- predict(envs[[ii]], maxEnt, ext=macaronesia.extent2,type="response")*w
      }
      }else {
        for (ii in 1:length(envs)){
          predictionsmX[[ii]][[1]] <- predictionsmX[[ii]][[1]] + predict(envs[[ii]], maxEnt, ext=macaronesia.extent,type="response")*w
          predictionsmX[[ii]][[2]] <- predictionsmX[[ii]][[2]] + predict(envs[[ii]], maxEnt, ext=macaronesia.extent2,type="response")*w
        }
      }
  }
  
  cat('MaxEnt model\n')
  w <- sum(tssmX)
  wmX <- sapply(tssmX, function(x) TSSweight(x))
  #tssmx <- sapply(tssmx, function(x) TSS(x))
  tssmX <- weighted.mean(sapply(tss$tssmX, function(x) TSS(x)), wmX, na.rm=TRUE)
  add_title <- 'MaxEnt'
  for (i in 1:length(env_names)){
    predictionsmX[[i]][[1]] <- predictionsmX[[i]][[1]]/w
    predictionsmX[[i]][[2]] <- predictionsmX[[i]][[2]]/w
    predictionsmX[[i]] <- do.call(merge,list(predictionsmX[[i]][[1]],predictionsmX[[i]][[2]]))
    mods <- predictionsmX[[i]]
    m <-  plotting(mods,species_name,tssmX,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
  }

  n <- c('_RCP26_','_RCP45_','_RCP60_','_RCP85_')
  cat('Different RCP scenarios\n')
  for (i in 1:4){
    rasterR <- predictionsmX[[i+5]] - predictionsmX[[1]]
    rasterG <- predictionsmX[[i+1]] - predictionsmX[[1]]
    m <- plottingCompare(rasteR,rasterG,species_name,tsss,add_title)
    saveWidget(m, file=paste(species_name,n[i],add_title,'.html', sep=''))
  }
  
}





main <- function(species_name,
                 oegopsids,
                 env,
                 future_env1,
                 future_env2,
                 future_env3,
                 future_env4,
                 future_env5,
                 future_env6,
                 future_env7,
                 future_env8,
                 background_macaronesia_test,
                 background_macaronesia_train.env,
                 macaronesia.extent){
  cat(paste(species_name),'-',dim(species_locations(oegopsids,species_name))[1])
  presence_locations <- species_locations(oegopsids,species_name)
  j <- na.omit(extract(env,presence_locations))
  presence_locations <- presence_locations[-c(attr(x=j, which='na.action')),]
  presence_locations <- spThin(presence_locations,dist=20)[[1]]
  cat('-->',paste(dim(presence_locations)[1],'\n'))
  
  envs <- list(
    env,
    future_env1,
    future_env2,
    future_env3,
    future_env4,
    future_env5,
    future_env6,
    future_env7,
    future_env8
    )
  env_names = list(
    '_Present_',
    '_2050.RCP26_',
    '_2050.RCP45_',
    '_2050.RCP60_',
    '_2050.RCP85_',
    '_2100.RCP26_',
    '_2100.RCP45_',
    '_2100.RCP60_',
    '_2100.RCP85_'
    )
  predictionsmX = list(
      predictions = c(),
      predictions1 = c(),
      predictions2 = c(),
      predictions3 = c(),
      predictions4 = c(),
      predictions5 = c(),
      predictions6 = c(),
      predictions7 = c(),
      predictions8 = c())
  tssmX <- c()#True Skill score

  threshesmX <- c()#thresholds
  for (i in 1:1){
    cat(paste(i,'\n'))
    group <- kfold(presence_locations, 4)
    presence_locations.pres_train <- presence_locations[group != 1, ]
    presence_locations.pres_test <- presence_locations[group == 1, ]
    sp.env <- extract(env,presence_locations.pres_train)
    presence_absence_env <- as.data.frame(rbind(sp.env,background_macaronesia_train.env))
    presence_absence <- c(rep(1,dim(sp.env)[1]),rep(0,dim(background_macaronesia_train.env)[1]))
    ENV_PA <- presence_absence_env
    ENV_PA$presence <- presence_absence
    

    maxEnt <- maxent(x=presence_absence_env,presence_absence,silent=TRUE)
    em <- evaluate(presence_locations.pres_test,background_macaronesia_test,maxEnt,env)
    tssmX <- c(tssmX,em)
    w <- TSSweight(em)
    if (i==1){
      for (ii in 1:length(envs)){
        predictionsmX[[ii]] <- predict(envs[[ii]], maxEnt, ext=macaronesia.extent,type="response")*w
      }
      }else {
        for (ii in 1:length(envs)){
          predictionsmX[[ii]] <- predictionsmX[[ii]] + predict(envs[[ii]], maxEnt, ext=macaronesia.extent,type="response")*w
        }
      }
  }
  
  cat('MaxEnt model\n')
  w <- sum(tssmX)
  wmX <- sapply(tssmX, function(x) TSSweight(x))
  #tssmx <- sapply(tssmx, function(x) TSS(x))
  tssmX <- weighted.mean(sapply(tss$tssmX, function(x) TSS(x)), wmX, na.rm=TRUE)
  add_title <- 'MaxEnt'
  for (i in 1:length(env_names)){
    predictionsmX[[i]] <- predictionsmX[[i]]/w
    mods <- predictionsmX[[i]]
    m <-  plotting(mods,species_name,tssmX,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
  }

  n <- c('_RCP26_','_RCP45_','_RCP60_','_RCP85_')
  cat('Different RCP scenarios\n')
  for (i in 1:4){
    rasterR <- predictionsmX[[i+5]] - predictionsmX[[1]]
    rasterG <- predictionsmX[[i+1]] - predictionsmX[[1]]
    m <- plottingCompare(rasteR,rasterG,species_name,tsss,add_title)
    saveWidget(m, file=paste(species_name,n[i],add_title,'.html', sep=''))
  }
}


#lst <- list.files(pattern='asc$')
lste <- list.files(pattern='^Present')
env <- stack(lste)
mask_depth <- load_layers(list_layers( datasets="Bio-ORACLE" )[25,]$layer_code)
mask_depth <- mask_depth[[1]]
mask_depth[mask_depth>-20] <- NA
mask_depth[mask_depth<-1200] <- NA
crs(env) <- crs(mask_depth)
env <- mask(env,mask_depth)

lst1 <- list.files(pattern='^2050AOGCM.RCP26.')
future_env1 <- stack(lst1)
crs(future_env1) <- crs(mask_depth)
future_env1 <- mask(future_env1,mask_depth)
names(future_env1) <- names(env)


lst2 <- list.files(pattern='^2050AOGCM.RCP45.')
future_env2 <- stack(lst2)
crs(future_env2) <- crs(mask_depth)
future_env2 <- mask(future_env2,mask_depth)
names(future_env2) <- names(env)

lst3 <- list.files(pattern='^2050AOGCM.RCP60.')
future_env3 <- stack(lst3)
crs(future_env3) <- crs(mask_depth)
future_env3 <- mask(future_env3,mask_depth)
names(future_env3) <- names(env)

lst4 <- list.files(pattern='^2050AOGCM.RCP85.')
future_env4 <- stack(lst4)
crs(future_env4) <- crs(mask_depth)
future_env4 <- mask(future_env4,mask_depth)
names(future_env4) <- names(env)

lst5 <- list.files(pattern='^2100AOGCM.RCP26.')
future_env5 <- stack(lst5)
crs(future_env5) <- crs(mask_depth)
future_env5 <- mask(future_env5,mask_depth)
names(future_env5) <- names(env)

lst6 <- list.files(pattern='^2100AOGCM.RCP45.')
future_env6 <- stack(lst6)
crs(future_env6) <- crs(mask_depth)
future_env6 <- mask(future_env6,mask_depth)
names(future_env6) <- names(env)


lst7 <- list.files(pattern='^2100AOGCM.RCP60.')
future_env7 <- stack(lst7)
crs(future_env7) <- crs(mask_depth)
future_env7 <- mask(future_env7,mask_depth)
names(future_env7) <- names(env)

lst8 <- list.files(pattern='^2100AOGCM.RCP85.')
future_env8 <- stack(lst8)
crs(future_env8) <- crs(mask_depth)
future_env8 <- mask(future_env8,mask_depth)
names(future_env8) <- names(env)









oegopsids <- read.csv("occorrences.csv")
# oegopsid_locations <- select(oegopsids,
#                              longitude,
#                              latitude)

background <- randomPoints(env, n=1000, extf = 1.25)
colnames(background) = c('longitude', 'latitude')
group <- kfold(background, 4)
background_train <- background[group != 1, ]
background_train.env <- extract(env, background_train)
background_test <- background[group == 1, ]
#total.extent <- extent(-35, -10, 12, 42)#lon,lon,lat,lat



nsloanii <- select(oegopsids[oegopsids$species=='Nototodarus sloanii',],
                   species,
                   longitude,
                   latitude)


# main180('Nototodarus sloanii',
#      nsloanii,env,
#      future_env1,
#      future_env2,
#      future_env3,
#      future_env4,
#      future_env5,
#      future_env6,
#      future_env7,
#      future_env8,
#      background_test,background_train.env,
#      extent(156, 180, -56, -20),extent(-180+0.008333333, -167, -56, -20))

species_name <- 'Nototodarus sloanii'
oegopsids <- nsloanii
background_macaronesia_test <-background_test
background_macaronesia_train.env<-background_train.env
macaronesia.extent<-extent(156, 180, -56, -20)
macaronesia.extent2<-extent(-180+0.008333333, -167, -56, -20)

oegopsids <- read.csv("occorrences.csv")
illex <- oegopsids[oegopsids$latitude<(12.0*-1),]
illex <- select(illex[illex$species=='Illex argentinus',],
                species,
                longitude,
                latitude)
# main('Illex argentinus',
#      illex,env,
#      future_env1,
#      future_env2,
#      future_env3,
#      future_env4,
#      future_env5,
#      future_env6,
#      future_env7,
#      future_env8,background_test,background_train.env,extent(-70, -23, -60, -10))

species_name <- 'Illex argentinus'#'Nototodarus sloanii'
oegopsids <- illex#nsloanii
macaronesia.extent<-extent(-70, -23, -55, -10)#extent(156, 180, -56, -20)


oegopsids <- read.csv("occorrences.csv")
tpacificus<- select(oegopsids[oegopsids$species=='Todarodes pacificus',],
                    species,
                    longitude,
                    latitude)
# main180('Todarodes pacificus',tpacificus,env,
#         future_env1,
#         future_env2,
#         future_env3,
#         future_env4,
#         future_env5,
#         future_env6,
#         future_env7,
#         future_env8,background_test,
#         background_train.env,extent(105, 180, 18, 65),extent(-180+0.008333333, -130, 42, 65))

species_name <- 'Todarodes pacificus'#'Illex argentinus'#'Nototodarus sloanii'
oegopsids <- tpacificus#illex#nsloanii
macaronesia.extent<-extent(105, 180, 18, 65)
macaronesia.extent2<-extent(-180+0.008333333, -130, 42, 65)



oegopsids <- read.csv("occorrences.csv")
bmagister<- select(oegopsids[oegopsids$species=='Berryteuthis magister',],
                   species,
                   longitude,
                   latitude)
# main180('Berryteuthis magister',bmagister,env,
#         future_env1,
#         future_env2,
#         future_env3,
#         future_env4,
#         future_env5,
#         future_env6,
#         future_env7,
#         future_env8,
#         background_test,background_train.env,
#         extent(126, 180, 35, 65),extent(-180+0.008333333, -120, 42, 65))

species_name <- 'Berryteuthis magister'#'Illex argentinus'#'Nototodarus sloanii'
oegopsids <- bmagister#illex#nsloanii
macaronesia.extent<-extent(126, 180, 35, 65)#extent(-100, -23, -80, 10)#extent(156, 180, -56, -20)
macaronesia.extent2<-extent(-180+0.008333333, -120, 42, 65)#extent(-180+0.008333333, -167, -56, -20)
