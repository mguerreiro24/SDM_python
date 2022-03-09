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
  if (TSS(eval)-0.5>0){
    (TSS(eval)-0.5)^2
  }else {
    0
  }
}

weights_tsss <- function(eval){
  if (eval-0.5>0){
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
    addCircles(points$x,points$y)%>%
    addLegend(position = 'topright',pal = palete, values(mods), title = paste(name,' - ',add_title,', TSS: ',TSS))
  return(m)
}

main180 <- function(species_name,
                    oegopsids,
                    env,
                    future_env,
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
  tssnn <- c()
  predictionsnn <- c()
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
  
  predictionsnn2 <- c()
  predictionsB2 <- c()
  predictionsa2 <- c()
  predictionsmX2 <- c()
  predictionsml2 <- c()
  predictionsma2 <- c()
  predictionsd2 <- c()
  predictionsb2 <- c()
  
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
    
    
    BRT <- gbm(presence ~ Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + Present.Benthic.Mean.Depth.Salinity.Lt.max + Present.Benthic.Mean.Depth.Salinity.Lt.min + Present.Benthic.Mean.Depth.Salinity.Mean + Present.Benthic.Mean.Depth.Temperature.Lt.max + Present.Benthic.Mean.Depth.Temperature.Lt.min + Present.Benthic.Mean.Depth.Temperature.Mean + Present.Surface.Current.Velocity.Mean.asc.BOv2_1 + Present.Surface.Salinity.Lt.max + Present.Surface.Salinity.Lt.min + Present.Surface.Salinity.Mean + Present.Surface.Temperature.Lt.max + Present.Surface.Temperature.Lt.min + Present.Surface.Temperature.Mean,
               cv.folds = 100,
               interaction.depth = 3,
               distribution = "bernoulli",
               shrinkage = 0.01,
               verbose = FALSE,
               n.trees = 500,
               data=ENV_PA)
    eB <- evaluate(presence_locations.pres_test,background_macaronesia_test,BRT,env)
    tB <- threshold(eB,stat='prevalence')
    pB <- predict(future_env, BRT, ext=macaronesia.extent,type="response")>tB 
    predictionsB <- c(predictionsB,pB)
    pB <- predict(future_env, BRT, ext=macaronesia.extent2,type="response")>tB 
    predictionsB2 <- c(predictionsB2,pB)
    tssB <- c(tssB,eB)
    
    
    GAM <- gam(presence ~ poly(Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1,3) + poly(Present.Benthic.Mean.Depth.Salinity.Lt.max,3) + poly(Present.Benthic.Mean.Depth.Salinity.Lt.min,3) + poly(Present.Benthic.Mean.Depth.Salinity.Mean,3) + poly(Present.Benthic.Mean.Depth.Temperature.Lt.max,3) + poly(Present.Benthic.Mean.Depth.Temperature.Lt.min,3) + poly(Present.Benthic.Mean.Depth.Temperature.Mean,3) + poly(Present.Surface.Current.Velocity.Mean.asc.BOv2_1,3) + poly(Present.Surface.Salinity.Lt.max,3) + poly(Present.Surface.Salinity.Lt.min,3) + poly(Present.Surface.Salinity.Mean,3) + poly(Present.Surface.Temperature.Lt.max,3) + poly(Present.Surface.Temperature.Lt.min,3) + poly(Present.Surface.Temperature.Mean,3),
               family=binomial(link='logit'),
               data=ENV_PA)
    ea <- evaluate(presence_locations.pres_test,background_macaronesia_test,GAM,env)
    ta <- threshold(ea,stat='prevalence')
    pa <- predict(future_env, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa <- c(predictionsa,pa)
    pa <- predict(future_env, GAM, ext=macaronesia.extent2,type="response")>ta
    predictionsa2 <- c(predictionsa2,pa)
    tssa <- c(tssa,ea)
    
    
    # #cat('neural')
    # NN <- nnet(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
    #            data=ENV_PA,
    #            size=10)
    # enn <- evaluate(presence_locations.pres_test,background_macaronesia_test,NN,env)
    # predictionsnn <- c(predictionsnn,predict(env, NN, ext=macaronesia.extent,type="response"))
    # tssnn <- c(tssnn,enn)
    
    
    #maxEnt <- maxent(env,presence_locations.pres_train,)
    maxEnt <- maxent(x=presence_absence_env,presence_absence,silent=TRUE)
    em <- evaluate(presence_locations.pres_test,background_macaronesia_test,maxEnt,env)
    tm <- threshold(em,stat='prevalence')
    pm <- predict(future_env, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX <- c(predictionsmX,pm)
    pm <- predict(future_env, maxEnt, ext=macaronesia.extent2,type="response")>tm
    predictionsmX2 <- c(predictionsmX2,pm)
    tssmx <- c(tssmx,em)
    
    
    glmm <- glm(presence ~ Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + Present.Benthic.Mean.Depth.Salinity.Lt.max + Present.Benthic.Mean.Depth.Salinity.Lt.min + Present.Benthic.Mean.Depth.Salinity.Mean + Present.Benthic.Mean.Depth.Temperature.Lt.max + Present.Benthic.Mean.Depth.Temperature.Lt.min + Present.Benthic.Mean.Depth.Temperature.Mean + Present.Surface.Current.Velocity.Mean.asc.BOv2_1 + Present.Surface.Salinity.Lt.max + Present.Surface.Salinity.Lt.min + Present.Surface.Salinity.Mean + Present.Surface.Temperature.Lt.max + Present.Surface.Temperature.Lt.min + Present.Surface.Temperature.Mean,
                family=binomial(link='logit'),
                data=ENV_PA)
    elm <- evaluate(presence_locations.pres_test,background_macaronesia_test,glmm,env)
    tlm <- threshold(elm,stat='prevalence')
    plm <- predict(future_env, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml <- c(predictionsml,plm)
    plm <- predict(future_env, glmm, ext=macaronesia.extent2,type="response")>tlm
    predictionsml2 <- c(predictionsml2,plm)
    tssml <- c(tssml,elm)
    
    
    maha <- mahal(sp.env)
    ema <- evaluate(presence_locations.pres_test,background_macaronesia_test,maha,env)
    tma <-  threshold(ema,stat='prevalence')
    pma <- predict(future_env, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma <- c(predictionsma, pma)
    pma <- predict(future_env, maha, ext=macaronesia.extent2,type="response")>tma
    predictionsma2 <- c(predictionsma2, pma)
    tssma <- c(tssma,ema)
    
    
    #d <- domain(env,presence_locations.pres_train)
    d <- domain(sp.env)
    ed <- evaluate(presence_locations.pres_test,background_macaronesia_test,d,env)
    td <- threshold(ed,stat='prevalence')
    pd <- predict(future_env, d, ext=macaronesia.extent,type="response")>td
    predictionsd <- c(predictionsd,pd)
    pd <- predict(future_env, d, ext=macaronesia.extent2,type="response")>td
    predictionsd2 <- c(predictionsd2,pd)
    tssd <- c(tssd,ed)
    
    
    # bio <- bioclim(env,presence_locations.pres_train)
    bio <- bioclim(sp.env)
    eb <- evaluate(presence_locations.pres_test,background_macaronesia_test,bio,env)
    tb <- threshold(eb,stat='prevalence')
    pb <- predict(future_env, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb <- c(predictionsb,pb)
    pb <- predict(future_env, bio, ext=macaronesia.extent2,type="response")>tb
    predictionsb2 <- c(predictionsb2,pb)
    tssb <- c(tssb,eb)
  }
  
  cat('BRT model\n')
  pB <- stack(predictionsB)
  wB <- sapply(tssB, function(x) TSSweight(x))
  tssB <- weighted.mean(sapply(tssB, function(x) TSS(x)), wB, na.rm=TRUE)
  pB <- weighted.mean(pB, wB, na.rm=TRUE)
  pB2 <- stack(predictionsB2)
  pB2 <- weighted.mean(pB2, wB, na.rm=TRUE)
  pB <- do.call(merge,list(pB,pB2))
  mods <- pB
  add_title <- 'BRT'
  m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))
  
  
  # cat("Neural Network model\n")
  # pN <- stack(predictionsnn)
  # wnn <- sapply(tssnn, function(x) TSSweight(x))
  # tssnn <- weighted.mean(sapply(tssnn, function(x) TSS(x)), wnn, na.rm=TRUE)
  # pN <- weighted.mean(pN, wnn, na.rm=TRUE)
  # mods <- pN
  # add_title <- 'Neural Network'
  # m <-  plotting(mods,species_name,tssnn,add_title,presence_locations)
  # saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))
  
  
  cat('GAM model\n')
  pa <- stack(predictionsa)
  wa <- sapply(tssa, function(x) TSSweight(x))
  tssa <- weighted.mean(sapply(tssa, function(x) TSS(x)), wa, na.rm=TRUE)
  pa <- weighted.mean(pa, wa, na.rm=TRUE)
  pa2 <- stack(predictionsa2)
  pa2<- weighted.mean(pa2, wa, na.rm=TRUE)
  pa <- do.call(merge,list(pa, pa2))
  mods <- pa
  add_title <- 'GAM'
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))
  
  
  cat('GLM model\n')
  pml <- stack(predictionsml)
  wml <- sapply(tssml, function(x) TSSweight(x))
  tssml <- weighted.mean(sapply(tssml, function(x) TSS(x)), wml, na.rm=TRUE)
  pml <- weighted.mean(pml, wml, na.rm=TRUE)
  pml2 <- stack(predictionsml2)
  pml2 <- weighted.mean(pml2, wml, na.rm=TRUE)
  pml <- do.call(merge,list(pml, pml2))
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
  pma2 <- stack(predictionsma2)
  pma2 <- weighted.mean(pma2, wma, na.rm=TRUE)
  pma <- do.call(merge,list(pma, pma2))
  #pma[pma==NA] <- 0
  mods <- pma
  add_title <- 'Mahalanobis'
  m <-  plotting(mods,species_name,tssma,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))
  
  
  cat('Bioclim model\n')#filter of .2 TSS is removing the layers
  pb <- stack(predictionsb)
  wb <- sapply(tssb, function(x) TSSweight(x))
  #tssb <- sapply(tssb, function(x) TSS(x))
  tssb <- weighted.mean(sapply(tssb, function(x) TSS(x)), wb, na.rm=TRUE)
  pb <- weighted.mean(pb, wb, na.rm=TRUE)
  pb2 <- stack(predictionsb2)
  pb2 <- weighted.mean(pb2, wb, na.rm=TRUE)
  pb <- do.call(merge,list(pb, pb2))
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
  pd2 <- stack(predictionsd2)
  pd2 <- weighted.mean(pd2, wd, na.rm=TRUE)
  pd <- do.call(merge,list(pd, pd2))
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
  pmX2 <- stack(predictionsmX2)
  pmX2 <- weighted.mean(pmX2, wmX, na.rm=TRUE)
  pmX <- do.call(merge,list(pmX, pmX2))
  mods <- pmX
  add_title <- 'MaxEnt'
  m <-  plotting(mods,species_name,tssmx,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))
  
  
  
  tsss <- c(tssB, tssa, tssml, tssma, tssb, tssd, tssmx)#tssnn
  cat('Ensemble\n')
  models <- stack(pB, pa, pml, pma, pb, pd, pmX)#pN
  names(models) <- c("BRT","GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  #w <- sapply(list(ea, elm, ema, eb, ed, em), function(x) TSSweight(x))
  #Mean <- weighted.mean( models[[c("GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], w)
  ws <- weights_tsss(tsss)
  Mean <- weighted.mean( models[[c("BRT", "GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  tsss <- weighted.mean(tsss,ws, na.rm=TRUE)
  mods <- Mean
  #mods <- mask(mods,mask_depth)
  add_title <- 'Ensemble'
  m <-  plotting(mods,species_name,tsss,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))
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
  tssnn <- c()
  predictionsnn <- c()
  tssB <- c()
  predictionsB <- c()
  predictionsB1 <- c()
  predictionsB2 <- c()
  predictionsB3 <- c()
  predictionsB4 <- c()
  predictionsB5 <- c()
  predictionsB6 <- c()
  predictionsB7 <- c()
  predictionsB8 <- c()
  tssa <- c()
  predictionsa <- c()
  predictionsa1 <- c()
  predictionsa2 <- c()
  predictionsa3 <- c()
  predictionsa4 <- c()
  predictionsa5 <- c()
  predictionsa6 <- c()
  predictionsa7 <- c()
  predictionsa8 <- c()
  tssmx <- c()
  predictionsmX <- c()
  predictionsmX1 <- c()
  predictionsmX2 <- c()
  predictionsmX3 <- c()
  predictionsmX4 <- c()
  predictionsmX5 <- c()
  predictionsmX6 <- c()
  predictionsmX7 <- c()
  predictionsmX8 <- c()
  tssml <- c()
  predictionsml <- c()
  predictionsml1 <- c()
  predictionsml2 <- c()
  predictionsml3 <- c()
  predictionsml4 <- c()
  predictionsml5 <- c()
  predictionsml6 <- c()
  predictionsml7 <- c()
  predictionsml8 <- c()
  tssma <- c()
  predictionsma <- c()
  predictionsma1 <- c()
  predictionsma2 <- c()
  predictionsma3 <- c()
  predictionsma4 <- c()
  predictionsma5 <- c()
  predictionsma6 <- c()
  predictionsma7 <- c()
  predictionsma8 <- c()
  tssd <- c()
  predictionsd <- c()
  predictionsd1 <- c()
  predictionsd2 <- c()
  predictionsd3 <- c()
  predictionsd4 <- c()
  predictionsd5 <- c()
  predictionsd6 <- c()
  predictionsd7 <- c()
  predictionsd8 <- c()
  tssb <- c()
  predictionsb <- c()
  predictionsb1 <- c()
  predictionsb2 <- c()
  predictionsb3 <- c()
  predictionsb4 <- c()
  predictionsb5 <- c()
  predictionsb6 <- c()
  predictionsb7 <- c()
  predictionsb8 <- c()
  
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
    
    
    #    BRT <- gbm(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
    BRT <- gbm(presence ~ Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + Present.Benthic.Mean.Depth.Salinity.Lt.max + Present.Benthic.Mean.Depth.Salinity.Lt.min + Present.Benthic.Mean.Depth.Salinity.Mean + Present.Benthic.Mean.Depth.Temperature.Lt.max + Present.Benthic.Mean.Depth.Temperature.Lt.min + Present.Benthic.Mean.Depth.Temperature.Mean + Present.Surface.Current.Velocity.Mean.asc.BOv2_1 + Present.Surface.Salinity.Lt.max + Present.Surface.Salinity.Lt.min + Present.Surface.Salinity.Mean + Present.Surface.Temperature.Lt.max + Present.Surface.Temperature.Lt.min + Present.Surface.Temperature.Mean,
    #BRT <- gbm(presence ~ X2100AOGCM.RCP85.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Current.Velocity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1,
               cv.folds = 100,
               interaction.depth = 3,
               distribution = "bernoulli",
               shrinkage = 0.01,
               verbose = FALSE,
               n.trees = 500,
               data=ENV_PA)
    eB <- evaluate(presence_locations.pres_test,background_macaronesia_test,BRT,env)
    tssB <- c(tssB,eB)
    tB <- threshold(eB,stat='prevalence')
    pB <- predict(env, BRT, ext=macaronesia.extent,type="response")>tB
    predictionsB <- c(predictionsB,pB)
    pB <- predict(future_env1, BRT, ext=macaronesia.extent,type="response")>tB
    predictionsB1 <- c(predictionsB1,pB)
    pB <- predict(future_env2, BRT, ext=macaronesia.extent,type="response")>tB
    predictionsB2 <- c(predictionsB2,pB)
    pB <- predict(future_env3, BRT, ext=macaronesia.extent,type="response")>tB
    predictionsB3 <- c(predictionsB3,pB)
    pB <- predict(future_env4, BRT, ext=macaronesia.extent,type="response")>tB
    predictionsB4 <- c(predictionsB4,pB)
    pB <- predict(future_env5, BRT, ext=macaronesia.extent,type="response")>tB
    predictionsB5 <- c(predictionsB5,pB)
    pB <- predict(future_env6, BRT, ext=macaronesia.extent,type="response")>tB
    predictionsB6 <- c(predictionsB6,pB)
    pB <- predict(future_env7, BRT, ext=macaronesia.extent,type="response")>tB
    predictionsB7 <- c(predictionsB7,pB)
    pB <- predict(future_env8, BRT, ext=macaronesia.extent,type="response")>tB
    predictionsB8 <- c(predictionsB8,pB)

    
    
    #    GAM <- gam(presence ~ poly(BO_sstmean,3)+poly(BO_dissox,3)+poly(BO2_tempmean_bdmin,3)+poly(BO2_salinitymean_bdmin,3)+poly(BO2_salinitymean_ss,3)+poly(BO2_carbonphytomean_ss,3)+poly(BO21_carbonphytomean_bdmin,3),
    GAM <- gam(presence ~ poly(Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1,3) + poly(Present.Benthic.Mean.Depth.Salinity.Lt.max,3) + poly(Present.Benthic.Mean.Depth.Salinity.Lt.min,3) + poly(Present.Benthic.Mean.Depth.Salinity.Mean,3) + poly(Present.Benthic.Mean.Depth.Temperature.Lt.max,3) + poly(Present.Benthic.Mean.Depth.Temperature.Lt.min,3) + poly(Present.Benthic.Mean.Depth.Temperature.Mean,3) + poly(Present.Surface.Current.Velocity.Mean.asc.BOv2_1,3) + poly(Present.Surface.Salinity.Lt.max,3) + poly(Present.Surface.Salinity.Lt.min,3) + poly(Present.Surface.Salinity.Mean,3) + poly(Present.Surface.Temperature.Lt.max,3) + poly(Present.Surface.Temperature.Lt.min,3) + poly(Present.Surface.Temperature.Mean,3),
    #GAM <- gam(presence ~ poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Current.Velocity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Salinity.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Salinity.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Salinity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1, 3),
               family=binomial(link='logit'),
               data=ENV_PA)
    ea <- evaluate(presence_locations.pres_test,background_macaronesia_test,GAM,env)
    tssa <- c(tssa,ea)
    ta <- threshold(ea,stat='prevalence')
    pa <- predict(env, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa <- c(predictionsa,pa)
    pa <- predict(future_env1, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa1 <- c(predictionsa1,pa)
    pa <- predict(future_env2, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa2 <- c(predictionsa2,pa)
    pa <- predict(future_env3, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa3 <- c(predictionsa3,pa)
    pa <- predict(future_env4, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa4 <- c(predictionsa4,pa)
    pa <- predict(future_env5, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa5 <- c(predictionsa5,pa)
    pa <- predict(future_env6, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa6 <- c(predictionsa6,pa)
    pa <- predict(future_env7, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa7 <- c(predictionsa7,pa)
    pa <- predict(future_env8, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa8 <- c(predictionsa8,pa)
    
    
    
    # #cat('neural')
    # NN <- nnet(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
    #            data=ENV_PA,
    #            size=10)
    # enn <- evaluate(presence_locations.pres_test,background_macaronesia_test,NN,env)
    # predictionsnn <- c(predictionsnn,predict(env, NN, ext=macaronesia.extent,type="response"))
    # tssnn <- c(tssnn,enn)
    
    
    #maxEnt <- maxent(env,presence_locations.pres_train,)
    maxEnt <- maxent(x=presence_absence_env,presence_absence,silent=TRUE)
    em <- evaluate(presence_locations.pres_test,background_macaronesia_test,maxEnt,env)
    tssmx <- c(tssmx,em)
    tm <- threshold(em,stat='prevalence')
    pm <- predict(env, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX <- c(predictionsmX,pm)
    pm <- predict(future_env1, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX1 <- c(predictionsmX1,pm)
    pm <- predict(future_env2, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX2 <- c(predictionsmX2,pm)
    pm <- predict(future_env3, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX3 <- c(predictionsmX3,pm)
    pm <- predict(future_env4, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX4 <- c(predictionsmX4,pm)
    pm <- predict(future_env5, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX5 <- c(predictionsmX5,pm)
    pm <- predict(future_env6, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX6 <- c(predictionsmX6,pm)
    pm <- predict(future_env7, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX7 <- c(predictionsmX7,pm)
    pm <- predict(future_env8, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX8 <- c(predictionsmX8,pm)
    
    
    
    #    glmm <- glm(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
    glmm <- glm(presence ~ Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + Present.Benthic.Mean.Depth.Salinity.Lt.max + Present.Benthic.Mean.Depth.Salinity.Lt.min + Present.Benthic.Mean.Depth.Salinity.Mean + Present.Benthic.Mean.Depth.Temperature.Lt.max + Present.Benthic.Mean.Depth.Temperature.Lt.min + Present.Benthic.Mean.Depth.Temperature.Mean + Present.Surface.Current.Velocity.Mean.asc.BOv2_1 + Present.Surface.Salinity.Lt.max + Present.Surface.Salinity.Lt.min + Present.Surface.Salinity.Mean + Present.Surface.Temperature.Lt.max + Present.Surface.Temperature.Lt.min + Present.Surface.Temperature.Mean,
    #glmm <- glm(presence ~ X2100AOGCM.RCP85.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Current.Velocity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1,
                family=binomial(link='logit'),
                data=ENV_PA)
    elm <- evaluate(presence_locations.pres_test,background_macaronesia_test,glmm,env)
    tssml <- c(tssml,elm)
    tlm <- threshold(elm,stat='prevalence')
    plm <- predict(env, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml <- c(predictionsml,plm)
    plm <- predict(future_env1, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml1 <- c(predictionsml1,plm)
    plm <- predict(future_env2, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml2 <- c(predictionsml2,plm)
    plm <- predict(future_env3, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml3 <- c(predictionsml3,plm)
    plm <- predict(future_env4, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml4 <- c(predictionsml4,plm)
    plm <- predict(future_env5, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml5 <- c(predictionsml5,plm)
    plm <- predict(future_env6, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml6 <- c(predictionsml6,plm)
    plm <- predict(future_env7, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml7 <- c(predictionsml7,plm)
    plm <- predict(future_env8, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml8 <- c(predictionsml8,plm)
    
    
    
    maha <- mahal(sp.env)
    ema <- evaluate(presence_locations.pres_test,background_macaronesia_test,maha,env)
    tssma <- c(tssma,ema)
    tma <-  threshold(ema,stat='prevalence')
    pma <- predict(env, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma <- c(predictionsma, pma)
    pma <- predict(future_env1, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma1 <- c(predictionsma1, pma)
    pma <- predict(future_env2, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma2 <- c(predictionsma2, pma)
    pma <- predict(future_env3, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma3 <- c(predictionsma3, pma)
    pma <- predict(future_env4, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma4 <- c(predictionsma4, pma)
    pma <- predict(future_env5, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma5 <- c(predictionsma5, pma)
    pma <- predict(future_env6, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma6 <- c(predictionsma6, pma)
    pma <- predict(future_env7, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma7 <- c(predictionsma7, pma)
    pma <- predict(future_env8, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma8 <- c(predictionsma8, pma)
    
    
    
    #d <- domain(env,presence_locations.pres_train)
    d <- domain(sp.env)
    ed <- evaluate(presence_locations.pres_test,background_macaronesia_test,d,env)
    tssd <- c(tssd,ed)
    td <- threshold(ed,stat='prevalence')
    pd <- predict(env, d, ext=macaronesia.extent,type="response")>td
    predictionsd <- c(predictionsd,pd)
    pd <- predict(future_env1, d, ext=macaronesia.extent,type="response")>td
    predictionsd1 <- c(predictionsd1,pd)
    pd <- predict(future_env2, d, ext=macaronesia.extent,type="response")>td
    predictionsd2 <- c(predictionsd2,pd)
    pd <- predict(future_env3, d, ext=macaronesia.extent,type="response")>td
    predictionsd3 <- c(predictionsd3,pd)
    pd <- predict(future_env4, d, ext=macaronesia.extent,type="response")>td
    predictionsd4 <- c(predictionsd4,pd)
    pd <- predict(future_env5, d, ext=macaronesia.extent,type="response")>td
    predictionsd5 <- c(predictionsd5,pd)
    pd <- predict(future_env6, d, ext=macaronesia.extent,type="response")>td
    predictionsd6 <- c(predictionsd6,pd)
    pd <- predict(future_env7, d, ext=macaronesia.extent,type="response")>td
    predictionsd7 <- c(predictionsd7,pd)
    pd <- predict(future_env8, d, ext=macaronesia.extent,type="response")>td
    predictionsd8 <- c(predictionsd8,pd)
    
    
    
    # bio <- bioclim(env,presence_locations.pres_train)
    bio <- bioclim(sp.env)
    eb <- evaluate(presence_locations.pres_test,background_macaronesia_test,bio,env)
    tssb <- c(tssb,eb)
    tb <- threshold(eb,stat='prevalence')
    pb <- predict(env, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb <- c(predictionsb,pb)
    pb <- predict(future_env1, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb1 <- c(predictionsb1,pb)
    pb <- predict(future_env2, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb2 <- c(predictionsb2,pb)
    pb <- predict(future_env3, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb3 <- c(predictionsb3,pb)
    pb <- predict(future_env4, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb4 <- c(predictionsb4,pb)
    pb <- predict(future_env5, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb5 <- c(predictionsb5,pb)
    pb <- predict(future_env6, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb6 <- c(predictionsb6,pb)
    pb <- predict(future_env7, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb7 <- c(predictionsb7,pb)
    pb <- predict(future_env8, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb8 <- c(predictionsb8,pb)
    
  }
  
  cat('BRT model\n')
  wB <- sapply(tssB, function(x) TSSweight(x))
  tssB <- weighted.mean(sapply(tssB, function(x) TSS(x)), wB, na.rm=TRUE)
  pB <- stack(predictionsB)
  pB <- weighted.mean(pB, wB, na.rm=TRUE)
  mods <- pB
  add_title <- 'BRT'
  m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  
  pB1 <- stack(predictionsB1)
  pB1 <- weighted.mean(pB1, wB, na.rm=TRUE)
  mods <- pB1
  m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP26_',add_title,'.html', sep=''))
  pB2 <- stack(predictionsB2)
  pB2 <- weighted.mean(pB2, wB, na.rm=TRUE)
  mods <- pB2
  m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP45_',add_title,'.html', sep=''))
  pB3 <- stack(predictionsB3)
  pB3 <- weighted.mean(pB3, wB, na.rm=TRUE)
  mods <- pB3
  m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP60_',add_title,'.html', sep=''))
  pB4 <- stack(predictionsB4)
  pB4 <- weighted.mean(pB4, wB, na.rm=TRUE)
  mods <- pB4
  m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP85_',add_title,'.html', sep=''))
  pB5 <- stack(predictionsB5)
  pB5 <- weighted.mean(pB5, wB, na.rm=TRUE)
  mods <- pB5
  m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP26_',add_title,'.html', sep=''))
  pB6 <- stack(predictionsB6)
  pB6 <- weighted.mean(pB6, wB, na.rm=TRUE)
  mods <- pB6
  m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP45_',add_title,'.html', sep=''))
  pB7 <- stack(predictionsB7)
  pB7 <- weighted.mean(pB7, wB, na.rm=TRUE)
  mods <- pB7
  m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP60_',add_title,'.html', sep=''))
  pB8 <- stack(predictionsB8)
  pB8 <- weighted.mean(pB8, wB, na.rm=TRUE)
  mods <- pB8
  m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP85_',add_title,'.html', sep=''))
  
  
  # cat("Neural Network model\n")
  # pN <- stack(predictionsnn)
  # wnn <- sapply(tssnn, function(x) TSSweight(x))
  # tssnn <- weighted.mean(sapply(tssnn, function(x) TSS(x)), wnn, na.rm=TRUE)
  # pN <- weighted.mean(pN, wnn, na.rm=TRUE)
  # mods <- pN
  # add_title <- 'Neural Network'
  # m <-  plotting(mods,species_name,tssnn,add_title,presence_locations)
  # saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))
  
  
  cat('GAM model\n')
  wa <- sapply(tssa, function(x) TSSweight(x))
  tssa <- weighted.mean(sapply(tssa, function(x) TSS(x)), wa, na.rm=TRUE)
  pa <- stack(predictionsa)
  pa <- weighted.mean(pa, wa, na.rm=TRUE)
  mods <- pa
  add_title <- 'GAM'
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  pa1 <- stack(predictionsa1)
  pa1 <- weighted.mean(pa1, wa, na.rm=TRUE)
  mods <- pa1
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP26_',add_title,'.html', sep=''))
  pa2 <- stack(predictionsa2)
  pa2 <- weighted.mean(pa2, wa, na.rm=TRUE)
  mods <- pa2
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP45_',add_title,'.html', sep=''))
  pa3 <- stack(predictionsa3)
  pa3 <- weighted.mean(pa3, wa, na.rm=TRUE)
  mods <- pa3
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP60_',add_title,'.html', sep=''))
  pa4 <- stack(predictionsa4)
  pa4 <- weighted.mean(pa4, wa, na.rm=TRUE)
  mods <- pa4
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP85_',add_title,'.html', sep=''))
  pa5 <- stack(predictionsa5)
  pa5 <- weighted.mean(pa5, wa, na.rm=TRUE)
  mods <- pa5
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP26_',add_title,'.html', sep=''))
  pa6 <- stack(predictionsa6)
  pa6 <- weighted.mean(pa6, wa, na.rm=TRUE)
  mods <- pa6
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP45_',add_title,'.html', sep=''))
  pa7 <- stack(predictionsa7)
  pa7 <- weighted.mean(pa7, wa, na.rm=TRUE)
  mods <- pa7
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP60_',add_title,'.html', sep=''))
  pa8 <- stack(predictionsa8)
  pa8 <- weighted.mean(pa8, wa, na.rm=TRUE)
  mods <- pa8
  m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP85_',add_title,'.html', sep=''))
  
  
  cat('GLM model\n')
  pml <- stack(predictionsml)
  wml <- sapply(tssml, function(x) TSSweight(x))
  tssml <- weighted.mean(sapply(tssml, function(x) TSS(x)), wml, na.rm=TRUE)
  pml <- weighted.mean(pml, wml, na.rm=TRUE)
  mods <- pml
  add_title <- 'GLM'
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  pml1 <-stack(predictionsml1)
  pml1 <- weighted.mean(pml1, wml, na.rm=TRUE)
  mods <- pml1
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP26_',add_title,'.html', sep=''))
  pml2 <-stack(predictionsml2)
  pml2 <- weighted.mean(pml2, wml, na.rm=TRUE)
  mods <- pml2
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP45_',add_title,'.html', sep=''))
  pml3 <-stack(predictionsml3)
  pml3 <- weighted.mean(pml3, wml, na.rm=TRUE)
  mods <- pml3
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP60_',add_title,'.html', sep=''))
  pml4 <-stack(predictionsml4)
  pml4 <- weighted.mean(pml4, wml, na.rm=TRUE)
  mods <- pml4
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP85_',add_title,'.html', sep=''))
  pml5 <-stack(predictionsml5)
  pml5 <- weighted.mean(pml5, wml, na.rm=TRUE)
  mods <- pml5
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP26_',add_title,'.html', sep=''))
  pml6 <-stack(predictionsml6)
  pml6 <- weighted.mean(pml6, wml, na.rm=TRUE)
  mods <- pml6
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP45_',add_title,'.html', sep=''))
  pml7 <-stack(predictionsml7)
  pml7 <- weighted.mean(pml7, wml, na.rm=TRUE)
  mods <- pml7
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP60_',add_title,'.html', sep=''))
  pml8 <-stack(predictionsml8)
  pml8 <- weighted.mean(pml8, wml, na.rm=TRUE)
  mods <- pml8
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP85_',add_title,'.html', sep=''))


  
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
  saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  pma1 <- stack(predictionsma1)
  pma1 <- weighted.mean(pma1, wma, na.rm=TRUE)
  mods <- pma1
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP26_',add_title,'.html', sep=''))
  pma2 <- stack(predictionsma2)
  pma2 <- weighted.mean(pma2, wma, na.rm=TRUE)
  mods <- pma2
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP45_',add_title,'.html', sep=''))
  pma3 <- stack(predictionsma3)
  pma3 <- weighted.mean(pma3, wma, na.rm=TRUE)
  mods <- pma3
  saveWidget(m, file=paste(species_name,'_2050.RCP60_',add_title,'.html', sep=''))
  pma4 <- stack(predictionsma4)
  pma4 <- weighted.mean(pma4, wma, na.rm=TRUE)
  mods <- pma4
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP85_',add_title,'.html', sep=''))
  pma5 <- stack(predictionsma5)
  pma5 <- weighted.mean(pma5, wma, na.rm=TRUE)
  mods <- pma5
  saveWidget(m, file=paste(species_name,'_2100.RCP26_',add_title,'.html', sep=''))
  pma6 <- stack(predictionsma6)
  pma6 <- weighted.mean(pma6, wma, na.rm=TRUE)
  mods <- pma6
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP45_',add_title,'.html', sep=''))
  pma7 <- stack(predictionsma7)
  pma7 <- weighted.mean(pma7, wma, na.rm=TRUE)
  mods <- pma7
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP60_',add_title,'.html', sep=''))
  pma8 <- stack(predictionsma8)
  pma8 <- weighted.mean(pma8, wma, na.rm=TRUE)
  mods <- pma8
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP85_',add_title,'.html', sep=''))



  cat('Bioclim model\n')#filter of .2 TSS is removing the layers
  pb <- stack(predictionsb)
  wb <- sapply(tssb, function(x) TSSweight(x))
  #tssb <- sapply(tssb, function(x) TSS(x))
  tssb <- weighted.mean(sapply(tssb, function(x) TSS(x)), wb, na.rm=TRUE)
  pb <- weighted.mean(pb, wb, na.rm=TRUE)
  mods <- pb
  add_title <- 'Bioclim'
  m <-  plotting(mods,species_name,tssb,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  pb <- stack(predictionsb)
  pb <- weighted.mean(pb, wb, na.rm=TRUE)
  mods <- pb
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  pb1 <- stack(predictionsb1)
  pb1 <- weighted.mean(pb1, wb, na.rm=TRUE)
  mods <- pb1
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP26_',add_title,'.html', sep=''))
  pb2 <- stack(predictionsb2)
  pb2 <- weighted.mean(pb2, wb, na.rm=TRUE)
  mods <- pb2
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP45_',add_title,'.html', sep=''))
  pb3 <- stack(predictionsb3)
  pb3 <- weighted.mean(pb3, wb, na.rm=TRUE)
  mods <- pb3
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP60_',add_title,'.html', sep=''))
  pb4 <- stack(predictionsb4)
  pb4 <- weighted.mean(pb4, wb, na.rm=TRUE)
  mods <- pb4
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP85_',add_title,'.html', sep=''))
  pb5 <- stack(predictionsb5)
  pb5 <- weighted.mean(pb5, wb, na.rm=TRUE)
  mods <- pb5
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP26_',add_title,'.html', sep=''))
  pb6 <- stack(predictionsb6)
  pb6 <- weighted.mean(pb6, wb, na.rm=TRUE)
  mods <- pb6
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP45_',add_title,'.html', sep=''))
  pb7 <- stack(predictionsb7)
  pb7 <- weighted.mean(pb7, wb, na.rm=TRUE)
  mods <- pb7
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP60_',add_title,'.html', sep=''))
  pb8 <- stack(predictionsb8)
  pb8 <- weighted.mean(pb8, wb, na.rm=TRUE)
  mods <- pb8
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP85_',add_title,'.html', sep=''))

  
  cat('Domain model\n')
  pd <- stack(predictionsd)
  wd <- sapply(tssd, function(x) TSSweight(x))
  #tssd <- sapply(tssd, function(x) TSS(x))
  tssd <- weighted.mean(sapply(tssd, function(x) TSS(x)), wd, na.rm=TRUE)
  pd <- weighted.mean(pd, wd, na.rm=TRUE)
  mods <- pd
  add_title <- 'Domain'
  m <-  plotting(mods,species_name,tssd,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  pd1 <- stack(predictionsd1)
  pd1 <- weighted.mean(pd1, wd, na.rm=TRUE)
  mods <- pd1
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP26_',add_title,'.html', sep=''))
  pd2 <- stack(predictionsd2)
  pd2 <- weighted.mean(pd2, wd, na.rm=TRUE)
  mods <- pd2
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP45_',add_title,'.html', sep=''))
  pd3 <- stack(predictionsd3)
  pd3 <- weighted.mean(pd3, wd, na.rm=TRUE)
  mods <- pd3
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP60_',add_title,'.html', sep=''))
  pd4 <- stack(predictionsd4)
  pd4 <- weighted.mean(pd4, wd, na.rm=TRUE)
  mods <- pd4
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP85_',add_title,'.html', sep=''))
  pd5 <- stack(predictionsd5)
  pd5 <- weighted.mean(pd5, wd, na.rm=TRUE)
  mods <- pd5
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP26_',add_title,'.html', sep=''))
  pd6 <- stack(predictionsd6)
  pd6 <- weighted.mean(pd6, wd, na.rm=TRUE)
  mods <- pd6
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP45_',add_title,'.html', sep=''))
  pd7 <- stack(predictionsd7)
  pd7 <- weighted.mean(pd7, wd, na.rm=TRUE)
  mods <- pd7
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP60_',add_title,'.html', sep=''))
  pd8 <- stack(predictionsd8)
  pd8 <- weighted.mean(pd8, wd, na.rm=TRUE)
  mods <- pd8
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP85_',add_title,'.html', sep=''))

  
  cat('MaxEnt model\n')
  pmX <- stack(predictionsmX)
  wmX <- sapply(tssmx, function(x) TSSweight(x))
  #tssmx <- sapply(tssmx, function(x) TSS(x))
  tssmx <- weighted.mean(sapply(tssmx, function(x) TSS(x)), wmX, na.rm=TRUE)
  pmX <- weighted.mean(pmX, wmX, na.rm=TRUE)
  mods <- pmX
  add_title <- 'MaxEnt'
  m <-  plotting(mods,species_name,tssd,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  pmX1 <- stack(predictionsmX1)
  pmX1 <- weighted.mean(pmX1, wmX, na.rm=TRUE)
  mods <- pmX1
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP26_',add_title,'.html', sep=''))
  pmX2 <- stack(predictionsmX2)
  pmX2 <- weighted.mean(pmX2, wmX, na.rm=TRUE)
  mods <- pmX2
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP45_',add_title,'.html', sep=''))
  pmX3 <- stack(predictionsmX3)
  pmX3 <- weighted.mean(pmX3, wmX, na.rm=TRUE)
  mods <- pmX3
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP60_',add_title,'.html', sep=''))
  pmX4 <- stack(predictionsmX4)
  pmX4 <- weighted.mean(pmX4, wmX, na.rm=TRUE)
  mods <- pmX4
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP85_',add_title,'.html', sep=''))
  pmX5 <- stack(predictionsmX5)
  pmX5 <- weighted.mean(pmX5, wmX, na.rm=TRUE)
  mods <- pmX5
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP26_',add_title,'.html', sep=''))
  pmX6 <- stack(predictionsmX6)
  pmX6 <- weighted.mean(pmX6, wmX, na.rm=TRUE)
  mods <- pmX6
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP45_',add_title,'.html', sep=''))
  pmX7 <- stack(predictionsmX7)
  pmX7 <- weighted.mean(pmX7, wmX, na.rm=TRUE)
  mods <- pmX7
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP60_',add_title,'.html', sep=''))
  pmX8 <- stack(predictionsmX8)
  pmX8 <- weighted.mean(pmX8, wmX, na.rm=TRUE)
  mods <- pmX8
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP85_',add_title,'.html', sep=''))

  
  tsss <- c(tssB, tssa, tssml, tssma, tssb, tssd, tssmx)#tssnn
  cat('Ensemble\n')
  models <- stack(pB, pa, pml, pma, pb, pd, pmX)#pN
  names(models) <- c("BRT","GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  #w <- sapply(list(ea, elm, ema, eb, ed, em), function(x) TSSweight(x))
  #Mean <- weighted.mean( models[[c("GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], w)
  ws <- weights_tsss(tsss)
  Mean <- weighted.mean( models[[c("BRT", "GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  tsss <- weighted.mean(tsss,ws, na.rm=TRUE)
  mods <- Mean
  #mods <- mask(mods,mask_depth)
  add_title <- 'Ensemble'
  m <-  plotting(mods,species_name,tssd,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  models <- stack(pB1, pa1, pml1, pma1, pb1, pd1, pmX1)
  names(models) <- c("BRT","GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  Mean <- weighted.mean( models[[c("BRT", "GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  mods <- Mean
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP26_',add_title,'.html', sep=''))
  models <- stack(pB2, pa2, pml2, pma2, pb2, pd2, pmX2)
  names(models) <- c("BRT","GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  Mean <- weighted.mean( models[[c("BRT", "GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  mods <- Mean
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP45_',add_title,'.html', sep=''))
  models <- stack(pB3, pa3, pml3, pma3, pb3, pd3, pmX3)
  names(models) <- c("BRT","GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  Mean <- weighted.mean( models[[c("BRT", "GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  mods <- Mean
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP60_',add_title,'.html', sep=''))
  models <- stack(pB4, pa4, pml4, pma4, pb4, pd4, pmX4)
  names(models) <- c("BRT","GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  Mean <- weighted.mean( models[[c("BRT", "GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  mods <- Mean
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2050.RCP85_',add_title,'.html', sep=''))
  models <- stack(pB5, pa5, pml5, pma5, pb5, pd5, pmX5)
  names(models) <- c("BRT","GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  Mean <- weighted.mean( models[[c("BRT", "GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  mods <- Mean
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP26_',add_title,'.html', sep=''))
  models <- stack(pB6, pa6, pml6, pma6, pb6, pd6, pmX6)
  names(models) <- c("BRT","GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  Mean <- weighted.mean( models[[c("BRT", "GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  mods <- Mean
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP45_',add_title,'.html', sep=''))
  models <- stack(pB7, pa7, pml7, pma7, pb7, pd7, pmX7)
  names(models) <- c("BRT","GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  Mean <- weighted.mean( models[[c("BRT", "GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  mods <- Mean
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP60_',add_title,'.html', sep=''))
  models <- stack(pB8, pa8, pml8, pma8, pb8, pd8, pmX8)
  names(models) <- c("BRT","GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")
  Mean <- weighted.mean( models[[c("BRT", "GAM", "GLM", "mahalanobis", "bioclim", "Domain", "maxent")]], ws, na.rm=TRUE)
  mods <- Mean
  m <-  plotting(mods,species_name,tssml,add_title,presence_locations)
  saveWidget(m, file=paste(species_name,'_2100.RCP85_',add_title,'.html', sep=''))
}

#lst <- list.files(pattern='asc$')
lste <- list.files(pattern='^Present')
env <- stack(lste)
mask_depth <- load_layers(list_layers( datasets="Bio-ORACLE" )[25,]$layer_code)
mask_depth <- mask_depth[[1]]
mask_depth[mask_depth>-10] <- NA
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
oegopsid_locations <- select(oegopsids,
                             longitude,
                             latitude)

background <- randomPoints(env, n=1000, extf = 1.25)
colnames(background) = c('longitude', 'latitude')
group <- kfold(background, 4)
background_train <- background[group != 1, ]
background_train.env <- extract(env, background_train)
background_test <- background[group == 1, ]
#total.extent <- extent(-35, -10, 12, 42)#lon,lon,lat,lat

illex <- oegopsids[oegopsids$latitude<(12.0*-1),]
illex <- select(illex[illex$species=='Illex argentinus',],
                species,
                longitude,
                latitude)
main('Illex argentinus',
     illex,env,
     future_env1,
     future_env2,
     future_env3,
     future_env4,
     future_env5,
     future_env6,
     future_env7,
     future_env8,background_test,background_train.env,extent(-100, -23, -80, 10))

nsloanii <- select(oegopsids[oegopsids$species=='Nototodarus sloanii',],
                   species,
                   longitude,
                   latitude)
main180('Nototodarus sloanii',nsloanii,env,future_env,background_test,background_train.env,extent(140, 180, -60, -10),extent(-180+0.008333333, -150, -60, -10))

tpacificus<- select(oegopsids[oegopsids$species=='Todarodes pacificus',],
                    species,
                    longitude,
                    latitude)
main180('Todarodes pacificus',tpacificus,env,future_env,background_test,background_train.env,extent(105, 180, 18, 80),extent(-180+0.008333333, -130, 18, 80))

bmagister<- select(oegopsids[oegopsids$species=='Berryteuthis magister',],
                   species,
                   longitude,
                   latitude)
main180('Berryteuthis magister',bmagister,env,future_env,background_test,background_train.env,extent(126, 180, 33, 80),extent(-180+0.008333333, -120, 33, 80))

