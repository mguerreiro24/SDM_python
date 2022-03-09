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
  background_macaronesia_test,
  background_macaronesia_train.env,
  macaronesia.extent,
  macaronesia.extent2,
  biomasses){
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


  BRT <- gbm(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
  cv.folds = 100,
  interaction.depth = 3,
  distribution = "bernoulli",
  shrinkage = 0.01,
  verbose = FALSE,
  n.trees = 500,
  data=ENV_PA)
  eB <- evaluate(presence_locations.pres_test,background_macaronesia_test,BRT,env)
  tB <- threshold(eB,stat='prevalence')
  pB <- predict(env, BRT, ext=macaronesia.extent,type="response")>tB 
  predictionsB <- c(predictionsB,pB)
  pB <- predict(env, BRT, ext=macaronesia.extent2,type="response")>tB 
  predictionsB2 <- c(predictionsB2,pB)
  tssB <- c(tssB,eB)


  GAM <- gam(presence ~ poly(BO_sstmean,3)+poly(BO_dissox,3)+poly(BO2_tempmean_bdmin,3)+poly(BO2_salinitymean_bdmin,3)+poly(BO2_salinitymean_ss,3)+poly(BO2_carbonphytomean_ss,3)+poly(BO21_carbonphytomean_bdmin,3),
              family=binomial(link='logit'),
              data=ENV_PA)
  ea <- evaluate(presence_locations.pres_test,background_macaronesia_test,GAM,env)
  ta <- threshold(ea,stat='prevalence')
  pa <- predict(env, GAM, ext=macaronesia.extent,type="response")>ta
  predictionsa <- c(predictionsa,pa)
  pa <- predict(env, GAM, ext=macaronesia.extent2,type="response")>ta
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
  pm <- predict(env, maxEnt, ext=macaronesia.extent,type="response")>tm
  predictionsmX <- c(predictionsmX,pm)
  pm <- predict(env, maxEnt, ext=macaronesia.extent2,type="response")>tm
  predictionsmX2 <- c(predictionsmX2,pm)
  tssmx <- c(tssmx,em)
  
  
  glmm <- glm(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
              family=binomial(link='logit'),
              data=ENV_PA)
  elm <- evaluate(presence_locations.pres_test,background_macaronesia_test,glmm,env)
  tlm <- threshold(elm,stat='prevalence')
  plm <- predict(env, glmm, ext=macaronesia.extent,type="response")>tlm
  predictionsml <- c(predictionsml,plm)
  plm <- predict(env, glmm, ext=macaronesia.extent2,type="response")>tlm
  predictionsml2 <- c(predictionsml2,plm)
  tssml <- c(tssml,elm)

  
  maha <- mahal(sp.env)
  ema <- evaluate(presence_locations.pres_test,background_macaronesia_test,maha,env)
  tma <-  threshold(ema,stat='prevalence')
  pma <- predict(env, maha, ext=macaronesia.extent,type="response")>tma
  predictionsma <- c(predictionsma, pma)
  pma <- predict(env, maha, ext=macaronesia.extent2,type="response")>tma
  predictionsma2 <- c(predictionsma2, pma)
  tssma <- c(tssma,ema)
  
  
  #d <- domain(env,presence_locations.pres_train)
  d <- domain(sp.env)
  ed <- evaluate(presence_locations.pres_test,background_macaronesia_test,d,env)
  td <- threshold(ed,stat='prevalence')
  pd <- predict(env, d, ext=macaronesia.extent,type="response")>td
  predictionsd <- c(predictionsd,pd)
  pd <- predict(env, d, ext=macaronesia.extent2,type="response")>td
  predictionsd2 <- c(predictionsd2,pd)
  tssd <- c(tssd,ed)
  
  
  # bio <- bioclim(env,presence_locations.pres_train)
  bio <- bioclim(sp.env)
  eb <- evaluate(presence_locations.pres_test,background_macaronesia_test,bio,env)
  tb <- threshold(eb,stat='prevalence')
  pb <- predict(env, bio, ext=macaronesia.extent,type="response")>tb
  predictionsb <- c(predictionsb,pb)
  pb <- predict(env, bio, ext=macaronesia.extent2,type="response")>tb
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
                 background_macaronesia_test,
                 background_macaronesia_train.env,
                 macaronesia.extent,
                 biomasses){
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
#    BRT <- gbm(presence ~ Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + Present.Benthic.Mean.Depth.Salinity.Lt.max + Present.Benthic.Mean.Depth.Salinity.Lt.min + Present.Benthic.Mean.Depth.Salinity.Mean + Present.Benthic.Mean.Depth.Temperature.Lt.max + Present.Benthic.Mean.Depth.Temperature.Lt.min + Present.Benthic.Mean.Depth.Temperature.Mean + Present.Surface.Current.Velocity.Mean.asc.BOv2_1 + Present.Surface.Salinity.Lt.max + Present.Surface.Salinity.Lt.min + Present.Surface.Salinity.Mean + Present.Surface.Temperature.Lt.max + Present.Surface.Temperature.Lt.min + Present.Surface.Temperature.Mean,
    BRT <- gbm(presence ~ X2100AOGCM.RCP85.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Current.Velocity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1,
               cv.folds = 100,
               interaction.depth = 3,
               distribution = "bernoulli",
               shrinkage = 0.01,
               verbose = FALSE,
               n.trees = 500,
               data=ENV_PA)
    eB <- evaluate(presence_locations.pres_test,background_macaronesia_test,BRT,env)
    tB <- threshold(eB,stat='prevalence')
    pB <- predict(env, BRT, ext=macaronesia.extent,type="response")>tB 
    predictionsB <- c(predictionsB,pB)
    tssB <- c(tssB,eB)
    
    
#    GAM <- gam(presence ~ poly(BO_sstmean,3)+poly(BO_dissox,3)+poly(BO2_tempmean_bdmin,3)+poly(BO2_salinitymean_bdmin,3)+poly(BO2_salinitymean_ss,3)+poly(BO2_carbonphytomean_ss,3)+poly(BO21_carbonphytomean_bdmin,3),
#    GAM <- gam(presence ~ poly(Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1,3) + poly(Present.Benthic.Mean.Depth.Salinity.Lt.max,3) + poly(Present.Benthic.Mean.Depth.Salinity.Lt.min,3) + poly(Present.Benthic.Mean.Depth.Salinity.Mean,3) + poly(Present.Benthic.Mean.Depth.Temperature.Lt.max,3) + poly(Present.Benthic.Mean.Depth.Temperature.Lt.min,3) + poly(Present.Benthic.Mean.Depth.Temperature.Mean,3) + poly(Present.Surface.Current.Velocity.Mean.asc.BOv2_1,3) + poly(Present.Surface.Salinity.Lt.max,3) + poly(Present.Surface.Salinity.Lt.min,3) + poly(Present.Surface.Salinity.Mean,3) + poly(Present.Surface.Temperature.Lt.max,3) + poly(Present.Surface.Temperature.Lt.min,3) + poly(Present.Surface.Temperature.Mean,3),
    GAM <- gam(presence ~ poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Current.Velocity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Salinity.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Salinity.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Salinity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1, 3),
               family=binomial(link='logit'),
               data=ENV_PA)
    ea <- evaluate(presence_locations.pres_test,background_macaronesia_test,GAM,env)
    ta <- threshold(ea,stat='prevalence')
    pa <- predict(env, GAM, ext=macaronesia.extent,type="response")>ta
    predictionsa <- c(predictionsa,pa)
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
    pm <- predict(env, maxEnt, ext=macaronesia.extent,type="response")>tm
    predictionsmX <- c(predictionsmX,pm)
    tssmx <- c(tssmx,em)
    
    
#    glmm <- glm(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
#    glmm <- glm(presence ~ Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + Present.Benthic.Mean.Depth.Salinity.Lt.max + Present.Benthic.Mean.Depth.Salinity.Lt.min + Present.Benthic.Mean.Depth.Salinity.Mean + Present.Benthic.Mean.Depth.Temperature.Lt.max + Present.Benthic.Mean.Depth.Temperature.Lt.min + Present.Benthic.Mean.Depth.Temperature.Mean + Present.Surface.Current.Velocity.Mean.asc.BOv2_1 + Present.Surface.Salinity.Lt.max + Present.Surface.Salinity.Lt.min + Present.Surface.Salinity.Mean + Present.Surface.Temperature.Lt.max + Present.Surface.Temperature.Lt.min + Present.Surface.Temperature.Mean,
    glmm <- glm(presence ~ X2100AOGCM.RCP85.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Current.Velocity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1,
                family=binomial(link='logit'),
                data=ENV_PA)
    elm <- evaluate(presence_locations.pres_test,background_macaronesia_test,glmm,env)
    tlm <- threshold(elm,stat='prevalence')
    plm <- predict(env, glmm, ext=macaronesia.extent,type="response")>tlm
    predictionsml <- c(predictionsml,plm)
    tssml <- c(tssml,elm)
    
    
    maha <- mahal(sp.env)
    ema <- evaluate(presence_locations.pres_test,background_macaronesia_test,maha,env)
    tma <-  threshold(ema,stat='prevalence')
    pma <- predict(env, maha, ext=macaronesia.extent,type="response")>tma
    predictionsma <- c(predictionsma, pma)
    tssma <- c(tssma,ema)
    
    
    #d <- domain(env,presence_locations.pres_train)
    d <- domain(sp.env)
    ed <- evaluate(presence_locations.pres_test,background_macaronesia_test,d,env)
    td <- threshold(ed,stat='prevalence')
    pd <- predict(env, d, ext=macaronesia.extent,type="response")>td
    predictionsd <- c(predictionsd,pd)
    tssd <- c(tssd,ed)
    
    
    # bio <- bioclim(env,presence_locations.pres_train)
    bio <- bioclim(sp.env)
    eb <- evaluate(presence_locations.pres_test,background_macaronesia_test,bio,env)
    tb <- threshold(eb,stat='prevalence')
    pb <- predict(env, bio, ext=macaronesia.extent,type="response")>tb
    predictionsb <- c(predictionsb,pb)
    tssb <- c(tssb,eb)
  }
  
  cat('BRT model\n')
  pB <- stack(predictionsB)
  wB <- sapply(tssB, function(x) TSSweight(x))
  tssB <- weighted.mean(sapply(tssB, function(x) TSS(x)), wB, na.rm=TRUE)
  pB <- weighted.mean(pB, wB, na.rm=TRUE)
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
  
  
  cat('Bioclim model\n')#filter of .2 TSS is removing the layers
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

#lst <- list.files(pattern='asc$')
lst <- list.files(pattern='^2100AOGCM.RCP85.')
env <- stack(lst)
#env <- load_layers(list_layers( datasets="Bio-ORACLE" )[c(25,21,12,158,212,310,298,335),]$layer_code)
#25,21,12,69,119,158,212,310,298,335
#bathi, SST,O2, "min deep O2",light bottom, temp bottom, sal bottom, sal surface,biomass,biomass deep
mask_depth <- load_layers(list_layers( datasets="Bio-ORACLE" )[25,]$layer_code)
mask_depth <- mask_depth[[1]]
#mask_depth <- crop(mask_depth,macaronesia.extent)
mask_depth[mask_depth>-10] <- NA
crs(env) <- crs(mask_depth)
env <- mask(env,mask_depth)
#macaronesia.extent <- extent(-35, -10, 12, 42)#lon,lon,lat,lat
#macaronesia.extent <- extent(-55, 0, -12, 50)
#macaronesia.env <- crop(env,macaronesia.extent)




oegopsids <- read.csv("occorrences.csv")
oegopsid_locations <- select(oegopsids,
                             longitude,
                             latitude)
biomasses <- read.csv("biomasses.csv")
# squid_occorrences <-  leaflet()%>%
#   addProviderTiles(provider = 'Esri.OceanBasemap') %>%
#   addScaleBar()%>%
#   addCircles(oegopsid_locations$longitude,oegopsid_locations$latitude)
# #saveWidget(m, file=paste(name,add_title,'.html',sep = ''))
# squid_occorrences
# 
# 
# oegopsid_locations <- spThin(oegopsid_locations,dist=500)[[1]]
# squid_occorrences <-  leaflet()%>%
#   addProviderTiles(provider = 'Esri.OceanBasemap') %>%
#   addScaleBar()%>%
#   addCircles(oegopsid_locations$x,oegopsid_locations$y)
#saveWidget(m, file=paste(name,add_title,'.html',sep = ''))
#squid_occorrences
#HB <- species_locations(oegopsids,'Histioteuthis bonnellii')
#squid_occorrences <-  leaflet()%>%
#  addProviderTiles(provider = 'Esri.OceanBasemap') %>%
#  addScaleBar()%>%
#  addCircles(HB$longitude,HB$latitude)
#saveWidget(m, file=paste(name,add_title,'.html',sep = ''))
#squid_occorrences



# background_macaronesia <- randomPoints(env, n=1000, ext=macaronesia.extent, extf = 1.25)
# 
# colnames(background_macaronesia) = c('longitude', 'latitude')
# group <- kfold(background_macaronesia, 4)
# background_macaronesia_train <- background_macaronesia[group != 1, ]
# background_macaronesia_train.env <- extract(env, background_macaronesia_train)
# background_macaronesia_test <- background_macaronesia[group == 1, ]


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
main('Illex argentinus',illex,env,background_test,background_train.env,extent(-100, -23, -80, 10),biomasses)

nsloanii <- select(oegopsids[oegopsids$species=='Nototodarus sloanii',],
                   species,
                   longitude,
                   latitude)
main180('Nototodarus sloanii',nsloanii,env,background_test,background_train.env,extent(140, 180, -60, -10),extent(-180+0.008333333, -150, -60, -10),biomasses)

tpacificus<- select(oegopsids[oegopsids$species=='Todarodes pacificus',],
                    species,
                    longitude,
                    latitude)
main180('Todarodes pacificus',tpacificus,env,background_test,background_train.env,extent(105, 180, 18, 80),extent(-180+0.008333333, -130, 18, 80),biomasses)

bmagister<- select(oegopsids[oegopsids$species=='Berryteuthis magister',],
                    species,
                    longitude,
                    latitude)
main180('Berryteuthis magister',bmagister,env,background_test,background_train.env,extent(126, 180, 33, 80),extent(-180+0.008333333, -120, 33, 80),biomasses)


#' 'Abralia veranyi',#)){#unique(oegopsids$species)
#'                        'Chiroteuthis veranii',
#'                        'Cranchia scabra',
#'                        'Histioteuthis bonnellii',
#'                        'Histioteuthis reversa',
#'                        'Liocranchia reinhardtii',
#'                        'Ommastrephes caroli',
#'                        'Onychoteuthis banksii',
#'                        #'Onykia robsoni',
#'                        'Pterygioteuthis gemmata',
#'                        #'Pterygioteuthis giardi',#not working, perhaps no need to remove NAs?
#'                        'Pyroteuthis margaritifera',
#'                        'Sthenoteuthis pterops',
#'                        #'Teuthowenia megalops',
#'                        #'Teuthowenia pellucida',
#'                        #'Todarodes filippovae',
#'                        'Todarodes sagittatus',
#'                        'Todaropsis eblanae'


species_name <- 'Liocranchia reinhardtii'

main('Liocranchia reinhardtii',oegopsids,env,background_macaronesia_test,background_macaronesia_train.env,macaronesia.extent,biomasses)


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
