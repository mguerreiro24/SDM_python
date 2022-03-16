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
  predictions <- list(
    predictionsB = list(
      predictions = list(c(),c()),
      predictions1 = list(c(),c()),
      predictions2 = list(c(),c()),
      predictions3 = list(c(),c()),
      predictions4 = list(c(),c()),
      predictions5 = list(c(),c()),
      predictions6 = list(c(),c()),
      predictions7 = list(c(),c()),
      predictions8 = list(c(),c())),
    predictionsa = list(
      predictions = list(c(),c()),
      predictions1 = list(c(),c()),
      predictions2 = list(c(),c()),
      predictions3 = list(c(),c()),
      predictions4 = list(c(),c()),
      predictions5 = list(c(),c()),
      predictions6 = list(c(),c()),
      predictions7 = list(c(),c()),
      predictions8 = list(c(),c())),
    predictionsmX = list(
      predictions = list(c(),c()),
      predictions1 = list(c(),c()),
      predictions2 = list(c(),c()),
      predictions3 = list(c(),c()),
      predictions4 = list(c(),c()),
      predictions5 = list(c(),c()),
      predictions6 = list(c(),c()),
      predictions7 = list(c(),c()),
      predictions8 = list(c(),c())),
    predictionsml = list(
      predictions = list(c(),c()),
      predictions1 = list(c(),c()),
      predictions2 = list(c(),c()),
      predictions3 = list(c(),c()),
      predictions4 = list(c(),c()),
      predictions5 = list(c(),c()),
      predictions6 = list(c(),c()),
      predictions7 = list(c(),c()),
      predictions8 = list(c(),c())),
    predictionsd = list(
      predictions = list(c(),c()),
      predictions1 = list(c(),c()),
      predictions2 = list(c(),c()),
      predictions3 = list(c(),c()),
      predictions4 = list(c(),c()),
      predictions5 = list(c(),c()),
      predictions6 = list(c(),c()),
      predictions7 = list(c(),c()),
      predictions8 = list(c(),c())),
    predictionsb = list(
      predictions = list(c(),c()),
      predictions1 = list(c(),c()),
      predictions2 = list(c(),c()),
      predictions3 = list(c(),c()),
      predictions4 = list(c(),c()),
      predictions5 = list(c(),c()),
      predictions6 = list(c(),c()),
      predictions7 = list(c(),c()),
      predictions8 = list(c(),c()))
    )
  tss <- list(
    tssB <- c(),
    tssa <- c(),
    tssmX <- c(),
    tssml <- c(),
    tssd <- c(),
    tssb <- c()
    )#True Skill score

  threshes <- list(
    B <- c(),
    a <- c(),
    mX <- c(),
    ml <- c(),
    d <- c(),
    b <- c()
    )#thresholds
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
    
     
    #BRT <- gbm(presence ~ Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + Present.Benthic.Mean.Depth.Salinity.Lt.max + Present.Benthic.Mean.Depth.Salinity.Lt.min + Present.Benthic.Mean.Depth.Salinity.Mean + Present.Benthic.Mean.Depth.Temperature.Lt.max + Present.Benthic.Mean.Depth.Temperature.Lt.min + Present.Benthic.Mean.Depth.Temperature.Mean + Present.Surface.Current.Velocity.Mean.asc.BOv2_1 + Present.Surface.Salinity.Lt.max + Present.Surface.Salinity.Lt.min + Present.Surface.Salinity.Mean + Present.Surface.Temperature.Lt.max + Present.Surface.Temperature.Lt.min + Present.Surface.Temperature.Mean,
    #            cv.folds = 100,
    #            interaction.depth = 3,
    #            distribution = "bernoulli",
    #            shrinkage = 0.01,
    #            verbose = FALSE,
    #            n.trees = 500,
    #            data=ENV_PA)
    #eB <- evaluate(presence_locations.pres_test,background_macaronesia_test,BRT,env)
    #tss$tssB <- c(tss$tssB,eB)
    #tB <- threshold(eB,stat='prevalence')
    #for (i in 1:length(envs)){
    #  predictions$predictionsB[[i]][[1]] <- c(predictions$predictionsB[[i]][[1]],predict(envs[[i]], BRT, ext=macaronesia.extent,type="response"))
    #  predictions$predictionsB[[i]][[2]] <- c(predictions$predictionsB[[i]][[2]],predict(envs[[i]], BRT, ext=macaronesia.extent2,type="response"))
    #}



    #GAM <- gam(presence ~ poly(Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1,3) + poly(Present.Benthic.Mean.Depth.Salinity.Lt.max,3) + poly(Present.Benthic.Mean.Depth.Salinity.Lt.min,3) + poly(Present.Benthic.Mean.Depth.Salinity.Mean,3) + poly(Present.Benthic.Mean.Depth.Temperature.Lt.max,3) + poly(Present.Benthic.Mean.Depth.Temperature.Lt.min,3) + poly(Present.Benthic.Mean.Depth.Temperature.Mean,3) + poly(Present.Surface.Current.Velocity.Mean.asc.BOv2_1,3) + poly(Present.Surface.Salinity.Lt.max,3) + poly(Present.Surface.Salinity.Lt.min,3) + poly(Present.Surface.Salinity.Mean,3) + poly(Present.Surface.Temperature.Lt.max,3) + poly(Present.Surface.Temperature.Lt.min,3) + poly(Present.Surface.Temperature.Mean,3),
    #           family=binomial(link='logit'),
    #           data=ENV_PA)
    #ea <- evaluate(presence_locations.pres_test,background_macaronesia_test,GAM,env)
    #tss$tssa <- c(tss$tssa,ea)
    #ta <- threshold(ea,stat='prevalence')
    #for (i in 1:length(envs)){
    #  predictions$predictionsa[[i]][[1]] <- c(predictions$predictionsa[[i]][[1]],predict(envs[[i]], GAM, ext=macaronesia.extent,type="response"))
    #  predictions$predictionsa[[i]][[2]] <- c(predictions$predictionsa[[i]][[2]],predict(envs[[i]], GAM, ext=macaronesia.extent2,type="response"))
    #}


    
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
    tss$tssmX <- c(tss$tssmX,em)
    # threshes$mX <- c(threshes$mX,threshold(em,stat='prevalence'))
    for (i in 1:length(envs)){
      predictions$predictionsmX[[i]][[1]] <- c(predictions$predictionsmX[[i]][[1]],predict(envs[[i]], maxEnt, ext=macaronesia.extent,type="response"))
      predictions$predictionsmX[[i]][[2]] <- c(predictions$predictionsmX[[i]][[2]],predict(envs[[i]], maxEnt, ext=macaronesia.extent2,type="response"))
    }

    
    
    glmm <- glm(presence ~ Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + Present.Benthic.Mean.Depth.Salinity.Lt.max + Present.Benthic.Mean.Depth.Salinity.Lt.min + Present.Benthic.Mean.Depth.Salinity.Mean + Present.Benthic.Mean.Depth.Temperature.Lt.max + Present.Benthic.Mean.Depth.Temperature.Lt.min + Present.Benthic.Mean.Depth.Temperature.Mean + Present.Surface.Current.Velocity.Mean.asc.BOv2_1 + Present.Surface.Salinity.Lt.max + Present.Surface.Salinity.Lt.min + Present.Surface.Salinity.Mean + Present.Surface.Temperature.Lt.max + Present.Surface.Temperature.Lt.min + Present.Surface.Temperature.Mean,
                family=binomial(link='logit'),
                data=ENV_PA)
    elm <- evaluate(presence_locations.pres_test,background_macaronesia_test,glmm,env)
    tss$tssml <- c(tss$tssml,elm)
    # threshes$ml <- c(threshes$ml,threshold(elm,stat='prevalence'))
    for (i in 1:length(envs)){
      predictions$predictionsml[[i]][[1]] <- c(predictions$predictionsml[[i]][[1]],predict(envs[[i]], glmm, ext=macaronesia.extent,type="response"))
      predictions$predictionsml[[i]][[2]] <- c(predictions$predictionsml[[i]][[2]],predict(envs[[i]], glmm, ext=macaronesia.extent2,type="response"))
    }

    
    
    # maha <- mahal(sp.env)
    # ema <- evaluate(presence_locations.pres_test,background_macaronesia_test,maha,env)
    # tssma <- c(tssma,ema)
    # tma <-  threshold(ema,stat='prevalence')
    # pma <- predict(env, maha, ext=macaronesia.extent,type="response")>tma
    # predictionsma1 <- c(predictionsma1, pma)
    # pma <- predict(env, maha, ext=macaronesia.extent2,type="response")>tma
    # predictionsma2 <- c(predictionsma2, pma)

    
    
    #d <- domain(env,presence_locations.pres_train)
    #d <- domain(sp.env)
    #ed <- evaluate(presence_locations.pres_test,background_macaronesia_test,d,env)
    #tss$tssd <- c(tss$tssd,ed)
    #td <- threshold(ed,stat='prevalence')
    #for (i in 1:length(envs)){
    #  predictions$predictionsd[[i]][[1]] <- c(predictions$predictionsd[[i]][[1]],predict(envs[[i]], d, ext=macaronesia.extent,type="response"))
    #  predictions$predictionsd[[i]][[2]] <- c(predictions$predictionsd[[i]][[2]],predict(envs[[i]], d, ext=macaronesia.extent2,type="response"))
    #}


    
    bio <- bioclim(env,presence_locations.pres_train)
    bio <- bioclim(sp.env)
    eb <- evaluate(presence_locations.pres_test,background_macaronesia_test,bio,env)
    tss$tssb <- c(tss$tssb,eb)
    # threshes$b <- c(threshes$b,threshold(eb,stat='prevalence'))
    for (i in 1:length(envs)){
      predictions$predictionsb[[i]][[1]] <- c(predictions$predictionsb[[i]][[1]],predict(envs[[i]], bio, ext=macaronesia.extent,type="response"))
      predictions$predictionsb[[i]][[2]] <- c(predictions$predictionsb[[i]][[2]],predict(envs[[i]], bio, ext=macaronesia.extent2,type="response"))
    }
  }
  
  #cat('BRT model\n')
  #wB <- sapply(tss$tssB, function(x) TSSweight(x))
  #tss$tssB <- weighted.mean(sapply(tss$tssB, function(x) TSS(x)), wB, na.rm=TRUE)
  #add_title <- 'BRT'
  #for (i in 1:length(env_names)){
  #  pB1 <- stack(predictions$predictionsB[[i]][[1]])
  #  pB1 <- weighted.mean(pB1, wB, na.rm=TRUE)
  #  pB2 <- stack(predictions$predictionsB[[i]][[2]])
  #  pB2 <- weighted.mean(pB2, wB, na.rm=TRUE)
  #  pB1 <- do.call(merge,list(pB1,pB2))
  #  mods <- pB1
  #  m <-  plotting(mods,species_name,tss$tssB,add_title,presence_locations)
  #  saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
  #  predictions$predictionsB[[i]] <- pB1
  #}


  
  # cat("Neural Network model\n")
  # pN <- stack(predictionsnn)
  # wnn <- sapply(tssnn, function(x) TSSweight(x))
  # tssnn <- weighted.mean(sapply(tssnn, function(x) TSS(x)), wnn, na.rm=TRUE)
  # pN <- weighted.mean(pN, wnn, na.rm=TRUE)
  # mods <- pN
  # add_title <- 'Neural Network'
  # m <-  plotting(mods,species_name,tssnn,add_title,presence_locations)
  # saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))
  
  
  #cat('GAM model\n')
  #wa <- sapply(tss$tssa, function(x) TSSweight(x))
  #tss$tssa <- weighted.mean(sapply(tss$tssa, function(x) TSS(x)), wa, na.rm=TRUE)
  #add_title <- 'GAM'
  #for (i in 1:length(env_names)){
  #  pB1 <- stack(predictions$predictionsa[[i]][[1]])
  #  pB1 <- weighted.mean(pB1, wa, na.rm=TRUE)
  #  pB2 <- stack(predictions$predictionsa[[i]][[2]])
  #  pB2 <- weighted.mean(pB2, wa, na.rm=TRUE)
  #  pB1 <- do.call(merge,list(pB1,pB2))
  #  mods <- pB1
  #  m <-  plotting(mods,species_name,tss$tssa,add_title,presence_locations)
  #  saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
  #  predictions$predictionsa[[i]] <- pB1
  #}


  
  cat('GLM model\n')
  wml <- sapply(tss$tssml, function(x) TSSweight(x))
  tss$tssml <- weighted.mean(sapply(tss$tssml, function(x) TSS(x)), wml, na.rm=TRUE)
  #threshes$ml <- weighted.mean(threshes$ml, wml, na.rm=TRUE)
  add_title <- 'GLM'
  for (i in 1:length(env_names)){
    pB1 <- stack(predictions$predictionsml[[i]][[1]])
    pB1 <- weighted.mean(pB1, wml, na.rm=TRUE)
    pB2 <- stack(predictions$predictionsml[[i]][[2]])
    pB2 <- weighted.mean(pB2, wml, na.rm=TRUE)
    pB1 <- do.call(merge,list(pB1,pB2))
    mods <- pB1
    m <-  plotting(mods,species_name,tss$tssml,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
    predictions$predictionsml[[i]] <- pB1
  }


  
  # cat('Mahalanobis model\n')
  # pma1 <- stack(predictionsma1)
  # wma <- sapply(tssma, function(x) TSSweight(x))
  # #tssma <- sapply(tssma, function(x) TSS(x))
  # tssma <- weighted.mean(sapply(tssma, function(x) TSS(x)), wma, na.rm=TRUE)
  # pma1 <- weighted.mean(pma1, wma, na.rm=TRUE)
  # pma1[pma1<0] <- 0
  # pma2 <- stack(predictionsma2)
  # pma2 <- weighted.mean(pma2, wma, na.rm=TRUE)
  # pma1 <- do.call(merge,list(pma1, pma2))
  # rm(pma2)
  # #pma[pma==NA] <- 0
  # mods <- pma1
  # add_title <- 'Mahalanobis'
  # m <-  plotting(mods,species_name,tssma,add_title,presence_locations)
  # saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  # 

  
  cat('Bioclim model\n')#filter of .2 TSS is removing the layers
  wb <- sapply(tss$tssb, function(x) TSSweight(x))
  tss$tssb <- weighted.mean(sapply(tss$tssb, function(x) TSS(x)), wb, na.rm=TRUE)
  #threshes$b <- weighted.mean(threshes$b, wb, na.rm=TRUE)
  add_title <- 'Bioclim'
  #tssb <- sapply(tssb, function(x) TSS(x))
  for (i in 1:length(env_names)){
    pB1 <- stack(predictions$predictionsb[[i]][[1]])
    pB1 <- weighted.mean(pB1, wb, na.rm=TRUE)
    pB2 <- stack(predictions$predictionsb[[i]][[2]])
    pB2 <- weighted.mean(pB2, wb, na.rm=TRUE)
    pB1 <- do.call(merge,list(pB1,pB2))
    mods <- pB1
    m <-  plotting(mods,species_name,tss$tssb,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
    predictions$predictionsb[[i]] <- pB1
  }


  
  #cat('Domain model\n')
  #wd <- sapply(tss$tssd, function(x) TSSweight(x))
  #tss$tssd <- weighted.mean(sapply(tss$tssd, function(x) TSS(x)), wd, na.rm=TRUE)
  #add_title <- 'Domain'
  #for (i in 1:length(env_names)){
  #  pB1 <- stack(predictions$predictionsd[[i]][[1]])
  #  pB1 <- weighted.mean(pB1, wd, na.rm=TRUE)
  #  pB2 <- stack(predictions$predictionsd[[i]][[2]])
  #  pB2 <- weighted.mean(pB2, wd, na.rm=TRUE)
  #  pB1 <- do.call(merge,list(pB1,pB2))
  #  mods <- pB1
  #  m <-  plotting(mods,species_name,tss$tssd,add_title,presence_locations)
  #  saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
  #  predictions$predictionsd[[i]] <- pB1
  #}


  
  cat('MaxEnt model\n')
  wmX <- sapply(tss$tssmX, function(x) TSSweight(x))
  #tssmx <- sapply(tssmx, function(x) TSS(x))
  tss$tssmX <- weighted.mean(sapply(tss$tssmX, function(x) TSS(x)), wmX, na.rm=TRUE)
  #threshes$mX <- weighted.mean(threshes$mX, wmX, na.rm=TRUE)
  add_title <- 'MaxEnt'
  for (i in 1:length(env_names)){
    pB1 <- stack(predictions$predictionsmX[[i]][[1]])
    pB1 <- weighted.mean(pB1, wmX, na.rm=TRUE)
    pB2 <- stack(predictions$predictionsmX[[i]][[2]])
    pB2 <- weighted.mean(pB2, wmX, na.rm=TRUE)
    pB1 <- do.call(merge,list(pB1,pB2))
    mods <- pB1
    m <-  plotting(mods,species_name,tss$tssmX,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
    predictions$predictionsmX[[i]] <- pB1
  }

  
  
  cat('Ensemble\n')
  tsss <- c(tss$tssml, tss$tssb, tss$tssmX)#tssB, tssnn, tssma,tss$tssa, tss$tssd, 
  ws <- weights_tsss(tsss)
  tsss <- weighted.mean(tsss,ws, na.rm=TRUE)
  add_title <- 'Ensemble'
  for (i in 1:length(env_names)){
    models <- stack(predictions$predictionsml[[i]],predictions$predictionsb[[i]],predictions$predictionsmX[[i]])
    names(models) <- c("GLM", "bioclim", "maxent")
    Mean <- weighted.mean( models[[c("GLM", "bioclim", "maxent")]], ws, na.rm=TRUE)
    mods <- Mean
    m <-  plotting(mods,species_name,tsss,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
    predictions$predictionsb[[i]] <- Mean
    # predictions$predictionsml[[i]] <- 0#
    # predictions$predictionsmX[[i]] <- 0#
  }
  n <- c('_RCP26_','_RCP45_','_RCP60_','_RCP85_')
  cat('Different RCP scenarios\n')
  for (i in 1:4){
    rasterR <- predictions$predictionsb[[i+5]] - predictions$predictionsb[[1]]
    rasterG <- predictions$predictionsb[[i+1]] - predictions$predictionsb[[1]]
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
  predictions <- list(
    predictionsB = list(
      predictions = c(),
      predictions1 = c(),
      predictions2 = c(),
      predictions3 = c(),
      predictions4 = c(),
      predictions5 = c(),
      predictions6 = c(),
      predictions7 = c(),
      predictions8 = c()),
    predictionsa = list(
      predictions = c(),
      predictions1 = c(),
      predictions2 = c(),
      predictions3 = c(),
      predictions4 = c(),
      predictions5 = c(),
      predictions6 = c(),
      predictions7 = c(),
      predictions8 = c()),
    predictionsmX = list(
      predictions = c(),
      predictions1 = c(),
      predictions2 = c(),
      predictions3 = c(),
      predictions4 = c(),
      predictions5 = c(),
      predictions6 = c(),
      predictions7 = c(),
      predictions8 = c()),
    predictionsml = list(
      predictions = c(),
      predictions1 = c(),
      predictions2 = c(),
      predictions3 = c(),
      predictions4 = c(),
      predictions5 = c(),
      predictions6 = c(),
      predictions7 = c(),
      predictions8 = c()),
    predictionsd = list(
      predictions = c(),
      predictions1 = c(),
      predictions2 = c(),
      predictions3 = c(),
      predictions4 = c(),
      predictions5 = c(),
      predictions6 = c(),
      predictions7 = c(),
      predictions8 = c()),
    predictionsb = list(
      predictions = c(),
      predictions1 = c(),
      predictions2 = c(),
      predictions3 = c(),
      predictions4 = c(),
      predictions5 = c(),
      predictions6 = c(),
      predictions7 = c(),
      predictions8 = c())
    )
  tss <- list(
    tssB <- c(),
    tssa <- c(),
    tssmX <- c(),
    tssml <- c(),
    tssd <- c(),
    tssb <- c()
    )
  threshes <- list(
    B <- c(),
    a <- c(),
    mX <- c(),
    ml <- c(),
    d <- c(),
    b <- c()
    )#thresholds
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
    # BRT <- gbm(presence ~ Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + Present.Benthic.Mean.Depth.Salinity.Lt.max + Present.Benthic.Mean.Depth.Salinity.Lt.min + Present.Benthic.Mean.Depth.Salinity.Mean + Present.Benthic.Mean.Depth.Temperature.Lt.max + Present.Benthic.Mean.Depth.Temperature.Lt.min + Present.Benthic.Mean.Depth.Temperature.Mean + Present.Surface.Current.Velocity.Mean.asc.BOv2_1 + Present.Surface.Salinity.Lt.max + Present.Surface.Salinity.Lt.min + Present.Surface.Salinity.Mean + Present.Surface.Temperature.Lt.max + Present.Surface.Temperature.Lt.min + Present.Surface.Temperature.Mean,
    # #BRT <- gbm(presence ~ X2100AOGCM.RCP85.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Current.Velocity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Salinity.Mean.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1 + X2100AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1,
    #            cv.folds = 100,
    #            interaction.depth = 3,
    #            distribution = "bernoulli",
    #            shrinkage = 0.01,
    #            verbose = FALSE,
    #            n.trees = 500,
    #            data=ENV_PA)
    # eB <- evaluate(presence_locations.pres_test,background_macaronesia_test,BRT,env)
    # tssB <- c(tssB,eB)
    # tB <- threshold(eB,stat='prevalence')
    # pB <- predict(env, BRT, ext=macaronesia.extent,type="response")>tB
    # predictionsB <- c(predictionsB,pB)

    
    
    # #    GAM <- gam(presence ~ poly(BO_sstmean,3)+poly(BO_dissox,3)+poly(BO2_tempmean_bdmin,3)+poly(BO2_salinitymean_bdmin,3)+poly(BO2_salinitymean_ss,3)+poly(BO2_carbonphytomean_ss,3)+poly(BO21_carbonphytomean_bdmin,3),
    # GAM <- gam(presence ~ poly(Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1,3) + poly(Present.Benthic.Mean.Depth.Salinity.Lt.max,3) + poly(Present.Benthic.Mean.Depth.Salinity.Lt.min,3) + poly(Present.Benthic.Mean.Depth.Salinity.Mean,3) + poly(Present.Benthic.Mean.Depth.Temperature.Lt.max,3) + poly(Present.Benthic.Mean.Depth.Temperature.Lt.min,3) + poly(Present.Benthic.Mean.Depth.Temperature.Mean,3) + poly(Present.Surface.Current.Velocity.Mean.asc.BOv2_1,3) + poly(Present.Surface.Salinity.Lt.max,3) + poly(Present.Surface.Salinity.Lt.min,3) + poly(Present.Surface.Salinity.Mean,3) + poly(Present.Surface.Temperature.Lt.max,3) + poly(Present.Surface.Temperature.Lt.min,3) + poly(Present.Surface.Temperature.Mean,3),
    # #GAM <- gam(presence ~ poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Salinity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Benthic.Mean.Depth.Temperature.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Current.Velocity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Salinity.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Salinity.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Salinity.Mean.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Temperature.Lt.max.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Temperature.Lt.min.asc.BOv2_1, 3) + poly(X2100AOGCM.RCP85.Surface.Temperature.Mean.asc.BOv2_1, 3),
    #            family=binomial(link='logit'),
    #            data=ENV_PA)
    # ea <- evaluate(presence_locations.pres_test,background_macaronesia_test,GAM,env)
    # tssa <- c(tssa,ea)
    # ta <- threshold(ea,stat='prevalence')
    # pa <- predict(env, GAM, ext=macaronesia.extent,type="response")>ta
    # predictionsa <- c(predictionsa,pa)
    
    
    
    # #cat('neural')
    # NN <- nnet(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
    #            data=ENV_PA,
    #            size=10)
    # enn <- evaluate(presence_locations.pres_test,background_macaronesia_test,NN,env)
    # predictionsnn <- c(predictionsnn,predict(env, NN, ext=macaronesia.extent,type="response"))
    # tssnn <- c(tssnn,enn)
    
    

    maxEnt <- maxent(x=presence_absence_env,presence_absence,silent=TRUE)
    em <- evaluate(presence_locations.pres_test,background_macaronesia_test,maxEnt,env)
    tss$tssmX <- c(tss$tssmX,em)
    # threshes$mX <- c(threshes$mX,threshold(em,stat='prevalence'))
    for (i in 1:length(envs)){
      predictions$predictionsmX[[i]] <- c(predictions$predictionsmX[[i]],predict(envs[[i]], maxEnt, ext=macaronesia.extent,type="response"))
    }
    
    

    glmm <- glm(presence ~ Present.Benthic.Mean.Depth.Current.Velocity.Mean.asc.BOv2_1 + Present.Benthic.Mean.Depth.Salinity.Lt.max + Present.Benthic.Mean.Depth.Salinity.Lt.min + Present.Benthic.Mean.Depth.Salinity.Mean + Present.Benthic.Mean.Depth.Temperature.Lt.max + Present.Benthic.Mean.Depth.Temperature.Lt.min + Present.Benthic.Mean.Depth.Temperature.Mean + Present.Surface.Current.Velocity.Mean.asc.BOv2_1 + Present.Surface.Salinity.Lt.max + Present.Surface.Salinity.Lt.min + Present.Surface.Salinity.Mean + Present.Surface.Temperature.Lt.max + Present.Surface.Temperature.Lt.min + Present.Surface.Temperature.Mean,
                family=binomial(link='logit'),
                data=ENV_PA)
    elm <- evaluate(presence_locations.pres_test,background_macaronesia_test,glmm,env)
    tss$tssml <- c(tss$tssml,elm)
    # threshes$ml <- c(threshes$ml,threshold(elm,stat='prevalence'))
    for (i in 1:length(envs)){
      predictions$predictionsml[[i]] <- c(predictions$predictionsml[[i]],predict(envs[[i]], glmm, ext=macaronesia.extent,type="response"))
    }
    
    
    
    # maha <- mahal(sp.env)
    # ema <- evaluate(presence_locations.pres_test,background_macaronesia_test,maha,env)
    # tssma <- c(tssma,ema)
    # tma <-  threshold(ema,stat='prevalence')
    # pma <- predict(env, maha, ext=macaronesia.extent,type="response")>tma
    # predictionsma <- c(predictionsma, pma)
    
    
    
    # #d <- domain(env,presence_locations.pres_train)
    # d <- domain(sp.env)
    # ed <- evaluate(presence_locations.pres_test,background_macaronesia_test,d,env)
    # tssd <- c(tssd,ed)
    # td <- threshold(ed,stat='prevalence')
    # pd <- predict(env, d, ext=macaronesia.extent,type="response")>td
    # predictionsd <- c(predictionsd,pd)
    

    
    bio <- bioclim(sp.env)
    eb <- evaluate(presence_locations.pres_test,background_macaronesia_test,bio,env)
    tss$tssb <- c(tss$tssb,eb)
    # threshes$b <- c(threshes$b,threshold(eb,stat='prevalence'))
    for (i in 1:length(envs)){
      predictions$predictionsb[[i]] <- c(predictions$predictionsb[[i]],predict(envs[[i]], bio, ext=macaronesia.extent,type="response"))
    }

  }
  
  # cat('BRT model\n')
  # wB <- sapply(tssB, function(x) TSSweight(x))
  # tssB <- weighted.mean(sapply(tssB, function(x) TSS(x)), wB, na.rm=TRUE)
  # pB <- stack(predictionsB)
  # pB <- weighted.mean(pB, wB, na.rm=TRUE)
  # mods <- pB
  # add_title <- 'BRT'
  # m <-  plotting(mods,species_name,tssB,add_title,presence_locations)
  # saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  

  
  # cat("Neural Network model\n")
  # pN <- stack(predictionsnn)
  # wnn <- sapply(tssnn, function(x) TSSweight(x))
  # tssnn <- weighted.mean(sapply(tssnn, function(x) TSS(x)), wnn, na.rm=TRUE)
  # pN <- weighted.mean(pN, wnn, na.rm=TRUE)
  # mods <- pN
  # add_title <- 'Neural Network'
  # m <-  plotting(mods,species_name,tssnn,add_title,presence_locations)
  # saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep=''))
  
  
  # cat('GAM model\n')
  # wa <- sapply(tssa, function(x) TSSweight(x))
  # tssa <- weighted.mean(sapply(tssa, function(x) TSS(x)), wa, na.rm=TRUE)
  # pa <- stack(predictionsa)
  # pa <- weighted.mean(pa, wa, na.rm=TRUE)
  # mods <- pa
  # add_title <- 'GAM'
  # m <-  plotting(mods,species_name,tssa,add_title,presence_locations)
  # saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))


  
  cat('GLM model\n')
  wml <- sapply(tss$tssml, function(x) TSSweight(x))
  tss$tssml <- weighted.mean(sapply(tss$tssml, function(x) TSS(x)), wml, na.rm=TRUE)
  add_title <- 'GLM'
  for (i in 1:length(env_names)){
    pB1 <- stack(predictions$predictionsml[[i]])
    pB1 <- weighted.mean(pB1, wml, na.rm=TRUE)
    mods <- pB1
    m <-  plotting(mods,species_name,tss$tssml,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
    predictions$predictionsml[[i]] <- pB1
  }


  
  # cat('Mahalanobis model\n')
  # pma <- stack(predictionsma)
  # wma <- sapply(tssma, function(x) TSSweight(x))
  # #tssma <- sapply(tssma, function(x) TSS(x))
  # tssma <- weighted.mean(sapply(tssma, function(x) TSS(x)), wma, na.rm=TRUE)
  # pma <- weighted.mean(pma, wma, na.rm=TRUE)
  # pma[pma<0] <- 0
  # #pma[pma==NA] <- 0
  # mods <- pma
  # add_title <- 'Mahalanobis'
  # m <-  plotting(mods,species_name,tssma,add_title,presence_locations)
  # saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  


  cat('Bioclim model\n')#filter of .2 TSS is removing the layers
  wb <- sapply(tss$tssb, function(x) TSSweight(x))
  tss$tssb <- weighted.mean(sapply(tss$tssb, function(x) TSS(x)), wb, na.rm=TRUE)
  add_title <- 'Bioclim'
  for (i in 1:length(env_names)){
    pB1 <- stack(predictions$predictionsb[[i]])
    pB1 <- weighted.mean(pB1, wb, na.rm=TRUE)
    mods <- pB1
    m <-  plotting(mods,species_name,tss$tssb,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
    predictions$predictionsb[[i]] <- pB1
  }
  

  
  # cat('Domain model\n')
  # pd <- stack(predictionsd)
  # wd <- sapply(tssd, function(x) TSSweight(x))
  # #tssd <- sapply(tssd, function(x) TSS(x))
  # tssd <- weighted.mean(sapply(tssd, function(x) TSS(x)), wd, na.rm=TRUE)
  # pd <- weighted.mean(pd, wd, na.rm=TRUE)
  # mods <- pd
  # add_title <- 'Domain'
  # m <-  plotting(mods,species_name,tssd,add_title,presence_locations)
  # saveWidget(m, file=paste(species_name,'_Present_',add_title,'.html', sep=''))
  


  cat('MaxEnt model\n')
  wmX <- sapply(tss$tssmX, function(x) TSSweight(x))
  tss$tssmX <- weighted.mean(sapply(tss$tssmX, function(x) TSS(x)), wmX, na.rm=TRUE)
  add_title <- 'MaxEnt'
  for (i in 1:length(env_names)){
    pB1 <- stack(predictions$predictionsmX[[i]])
    pB1 <- weighted.mean(pB1, wmX, na.rm=TRUE)
    mods <- pB1
    m <-  plotting(mods,species_name,tss$tssmX,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
    predictions$predictionsmX[[i]] <- pB1
  }
  


  cat('Ensemble\n')
  tsss <- c(tss$tssml, tss$tssb, tss$tssmX)
  ws <- weights_tsss(tsss)
  tsss <- weighted.mean(tsss,ws, na.rm=TRUE)
  add_title <- 'Ensemble'
  for (i in 1:length(env_names)){
    models <- stack(predictions$predictionsml[[i]],predictions$predictionsb[[i]],predictions$predictionsmX[[i]])
    names(models) <- c("GLM", "bioclim", "maxent")
    Mean <- weighted.mean( models[[c("GLM", "bioclim", "maxent")]], ws, na.rm=TRUE)
    mods <- Mean
    m <-  plotting(mods,species_name,tsss,add_title,presence_locations)
    saveWidget(m, file=paste(species_name,env_names[[i]],add_title,'.html', sep=''))
    predictions$predictionsb[[i]] <- Mean
  }
  cat('Different RCP scenarios\n')
  n <- c('_RCP26_','_RCP45_','_RCP60_','_RCP85_')
  for (i in 1:4){
    rasterR <- predictions$predictionsb[[i+5]] - predictions$predictionsb[[1]]
    rasterG <- predictions$predictionsb[[i+1]] - predictions$predictionsb[[1]]
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
