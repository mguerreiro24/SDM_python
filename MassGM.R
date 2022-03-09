library(dismo)
library(ncdf4)
library(dplyr)
library(raster)
library(leaflet)
library(htmlwidgets)
library(sdmpredictors)
library(gbm)#BRT
library(mgcv)#GAM GARP?

setwd("D:/PhD/articles/article SDM/")


'--------------------------------Functions----------------------------------------'


# spThin <- function(xy,dist,rep=1) {
#   #dist -> distance between points in kilometers
#   #rep -> replicates
#   .r <- raster(resolution=0.008333333) # empty raster with 1km resolution
#   .a <- raster(ext=extent(.r),resolution=res(.r)[1]*dist)
#   .ac <- cellFromXY(.a,xy)
#   .tbl <- data_frame(cell=.ac,x=xy[,1],y=xy[,2]) %>% group_by(cell)
#   o <- list()
#   for (i in 1:rep) {
#     .s <- .tbl %>% sample_n(1) %>% ungroup()
#     o[[i]] <- data.frame(.s %>% dplyr::select(x,y))
#   }
#   o
# }

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

# species_locations <- function(locations,species_name){
#    cranch_locations <- select(locations[locations$species==species_name,],
#                               longitude,
#                               latitude)
#    return(cranch_locations)
#  }
# 
# genus_locations <- function(locations,species_name){
#   cranch_locations <- select(locations[sub(" .*", "", oegopsids$species)==species_name,],
#                              longitude,
#                              latitude)
#   return(cranch_locations)
# }
# 
# family_locations <- function(locations,species_name){
#   cranch_locations <- select(locations[oegopsids$family==species_name,],
#                              longitude,
#                              latitude)
#   return(cranch_locations)
# }

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
    addCircles(points$longitude,points$latitude)%>%
    addLegend(position = 'topright',pal = palete, values(mods), title = paste(name,' - ',add_title,', TSS: ',TSS))
  return(m)
}


'--------------------------------Script----------------------------------------'


taxon_name <- 'Cranchiidae'

env <- load_layers(list_layers( datasets="Bio-ORACLE" )[c(21,12,158,212,310,298,335),]$layer_code)
#25,21,12,69,119,158,212,310,298,335
#bathi, SST,O2, "min deep O2",light bottom, temp bottom, sal bottom, sal surface,biomass,biomass deep

mask_depth <- load_layers(list_layers( datasets="Bio-ORACLE" )[25,]$layer_code)
mask_depth <- mask_depth[[1]]
mask_depth[mask_depth>-200] <- NA
env <- mask(env,mask_depth)

#macaronesia.extent <- extent(-35, -10, 12, 42)#lon,lon,lat,lat
macaronesia.extent <- extent(-55, 0, -14, 50)
macaronesia.env <- crop(env,macaronesia.extent)


biomasses <- read.csv("biomasses.csv")
net_locations <- select(biomasses[biomasses$species==taxon_name,],longitude,latitude,BIOMASS)
background_macaronesia <- select(net_locations[net_locations$BIOMASS==0.0,],longitude,latitude)
group <- kfold(background_macaronesia, 4)
background_macaronesia_train <- background_macaronesia[group != 1, ]
background_macaronesia_train.env <- extract(env, background_macaronesia_train)
background_macaronesia_test <- background_macaronesia[group == 1, ]

presence_locations <- select(net_locations[net_locations$BIOMASS>0.0,],longitude,latitude)


#predictionsB <- c()
#tssB <- c()
predictionsa <- c()
tssa <- c()
#predictionsml <- c()
#tssml <- c()

for (i in 1:200){
group <- kfold(presence_locations, 4)
presence_locations.pres_train <- presence_locations[group != 1, ]
presence_locations.pres_test <- presence_locations[group == 1, ]
sp.env <- extract(env,presence_locations.pres_train)
presence_absence_env <- as.data.frame(rbind(sp.env,background_macaronesia_train.env))
presence_absence <- c(select(net_locations[net_locations$BIOMASS>0.0,], BIOMASS)[group != 1, ]/max(select(net_locations[net_locations$BIOMASS>0.0,], BIOMASS)), rep(0,dim(background_macaronesia_train.env)[1]))
ENV_PA <- presence_absence_env
ENV_PA$presence <- presence_absence




cat(paste(i,'\n'))
# BRT <- gbm(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
# cv.folds = 100,
# interaction.depth = 3,
# distribution = "bernoulli",
# shrinkage = 0.01,
# verbose = FALSE,
# n.trees = 500,
# data=ENV_PA)
# eB <- evaluate(presence_locations.pres_test,background_macaronesia_test,BRT,env)
# tB <- threshold(eB,stat='prevalence')
# pB <- predict(env, BRT, ext=macaronesia.extent,type="response")>tB 
# predictionsB <- c(predictionsB,pB)
# tssB <- c(tssB,eB)


GAM <- gam(presence ~ poly(BO_sstmean,3)+poly(BO_dissox,3)+poly(BO2_tempmean_bdmin,3)+poly(BO2_salinitymean_bdmin,3)+poly(BO2_salinitymean_ss,3)+poly(BO2_carbonphytomean_ss,3),
            family=binomial(link='logit'),
            data=ENV_PA)
ea <- evaluate(presence_locations.pres_test,background_macaronesia_test,GAM,env)
ta <- threshold(ea,stat='prevalence')
pa <- predict(env, GAM, ext=macaronesia.extent,type="response")>ta
predictionsa <- c(predictionsa,pa)
tssa <- c(tssa,ea)


# glmm <- glm(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
#             family=binomial(link='logit'),
#             data=ENV_PA)
# elm <- evaluate(presence_locations.pres_test,background_macaronesia_test,glmm,env)
# tlm <- threshold(elm,stat='prevalence')
# plm <- predict(env, glmm, ext=macaronesia.extent,type="response")>tlm
# predictionsml <- c(predictionsml,plm)
# tssml <- c(tssml,elm)
}


#cat('BRT model\n')
#pB <- stack(predictionsB)
#wB <- sapply(tssB, function(x) TSSweight(x))
#tssB <- weighted.mean(sapply(tssB, function(x) TSS(x)), wB, na.rm=TRUE)
#pB <- weighted.mean(pB, wB, na.rm=TRUE)
#mods <- pB
#add_title <- 'BRT'
#m <-  plotting(mods,taxon_name,tssB,add_title,presence_locations)
#saveWidget(m, file=paste(taxon_name,'_',add_title,'.html', sep=''))

cat('GAM model\n')
pa <- stack(predictionsa)
wa <- sapply(tssa, function(x) TSSweight(x))
tssa <- weighted.mean(sapply(tssa, function(x) TSS(x)), wa, na.rm=TRUE)
pa <- weighted.mean(pa, wa, na.rm=TRUE)
mods <- pa
add_title <- 'GAM'
m <-  plotting(mods,taxon_name,tssa,add_title,presence_locations)
saveWidget(m, file=paste(taxon_name,'_',add_title,'.html', sep=''))

# cat('GLM model\n')
# pml <- stack(predictionsml)
# wml <- sapply(tssml, function(x) TSSweight(x))
# tssml <- weighted.mean(sapply(tssml, function(x) TSS(x)), wml, na.rm=TRUE)
# pml <- weighted.mean(pml, wml, na.rm=TRUE)
# mods <- pml
# add_title <- 'GLM'
# m <-  plotting(mods,taxon_name,tssml,add_title,presence_locations)
# saveWidget(m, file=paste(taxon_name,'_',add_title,'.html', sep=''))


# tsss <- c(tssa, tssml)
# cat('Ensemble\n')
# models <- stack(pa, pml)
# names(models) <- c("GAM", "GLM")
# ws <- weights_tsss(tsss)
# Mean <- weighted.mean( models[[c("GAM", "GLM")]], ws, na.rm=TRUE)
# tsss <- weighted.mean(tsss,ws, na.rm=TRUE)
# mods <- Mean
mods <- mods*max(select(net_locations[net_locations$BIOMASS>0.0,], BIOMASS))
add_title <- 'Ensemble'

palete = colorNumeric(c("#5E85B8","#EDF0C0","#C13127"), c(0,max(select(net_locations[net_locations$BIOMASS>0.0,], BIOMASS))), na.color = 'transparent')
m <-  leaflet()%>%
  addProviderTiles(provider = 'Esri.OceanBasemap') %>%
  addScaleBar()%>%
  addRasterImage(colors = palete, mods, opacity = .6)%>%
  addCircles(presence_locations$longitude,presence_locations$latitude)%>%
  addCircles(background_macaronesia$longitude,background_macaronesia$latitude, color = "#FF0000")%>%
  addLegend(position = 'topright',pal = palete, values(mods), title = paste(taxon_name,' - ',add_title,', TSS: ',tssa))

saveWidget(m, file=paste(taxon_name,'_',add_title,'.html', sep=''))
m
