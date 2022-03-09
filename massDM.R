library(dismo)
library(ncdf4)
library(dplyr)
library(raster)
library(leaflet)
library(htmlwidgets)
library(sdmpredictors)
#library(ENMTools)
#library(nnet)#neural networks
library(gbm)#BRT
library(mgcv)#GAM GARP?
library(quantreg)


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

# species_locations <- function(locations,species_name){
#   cranch_locations <- select(locations[locations$species==species_name,],
#                              longitude,
#                              latitude)
#   return(cranch_locations)
# }

species_locations <- function(locations,species_name){
  cranch_locations <- select(locations[sub(" .*", "", oegopsids$species)==species_name,],
                             longitude,
                             latitude)
  return(cranch_locations)
}

family_locations <- function(locations,species_name){
  cranch_locations <- select(locations[oegopsids$family==species_name,],
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

env <- load_layers(list_layers( datasets="Bio-ORACLE" )[c(21,12,158,212,310,298,335),]$layer_code)
#25,21,12,69,119,158,212,310,298,335
#bathi, SST,O2, "min deep O2",light bottom, temp bottom, sal bottom, sal surface,biomass,biomass deep
mask_depth <- load_layers(list_layers( datasets="Bio-ORACLE" )[25,]$layer_code)
mask_depth <- mask_depth[[1]]
#mask_depth <- crop(mask_depth,macaronesia.extent)
mask_depth[mask_depth>-200] <- NA
env <- mask(env,mask_depth)
#macaronesia.extent <- extent(-35, -10, 12, 42)#lon,lon,lat,lat
macaronesia.extent <- extent(-55, 0, -14, 50)
macaronesia.env <- crop(env,macaronesia.extent)




oegopsids <- read.csv("occorrences.csv")
oegopsid_locations <- select(oegopsids,
                             longitude,
                             latitude)
biomasses <- read.csv("biomasses.csv")


background_macaronesia <- randomPoints(env, n=1000, ext=macaronesia.extent, extf = 1.25)

colnames(background_macaronesia) = c('longitude', 'latitude')
group <- kfold(background_macaronesia, 4)
background_macaronesia_train <- background_macaronesia[group != 1, ]
background_macaronesia_train.env <- extract(env, background_macaronesia_train)
background_macaronesia_test <- background_macaronesia[group == 1, ]

species_name <- 'Histioteuthidae'#'Histioteuthis reversa','Liocranchia reinhardtii','Cranchia scabra','Pterygioteuthis gemmata'

cat(paste(species_name),'-',dim(family_locations(oegopsids,species_name))[1])
#presence_locations <- species_locations(oegopsids,species_name)
presence_locations <- family_locations(oegopsids,species_name)
presence_locations <- presence_locations[presence_locations$latitude>-20,]
presence_locations <- presence_locations[presence_locations$longitude>-60,]
presence_locations <- presence_locations[presence_locations$longitude<10,]


j <- na.omit(extract(env,presence_locations))
presence_locations <- presence_locations[-c(attr(x=j, which='na.action')),]
presence_locations <- spThin(presence_locations,dist=500)[[1]]
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

for (i in 1:2){
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
tssB <- c(tssB,eB)


GAM <- gam(presence ~ poly(BO_sstmean,3)+poly(BO_dissox,3)+poly(BO2_tempmean_bdmin,3)+poly(BO2_salinitymean_bdmin,3)+poly(BO2_salinitymean_ss,3)+poly(BO2_carbonphytomean_ss,3)+poly(BO21_carbonphytomean_bdmin,3),
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


glmm <- glm(presence ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
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
#predictionsma[predictionsma<0] <- 0
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
m

#net_locations <- select(biomasses[sub(" .*", "", biomasses$species)==species_name,],longitude,latitude)
#biomass <- select(biomasses[sub(" .*", "", biomasses$species)==species_name,],BIOMASS)
net_locations <- select(biomasses[biomasses$species==species_name,],longitude,latitude)
biomass <- select(biomasses[biomasses$species==species_name,],BIOMASS)
lineb <- extract(Mean,net_locations)
lineb <- cbind(biomass, lineb)
lineb <- unique(na.omit(lineb))
plot(lineb$lineb,lineb$BIOMASS)
#probs <- extract(env,net_locations)
#probs <- na.omit(cbind(probs,biomass))
#probs <- unique(probs)
#LM <- glm(BIOMASS ~ BO_sstmean+BO_dissox+BO2_tempmean_bdmin+BO2_salinitymean_bdmin+BO2_salinitymean_ss+BO2_carbonphytomean_ss+BO21_carbonphytomean_bdmin,
#               data=probs)
#BIOMASS <- predict(env, LM, ext=macaronesia.extent,type="response")>0
LM <- rq(BIOMASS ~ lineb, data=lineb, tau = seq(0, 1, by = 0.05))
LM <- lm(BIOMASS ~ lineb, data=lineb)
summary(LM)
coef <- summary(LM)$coefficients[,1]
##############################################


BIOMASS <- Mean*coef[2] + coef[1]
mods <- BIOMASS#*Mean
add_title <- 'Biomass'
palete = colorNumeric(c("#5E85B8","#EDF0C0","#C13127"), values(BIOMASS), na.color = 'transparent')
name <- species_name
m <-  leaflet()%>%
   addProviderTiles(provider = 'Esri.OceanBasemap') %>%
   addScaleBar()%>%
   addRasterImage(colors = palete, mods, opacity = .6)%>%
   addCircles(net_locations$longitude,net_locations$latitude)%>%
   addLegend(position = 'topright',pal = palete, values(mods), title = paste(name,' - ',add_title,', TSS: ',tsss))
m
saveWidget(m, file=paste(species_name,'_',add_title,'.html', sep='')) 
