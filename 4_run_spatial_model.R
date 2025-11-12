##############################################################
## Spatial model -> estimating FOI for whole of Timor Leste ##
##############################################################

# key: FOI = force of infection (also called lambda)

library(dplyr)
library(tidyr)
library(ggplot2)
library(INLAutils)
library(INLA)
library(readxl)
library(sf)
library(tibble)
library(scales)

############################################################################
## Sample the catalytic model FOI estimates and match to the climate data ##
############################################################################

## read in posteriors from catalytic model
lambda <- as.data.frame(readRDS("spatial_model_data/FOI_catalytic_model_posteriors.RDS"))
tot_posterior <- dim(lambda)[1]


## sample n=100 from the catalytic model posteriors
int <- sample.int(size = 100, n = tot_posterior)
lambda |> dplyr::slice(int) -> lambda_samples


## check samples are representative the mean and median of posterior
lambda |> 
  pivot_longer(1:102) |> # number of EAs with catalytic model estimates = 102
  group_by(name) |> reframe(median=quantile(value, probs=c(0.5)),
                           ciL=quantile(value, probs=c(0.025)),
                           ciU=quantile(value, probs=c(0.975))) -> posterior_medians
saveRDS(posterior_medians, "catalytic_model_posterior_medians.RDS")

lambda_samples |> 
  pivot_longer(1:102) |>
  group_by(name)|>reframe(median=quantile(value, probs=c(0.5)),
                           ciL=quantile(value, probs=c(0.025)),
                           ciU=quantile(value, probs=c(0.975))) -> sample_medians

## plot to check the samples are representative of the full posterior
posterior_medians |> 
  mutate(type="full") |> 
  rbind(sample_medians |> mutate(type="samples")) |>
  ggplot() +
  geom_point(aes(x=as.character(name), y=median, col=type), position="dodge") +
  geom_errorbar(aes(x=as.character(name), ymin=ciL, ymax=ciU, col=type), position="dodge") +
  labs(x="EA", y="FOI", col=NULL) +
  theme_classic() +
  theme(axis.text.x = element_blank())


## read in climate data and coordinates for all EAs in Timor-Leste
climate_data <- read.csv("spatial_model_data/climate_FOI.csv")
coords <- read.csv("spatial_model_data/EA_centroids.csv")

## scale climate variables
standard <- function(x){ (x-mean(x,na.rm=T))/ sd(x, na.rm=T)}
climate_data |> 
  dplyr::select(ea_codes, ele_eemean, pop_eemean, pre_eemean, tmax_eemean, tmin_eemean,
                               rh_tmax_eemean, rh_tmin_eemean, mir_eemean, evi_eemean) |>
  mutate(pop_eemean = log(pop_eemean)) |>
  mutate(across(c(2:10), standard)) |>
  cbind(coords) |> rename(EA.name=ea_codes) -> climate_data


## reformat lambda samples
lambda_samples |> 
  mutate(index=1:100) |> # n sample
  pivot_longer(1:102) |> # n EA
  pivot_wider(values_from = value, names_from = index) |> 
  rename(EA.name=name) |> mutate(EA.name=as.integer(EA.name)) -> lambda_samples


## match climate data to lambdas
climate_data |> left_join(lambda_samples) -> data_all


#########################################################
## set up the spatial model runs for the n=100 samples ##
#########################################################

## spde and mesh
mesh <- readRDS("spatial_model_data/saved_mesh.RDS") # the mesh is a triangulated grid connecting the EA centroids 
spde <- readRDS("spatial_model_data/saved_spde.RDS") # SPDE is made using the mesh (see setup_spatial_model_sensitivity.R)

## spatial field indices
field.indices = inla.spde.make.index("field", n.spde=mesh$n)

## select and remove EAs to save for cross validation
prep_cv <- function(n){
  pred.ea <- sample.int(size=10, n=102) 
  n[pred.ea] <-NA
  return(n)
}
lambda_samples |>
  mutate(across(c(2:101), prep_cv)) -> lambda_samples_cv

## variables 
variables <- colnames(data_all)[c(2:10)]

## run and save the baseline spatial model and select the variables for the final models
dir.create("temp")
dir.create("temp/data")
dir.create("temp/model")
for(i in 1:100){
  
  lambda_samples_cv |> dplyr::select(c(1,i+1)) |> `colnames<-` (c("EA.name","foi")) -> n_sample
  
  ## for each n remove 10% of EAs randomly
  foi_est <- n_sample[!is.na(n_sample$foi),]
  data_est <- data_all[data_all$EA.name %in% foi_est$EA.name, c(1:12)]
  data_est|>left_join(foi_est)->data_est
  data_pred <- data_all[!data_all$EA.name %in% foi_est$EA.name, c(1:12)]
  coords_est <- as.matrix(data_est[,c("lon","lat")])
  coords_pred <- as.matrix(data_all[,c("lon","lat")])
  
  # make projector matrix A
  A.est = inla.spde.make.A(mesh, loc=as.matrix(coords_est))
  A.pred = inla.spde.make.A(mesh, loc=as.matrix(coords_pred))
  
  # create the stack objects
  stack.pred = inla.stack(data=list(y=NA), 
                          A=list(A.pred,1),
                          effects=list(c(field.indices,
                                         list(Intercept=1)),
                                       list(pre_eemean = data_all$pre_eemean,
                                            ele_eemean = data_all$ele_eemean,
                                            tmax_eemean = data_all$tmax_eemean,
                                            tmin_eemean = data_all$tmin_eemean,
                                            rh_tmax_eemean = data_all$rh_tmax_eemean,
                                            rh_tmin_eemean = data_all$rh_tmin_eemean,
                                            mir_eemean = data_all$mir_eemean,
                                            evi_eemean = data_all$evi_eemean,
                                            pop_eemean = data_all$pop_eemean)),
                          tag="pred")
  
  stack.est<- inla.stack(data=list(y = log(data_est$foi)),
                         A=list(A.est,1),
                         effects=list(c(field.indices, 
                                        list(Intercept=1)),
                                      list(pre_eemean = data_est$pre_eemean,
                                           ele_eemean = data_est$ele_eemean,
                                           tmax_eemean = data_est$tmax_eemean,
                                           tmin_eemean = data_est$tmin_eemean,
                                           rh_tmax_eemean = data_est$rh_tmax_eemean,
                                           rh_tmin_eemean = data_est$rh_tmin_eemean,
                                           mir_eemean = data_est$mir_eemean,
                                           evi_eemean = data_est$evi_eemean,
                                           pop_eemean = data_est$pop_eemean)), 
                         tag="est",
                         remove.unused=TRUE,
                         compress=TRUE)
  
  stack = inla.stack(stack.est, stack.pred)
  temp_data <- inla.stack.data(stack)
  
  saveRDS(stack, paste0("temp/data/stack_", i, ".RDS"))
  saveRDS(n_sample, paste0("temp/data/nsample_", i, ".RDS"))
  saveRDS(data_est, paste0("temp/data/data_est_", i, ".RDS"))
  saveRDS(data_pred, paste0("temp/data/data_pred_", i, ".RDS"))
  
  ## baseline model
  formula_base <- y ~ -1 + Intercept + f(field, model=spde) 
  model_base <- inla(formula_base,
                     data=inla.stack.data(stack),
                     family="gaussian",
                     control.predictor = list(A=inla.stack.A(stack), compute=TRUE),
                     control.compute=list(waic =TRUE, cpo=TRUE, dic=TRUE),
                     keep=FALSE, verbose=FALSE)
  
  
  saveRDS(model_base, paste0("temp/model/baseline_model_", i, ".RDS"))
  
  # variable selection using forward selection process
  result <- INLAstep(fam1 = "gaussian", 
                     data_all,
                     in_stack = stack,
                     spde=spde,
                     invariant = "-1 + Intercept + f(field, model=spde)",
                     direction = 'forwards',
                     include = 2:10,
                     y = 'y',
                     powerl = 1,
                     inter = 1,
                     thresh = 4)
  
  saveRDS(result$best_model, paste0("temp/model/final_model_", i, ".RDS"))
  saveRDS(result$best_model$summary.fixed, paste0("temp/model/model_fixed_", i, ".RDS"))
  saveRDS(result$progress, paste0("temp/model/model_steps_", i, ".RDS"))
  
  print(i)
}


## run final models for each n sample with the chosen variables 
for(i in 1:100){ 
  var <- readRDS(paste0("temp/model/model_fixed_", i, ".RDS")) |> tibble::rownames_to_column(var = "RowNames")
  var <- var$RowNames[-1]
  stack <- readRDS(paste0("temp/data/stack_", i, ".RDS"))
  
  # formula of the final model with selected variables
  if(length(var)==1){
    formula <-y ~ -1 + Intercept + f(field, model=spde) + inla.stack.data(stack)[[`var`]]
  }
  if(length(var)==2){
    var1<-var[1]
    var2<-var[2]
    formula <-y ~ -1 + Intercept + f(field, model=spde) + inla.stack.data(stack)[[`var1`]]+ inla.stack.data(stack)[[`var2`]]
  }
  if(length(var)==3){
    var1<-var[1]
    var2<-var[2]
    var3<-var[3]
    formula <-y ~ -1 + Intercept + f(field, model=spde) + inla.stack.data(stack)[[`var1`]]+ inla.stack.data(stack)[[`var2`]]+ inla.stack.data(stack)[[`var3`]]
  }
  
  model_final <- inla(formula,
                      data=inla.stack.data(stack),
                      family="gaussian",
                      control.predictor = list(A=inla.stack.A(stack), compute=TRUE),
                      control.compute=list(waic =TRUE, cpo=T, dic=T, return.marginals.predictor=TRUE),
                      keep=FALSE, verbose=FALSE)
  
  saveRDS(model_final, paste0("temp/model/final_marg_model_", i, ".RDS"))
  print(i)
}

## FOI posteriors
model_list <- list()
files<-list.files("temp/model/", pattern="final_marg_model_")
for(i in 1:length(files)){
  model_list[[i]]<-readRDS(paste0("temp/model/",files[i]))
}

for(i in 1:100){ # for each sample
  model <- model_list[[i]]
  stack_data <- readRDS(paste0("temp/data/stack_",i, ".RDS"))
  pred_index <- inla.stack.index(stack_data,"pred")$data
  est_index <- inla.stack.index(stack_data,"est")$data
  ea_in_sample <- readRDS(paste0("temp/data/data_est_",i, ".RDS"))$EA.name
  ea_all <- readRDS(paste0("temp/data/data_pred_",i, ".RDS"))$EA.name
  
  list <- lapply(model$marginals.linear.predictor[pred_index], function(x){exp(unlist(x)[,1])})
  names(list) <- climate_data$EA.name 
  
  if(!exists("marg_foi")){
    marg_foi <- list(list)
  }else{
    marg_foi <- list(list, marg_foi)
  }
  keys <- unique(unlist(lapply(marg_foi, names)))
  marg_foi <- setNames(do.call(mapply, c(FUN=c, lapply(marg_foi, `[`, keys))), keys)
}

# calculate seroprevalence at age 9 from FOI
marg_sp <- lapply(marg_foi, function(x) 1-exp(-9*x) )

# extract quantiles of estimates per EA
quant_foi <- lapply(marg_foi, function(x)c((quantile(x,probs=c(0.5))), sd(x)))
quant_sp <- lapply(marg_sp, function(x)c((quantile(x,probs=c(0.5))), sd(x)))

pred_foi <- data.frame(EA=names(quant_foi),
                     pred=unlist(sapply(quant_foi, function(x)x[1])),
                     sd=unlist(sapply(quant_foi, function(x)x[2])),
                     pred_sp=unlist(sapply(quant_sp, function(x)x[1])),
                     sd_sp=unlist(sapply(quant_sp, function(x)x[2])))
pred_foi |> left_join(climate_data |> 
                        dplyr::select(lat,lon,EA.name)|> 
                        mutate(EA=as.character(EA.name))) -> pred_foi

# quantiles of the estimated total FOI and seroprevalence
quantile(pred_foi$pred, na.rm=T)
quantile(pred_foi$pred_sp, na.rm=T, probs=c(0.025,0.5,0.975))


########################################
## generating maps of FOI predictions ##
########################################

# read in and process shapefile for the whole of Timor Leste
codes <- read.csv("spatial_model_data/ea_codes.csv")
shapefile <- st_read("spatial_model_data/EAs_2015.kml")
shapefile$EA <- codes$ea_codes
shapefile$District <- codes$district_n                     

pred_foi|>left_join(codes |> mutate(EA=as.character(ea_codes))) -> pred_foi
pred_foi$district[pred_foi$lat>-8.4&pred_foi$lon<126]<-"Atauro"
saveRDS(pred_foi, "spatial_model_predicted_foi.RDS")

# join together the shapefile and the estimated FOI for maps
sf <- left_join(shapefile, pred_foi |> mutate(EA=as.integer(EA)))

# make maps of the central FOI estimates (continuous scale) and standard deviation
A <- ggplot(sf)+geom_sf(aes(fill=pred),colour=NA)+
  labs(fill="Predicted\nDENV\nFOI")+
  theme(rect = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())+
  scale_fill_viridis_c(option = "magma",begin = 0)

B <- ggplot(sf)+geom_sf(aes(fill=sd),colour=NA)+
  labs(fill="Standard\ndeviation")+
  theme(rect = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())+
  scale_fill_viridis_c(option = "magma",begin = 0)

cowplot::plot_grid(A,B, labels=c("A", "B"), ncol=1)

# make maps of the central FOI estimates (binned) and seroprevalence
pred_foi |> mutate(exceed = c(case_when(pred < 0.05 ~ "<0.05",
                                         pred >=0.05 & pred <0.1 ~ "0.05-0.1",
                                         pred >=0.1 & pred <0.15 ~ "0.1-0.15",
                                         pred >=0.15 & pred <0.2 ~ "0.15-0.2",
                                         pred >=0.2 & pred <0.25 ~ "0.2-0.25",
                                         pred>0.25 ~ ">0.25")))|>
  mutate(SP_bins = case_when(pred_sp<0.4 ~ "Under 40%",
                             pred_sp > 0.6 ~ "Above 60%",
                             pred_sp >=0.4 & pred_sp <= 0.6 ~ "40-60%")) -> pred_foi_exceed

sf <- left_join(shapefile, pred_foi_exceed |> mutate(EA=as.integer(EA)))

colours <-  scales::seq_gradient_pal("orange", "#440154")(seq(0,1,length.out=6))
colours_sp <- c("#FFEE99","#FFA700","#ED2939")

map <- ggplot(sf)+geom_sf(aes(fill=exceed),colour=NA)+
  labs(fill="Predicted\nDENV\nFOI")+
  theme(rect = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=colours, breaks = c("<0.05","0.05-0.1","0.1-0.15","0.15-0.2","0.2-0.25",">0.25"))

map_sp <- ggplot(sf)+geom_sf(aes(fill=SP_bins),colour=NA)+
  labs(fill="Predicted\nDENV\nseroprevalence\nat age 9")+
  theme(rect = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())+
  scale_fill_manual(values=colours_sp, breaks = c("Under 40%","40-60%","Above 60%"))

# aggregate the FOI predictions by municipality
median.quartile <- function(x){
  out <- quantile(x, probs = c(0.025,0.5,0.975))
  names(out) <- c("ymin","y","ymax")
  return(out) 
}
ggplot(sf)+
  stat_summary(aes(x=district,y=pred),fun.y="median",geom='point')+
  stat_summary(aes(x=district,y=pred),fun=median.quartile,geom='line')+
  theme_bw()+
  labs(x="Municipality",y="DENV FOI")+
  theme(text=element_text(size=14),
        axis.text.x = element_text(angle=45, hjust=1)) -> violin

cowplot::plot_grid(map,map_sp,violin,ncol=1,labels="AUTO")

# median FOI and sp per municipality
sf |> group_by(district) |> reframe(n=median(pred), ci_low=median.quartile(pred)[1], ci_upp=median.quartile(pred)[3])
1-exp(-9*0.05)
1-exp(-9*0.15)

## make probability of exceeding seroprevalence threshold maps 
quant_foi<-lapply(marg_foi, function(x) c(length(x[x>0.057])/length(x),
                                          length(x[x>0.102])/length(x)
)) 

pred_foi <- data.frame(EA=names(quant_foi),
                     T1=as.vector(unlist(sapply(quant_foi, function(x)x[1]))),
                     T2=as.vector(unlist(sapply(quant_foi, function(x)x[2]))))

pred_foi |> left_join(climate_data |> 
                        dplyr::select(lat,lon,EA.name)|> 
                        mutate(EA=as.character(EA.name))) -> pred_foi
pred_foi|>left_join(codes|>mutate(EA=as.character(ea_codes))) -> pred_foi
pred_foi$district[pred_foi$lat>-8.4 & pred_foi$lon<126]<-"Atauro"

sf <- left_join(shapefile, pred_foi |> mutate(EA=as.integer(EA)))

#P(serop9)>40% and P(serop9)>60%
colours <- scales::seq_gradient_pal("orange", "#440154")(seq(0,1,length.out=6))
AA <- ggplot(sf)+geom_sf(aes(fill=T1),colour=NA)+
  labs(fill="P(serop9)>40%")+
  theme(rect = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())+
  scale_fill_gradientn(colours=colours, limits=c(0,1))

A <- ggplot(sf)+geom_sf(aes(fill=T2),colour=NA)+
  labs(fill="P(serop9)>60%")+
  theme(rect = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())+
  scale_fill_gradientn(colours=colours, limits=c(0,1))

cowplot::plot_grid(A,AA, labels="AUTO", ncol=1)


######################
## model assessment ##
######################

## how many times was each variable chosen in the final model out of 100 iterations?
var_list<-list()
files<-list.files("temp/model/", pattern="model_fixed_")
for(i in 1:length(files)){
  var_list[[i]]<-readRDS(paste0("temp/model/",files[i]))
  var_list[[i]] <- var_list[[i]] |> tibble::rownames_to_column() |>
                    dplyr::mutate(order=1:n(), num=i) |> rename(var=rowname)
}

# filter to be the variables kept (order<3)
var_list |> bind_rows() |> filter(order==4) -> var4
var_list |> bind_rows() |> filter(order==3) -> var3
100 - nrow(var4) - nrow(var3)
var_list |> bind_rows() |> filter(!num %in% var4$num & !num %in% var3$num) |>
  group_by(var) |> summarise(n()) 
var_list |> bind_rows() |> filter(num %in% var3$num) |>
  group_by(var) |> summarise(n()) 
var_list |> bind_rows() |> 
  group_by(var) |> summarise(n())

## variable regression coefficients 
marg_rhmin<-marg_rhmax<-marg_evi<-list()
for(i in 1:100){
  model<- model_list[[i]]
  vars <- var_list[[i]]$var[-1]
  
  if(length(vars)==1){
    if("rh_tmin_eemean" == vars){
      marg_rhmin[[i]] <- model$marginals.fixed$`inla.stack.data(stack)[[var]]`[,1]
    }
    if("rh_tmax_eemean" == vars){
      marg_rhmax[[i]] <- model$marginals.fixed$`inla.stack.data(stack)[[var]]`[,1]
    }
    if("evi_eemean" == vars){
      marg_evi[[i]] <- model$marginals.fixed$`inla.stack.data(stack)[[var]]`[,1]
    }
  }
  
  if(length(vars)==2){
    if("rh_tmin_eemean" == vars[1]){
      marg_rhmin[[i]] <- model$marginals.fixed$`inla.stack.data(stack)[[var1]]`[,1]
    }
    if("rh_tmax_eemean" == vars[1]){
      marg_rhmax[[i]] <- model$marginals.fixed$`inla.stack.data(stack)[[var1]]`[,1]
    }
    if("evi_eemean" == vars[1]){
      marg_evi[[i]] <- model$marginals.fixed$`inla.stack.data(stack)[[var1]]`[,1]
    }
    if("rh_tmin_eemean" == vars[2]){
      marg_rhmin[[i]] <- model$marginals.fixed$`inla.stack.data(stack)[[var2]]`[,1]
    }
    if("rh_tmax_eemean" == vars[2]){
      marg_rhmax[[i]] <- model$marginals.fixed$`inla.stack.data(stack)[[var2]]`[,1]
    }
    if("evi_eemean" == vars[2]){
      marg_evi[[i]] <- model$marginals.fixed$`inla.stack.data(stack)[[var2]]`[,1]
    }
  }
  
}

marg_rhtmin_tot <- unlist(marg_rhmin)
marg_rhtmax_tot <- unlist(marg_rhmax)
marg_evi_tot <- unlist(marg_evi)

hist(marg_evi_tot)
hist(marg_rhtmin_tot)
hist(marg_rhtmax_tot)
quantile(marg_evi_tot, probs=c(0.5,0.025,0.975))
quantile(marg_rhtmin_tot, probs=c(0.5,0.025,0.975))
quantile(marg_rhtmax_tot, probs=c(0.5,0.025,0.975))


## Cross validation
posterior_medians <- readRDS("catalytic_model_posterior_medians.RDS")

marg_foi <- vector(mode='list', length=102)
names(marg_foi)<-colnames(lambda)

for(i in 1:100){ # for each sample
  model<- model_list[[i]]
  stack_data<-readRDS(paste0("temp/data/stack_",i, ".RDS"))
  pred_index<-inla.stack.index(stack_data,"pred")$data
  
  with_cv_in_sample <- readRDS(paste0("temp/data/data_est_",i, ".RDS"))
  posterior_medians |> filter(name %in% with_cv_in_sample$EA.name == F) -> cv_in_sample
  
  list <- lapply(model$marginals.linear.predictor[pred_index], function(x){exp(unlist(x)[,1])})
  names(list)<- climate_data$EA.name 
  list <- list[names(list) %in% cv_in_sample$name]
  
  marg_foi <- list(list, marg_foi)
  
  keys <- unique(unlist(lapply(marg_foi, names)))
  marg_foi <- setNames(do.call(mapply, c(FUN=c, lapply(marg_foi, `[`, keys))), keys)
}

quant_foi<-lapply(marg_foi, function(x)c((quantile(x[x<=1],probs=c(0.5, 0.025,0.975)))))
pred_foi<-data.frame(name=names(quant_foi),
                     pred=unlist(sapply(quant_foi, function(x)x[1])),
                     ciL_pred=unlist(sapply(quant_foi, function(x)x[2])),
                     ciU_pred=unlist(sapply(quant_foi, function(x)x[3])))

posterior_medians |> left_join(pred_foi) -> pred_cv

pred_cv |>
  ggplot()+
  geom_errorbar(aes(x=median, ymin= ciL_pred, ymax=ciU_pred), col="darkgrey")+
  geom_errorbar(aes(y=pred, xmin= ciL, xmax=ciU), col="darkgrey")+
  geom_point(aes(x=median, y=pred))+
  geom_abline(slope=1, col="navy", linetype="dashed", linewidth=1)+
  labs(y="Spatial model estimate",x="Sero-catalytic model estimate") + theme_bw() +
  coord_cartesian(c(0,1), c(0,1))+
  theme(text=element_text(size=15))
