#########################################
## baseline spatial model preparation ###
#########################################

## packages
library(INLA)
library(ggplot2)
library(dplyr)
library(sf)
library(RColorBrewer)
library(lattice)
library(inlabru)
library(raster)
library(cowplot)

#################
### DATA PREP ###
#################

## scaling
standard <- function(x){ (x-mean(x,na.rm=T))/ sd(x, na.rm=T)}

## READ IN DATA
data <- read.csv("spatial_model_data/climate_FOI.csv")

## READ IN COORDS
coords <- read.csv("spatial_model_data/EA_centroids.csv")
coords$ea_codes <- data$ea_codes

## MATCH COORDS AND EAS
data_all <- data |> left_join(coords)
data_est <- data_all[!is.na(data_all$sampled),]
data_pred <- data_all[is.na(data_all$sampled),]


#####################
## CREATE THE MESH ##
#####################

coords_all <- as.matrix(data_all[,19:20])
coords <- as.matrix(data_est[,19:20])

## boundary of shapefile
shapefile <- st_read("spatial_model_data/EAs_2015.kml")
shapefile_outline <- readRDS("spatial_model_data/TimorLeste_2015_outline.RDS")
plot(shapefile_outline)

## mesh
mesh_n1 = inla.mesh.2d(loc=coords,
                       cutoff=0.03, offset=c(0.2,0.3),max.edge=c(0.2,0.4),min.angle=c(25,25),
                       boundary = NULL
)
mesh_n2 = inla.mesh.2d(loc=coords,
                       cutoff=0.03, offset=c(0.2,0.3),max.edge=c(0.2,0.4),min.angle=c(35,35),
                       boundary = NULL
)
mesh_n3 = inla.mesh.2d(loc=coords,
                       cutoff=0.03, offset=c(0.2,0.3),max.edge=c(0.2,0.4),min.angle=c(15,15),
                       boundary = NULL
)
mesh_n4 = inla.mesh.2d(loc=coords,
                       cutoff=0.03, offset=c(0.2,0.3),max.edge=c(0.5,0.5),min.angle=c(25,25),
                       boundary = NULL
)
mesh_n5 = inla.mesh.2d(loc=coords,
                       cutoff=0.03, offset=c(0.2,0.6),max.edge=c(0.2,0.4),min.angle=c(25,25),
                       boundary = NULL
)
mesh_n6 = inla.mesh.2d(loc=coords,
                       cutoff=0.03, offset=c(0.6,0.6),max.edge=c(0.2,0.4),min.angle=c(25,25),
                       boundary = NULL
)


mydata_sf<-st_as_sf(data_all, coords = c("lon", "lat"), crs = 2100)
dist_df<-as.data.frame(st_distance(mydata_sf))
range(dist_df[lower.tri(dist_df)])
max.edge = diff(range(dist_df))/(3*5)
bound.outer = diff(range(dist_df))/5

country.bdry <-  shapefile_outline %>% inla.sp2segment()

mesh_n7 <- inla.mesh.2d(boundary = country.bdry,
                        min.angle = c(25,25),
                        max.edge = c(0.2,0.4),
                        offset=c(0.2,0.6),
                        cutoff = 0.03) 
mesh_n8 <- inla.mesh.2d(boundary = country.bdry,
                        min.angle = c(25,25),
                        max.edge = c(0.2,0.4),
                        offset=c(0.2,0.6),
                        cutoff = 0.02) 

A<-ggplot() +
  geom_sf(data=shapefile,fill='lightblue',col='transparent', alpha=0.5)+  
  gg(mesh_n1) +
  geom_point(aes(x=coords[,1],y=coords[,2]),col='black',size=1.7,alpha=0.5) + theme_bw()
B<-ggplot() +
  geom_sf(data=shapefile,fill='lightblue',col='transparent', alpha=0.5)+  
  gg(mesh_n2) +
  geom_point(aes(x=coords[,1],y=coords[,2]),col='black',size=1.7,alpha=0.5) + theme_bw()
C<-ggplot() +
  geom_sf(data=shapefile,fill='lightblue',col='transparent', alpha=0.5)+  
  gg(mesh_n3) +
  geom_point(aes(x=coords[,1],y=coords[,2]),col='black',size=1.7,alpha=0.5) + theme_bw()
D<-ggplot() +
  geom_sf(data=shapefile,fill='lightblue',col='transparent', alpha=0.5)+  
  gg(mesh_n4) +
  geom_point(aes(x=coords[,1],y=coords[,2]),col='black',size=1.7,alpha=0.5) + theme_bw()
E<-ggplot() +
  geom_sf(data=shapefile,fill='lightblue',col='transparent', alpha=0.5)+  
  gg(mesh_n5) +
  geom_point(aes(x=coords[,1],y=coords[,2]),col='black',size=1.7,alpha=0.5) + theme_bw()
FF<-ggplot() +
  geom_sf(data=shapefile,fill='lightblue',col='transparent', alpha=0.5)+  
  gg(mesh_n6) +
  geom_point(aes(x=coords[,1],y=coords[,2]),col='black',size=1.7,alpha=0.5) + theme_bw()
G<-ggplot() +
  geom_sf(data=shapefile,fill='lightblue',col='transparent', alpha=0.5)+  
  gg(mesh_n7) +
  geom_point(aes(x=coords[,1],y=coords[,2]),col='black',size=1.7,alpha=0.5) + theme_bw()
H<-ggplot() +
  geom_sf(data=shapefile,fill='lightblue',col='transparent', alpha=0.5)+  
  gg(mesh_n8) +
  geom_point(aes(x=coords[,1],y=coords[,2]),col='black',size=1.7,alpha=0.5) + theme_bw()


plot_grid(A,B,C,D,E,FF,G,H, ncol=3, nrow=3, labels=c("A","B","C","D","E","F","G","H"))

mesh_list<-list(mesh_n1, mesh_n2,mesh_n3,mesh_n4,mesh_n5,mesh_n6,mesh_n7,mesh_n8)

##################
## set priors ####
##################

mesh <- mesh_list[[8]]
spde <- inla.spde2.pcmatern(mesh,
                          alpha=2,# smoothness parameter
                          prior.range = c(0.001,0.1), # probability of <10% that range less than minimum distance between EA centroids 
                          prior.sigma = c(1,0.1)) # probability of <5% that sd of range is greater than 1 

########################
## save mesh & spde ####
########################

saveRDS(spde, "spatial_model_data/saved_spde.RDS")
saveRDS(mesh_n9, "spatial_model_data/saved_mesh.RDS")


##################################################
## test baseline models with different meshes ####
##################################################

## read in posteriors from catalytic model
lambda <- as.data.frame(readRDS("spatial_model_data/FOI_catalytic_model_posteriors.RDS"))
tot_posterior <- dim(lambda)[1]

## sample n=100 from the catalytic model posteriors
int <- sample.int(size = 100, n = tot_posterior)
lambda |> dplyr::slice(int) -> lambda_samples

lambda_samples |> 
  pivot_longer(1:102) |>
  group_by(name)|>reframe(median=quantile(value, probs=c(0.5)),
                          ciL=quantile(value, probs=c(0.025)),
                          ciU=quantile(value, probs=c(0.975))) -> sample_medians

sample_medians |> mutate(name=as.integer(name)) |>
  rename(ea_codes ="name") |> left_join(data_est) -> data_est

store_waic <- list()
for(i in 1:8){
mesh<- mesh_list[[i]]
spde<-inla.spde2.pcmatern(mesh,
                          alpha=2,
                          prior.range = c(0.001,0.1), 
                          prior.sigma = c(1,0.1)) 
A.est = inla.spde.make.A(mesh, loc=as.matrix(coords))
A.pred = inla.spde.make.A(mesh, loc=as.matrix(coords_all))

# spatial field indices
field.indices =inla.spde.make.index("field", n.spde=mesh$n)

# create the stack object
stack.pred = inla.stack(data=list(y=NA), 
                        A=list(A.pred,1),
                        effects=
                          list(c(field.indices,
                                 list(Intercept=1)),
                               list(pre_eemean = standard(data_all$pre_eemean),
                                    ele_eemean = standard(data_all$ele_eemean),
                                    tmax_eemean = standard(data_all$tmax_eemean),
                                    tmin_eemean = standard(data_all$tmin_eemean),
                                    rh_tmax_eemean = standard(data_all$rh_tmax_eemean),
                                    rh_tmin_eemean = standard(data_all$rh_tmin_eemean),
                                    mir_eemean = standard(data_all$mir_eemean),
                                    evi_eemean = standard(data_all$evi_eemean),
                                    pop_eemean = standard(log(data_all$pop_eemean)))),
                        tag="pred")

stack.est<- inla.stack(data=list(y = log(data_est$median)),
                       A=list(A.est,1),
                       effects=list(c(field.indices, 
                                      list(Intercept=1)),
                                    list(pre_eemean = standard(data_est$pre_eemean),
                                         ele_eemean = standard(data_est$ele_eemean),
                                         tmax_eemean = standard(data_est$tmax_eemean),
                                         tmin_eemean = standard(data_est$tmin_eemean),
                                         rh_tmax_eemean = standard(data_est$rh_tmax_eemean),
                                         rh_tmin_eemean = standard(data_est$rh_tmin_eemean),
                                         mir_eemean = standard(data_est$mir_eemean),
                                         evi_eemean = standard(data_est$evi_eemean),
                                         pop_eemean=standard(log(data_est$pop_eemean)))), 
                       tag="est",
                       remove.unused=TRUE,
                       compress=TRUE)

stack = inla.stack(stack.est, stack.pred)

formula_base <- y ~ -1 + Intercept + f(field, model=spde) 
model_base <- inla(formula_base,
                   data=inla.stack.data(stack),
                   family="gaussian",
                   control.predictor = list(A=inla.stack.A(stack), compute=TRUE),
                   control.compute=list(waic =TRUE, cpo=T, dic=T),
                   keep=FALSE, verbose=F)
store_waic[[i]] <- model_base$waic$waic
}

store_waic

