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

###############
### DATA PREP #
###############

# scaling
standard <- function(x){ (x-mean(x,na.rm=T))/ sd(x, na.rm=T)}

## READ IN DATA
data <- read.csv("spatial_model_data/climate_FOI.csv")

## READ IN COORDS
coords<- read.csv("spatial_model_data/centroids.csv")
coords$ea_codes <- data$ea_codes

## MATCH COORDS AND EAS
data_all <- data |> left_join(coords)
data_est <- data_all[!is.na(data_all$median),]
data_pred <- data_all[is.na(data_all$median),]


#####################
## CREATE THE MESH ##
#####################

coords_all <- as.matrix(data_all[,23:24])
coords <- as.matrix(data_est[,23:24])

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
bnd <- inla.nonconvex.hull(coords, resolution=c(100), convex = 0.25)
mesh_n7 = inla.mesh.2d(
  cutoff=0.03, offset=c(0.2,0.6),max.edge=c(0.2,0.4),min.angle=c(25,25),
  boundary = bnd
)

mydata_sf<-st_as_sf(data_all, coords = c("lon", "lat"), crs = 2100)
dist_df<-as.data.frame(st_distance(mydata_sf))
range(dist_df[lower.tri(dist_df)])
max.edge = diff(range(dist_df))/(3*5)
bound.outer = diff(range(dist_df))/5

country.bdry <-  shapefile_outline %>% inla.sp2segment()
mesh_n7 <- inla.mesh.2d(boundary = country.bdry,
                         min.angle = c(25,25),
                      max.edge = c(1,2)*max.edge,
                      offset=c(max.edge,bound.outer),
                      cutoff = max.edge/8) # reduce cutoff so smaller triangles at boundary edge
mesh_n8 <- inla.mesh.2d(boundary = country.bdry,
                        min.angle = c(25,25),
                        max.edge = c(0.2,0.4),
                        offset=c(0.2,0.6),
                        cutoff = 0.03) 
mesh_n9 <- inla.mesh.2d(boundary = country.bdry,
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
I<-ggplot() +
  geom_sf(data=shapefile,fill='lightblue',col='transparent', alpha=0.5)+  
  gg(mesh_n9) +
  geom_point(aes(x=coords[,1],y=coords[,2]),col='black',size=1.7,alpha=0.5) + theme_bw()


plot_grid(A,B,C,D,E,FF,G,H, I,ncol=3, nrow=3, labels=c("A","B","C","D","E","F","G","H","J"))

mesh_list<-list(mesh_n1, mesh_n2,mesh_n3,mesh_n4,mesh_n5,mesh_n6,mesh_n7,mesh_n8,mesh_n9)

##################
## set priors ####
##################

mesh <- mesh_list[[9]]
spde <- inla.spde2.pcmatern(mesh,
                          alpha=2,# smoothness parameter
                          prior.range = c(0.001,0.1), # probability of <10% that range less than minimum distance between EA centroids 
                          prior.sigma = c(1,0.1)) # probability of <5% that sd of range is greater than 1 

########################
## save mesh & spde ####
########################

saveRDS(spde, "spatial_model_data/saved_spde.RDS")
saveRDS(mesh_n9, "spatial_model_data/saved_mesh.RDS")
