########################################################
#
# SCRIPT TO RUN HIERARCHICAL BASED-MODEL IN TIMOR LESTE
#
########################################################

#Install and download the required libraries
library(rstan)
library(dplyr)
library(tidyverse)
library(truncnorm)
library(binom)
library(ggplot2)
library(purrr)

#Set working directory to you folder
setwd() #fill in
wd <- getwd()

# Load dataset
df_org <- read.csv(paste0(wd, "/data/Timor_dataset.csv"))
rur_urb <- read.csv(paste0(wd, "/data/Dengue_EA_urban_rural_assignment.csv"))

#Tidy information
colnames(rur_urb)[2] <- "EA"
colnames(rur_urb)[6] <- "class"
rur_urb$class <- ifelse(rur_urb$class == "Urban", 1, 0)

df_org <- left_join(df_org, rur_urb, by = "EA")

set.seed(1234)

EA <- unique(df_org$EA) #list of EA, 102
n_ea <- length(unique(df_org$EA)) #number of EA

df <- df_org 

N <- nrow(df) #number of individuals
df$household_id <- as.integer(as.factor(df$household)) #household IDs
H <- length(unique(df$household_id)) #number of households (HH)

df$village_id <- as.integer(as.factor(df$EA)) #village IDs
V <- length(unique(df$village_id)) #number of villages (EA)

dir.create(paste0(wd, "/optimised_Timor/")) #create output directory

source(paste0(wd,"/R/functions_sesp.R")) #R script with functions

run_foi_model(df=df, #dataframe
              wd=wd,  #input working directory
              n_ea = n_ea, #number of EA
              mean_sigma_n = 0.1, #mean distribution of sigma_n for prior
              mean_sigma_v = 0.1, #mean distribution of sigma_v for prior
              sd_sigma_n = 0.1, #standard deviation distribution of sigma_n for prior
              sd_sigma_v = 0.1, #standard deviation distribution of sigma_v for prior
              wd_out = paste0(wd, "/optimised_Timor/", "EA_102_sesp_log", "/"),   #output working directory
              n_it = 15000, #number of iterations
              model = "Model_TL.stan") #rstan model

source(paste0(wd,"/R/functions_rur.R"))

run_foi_model(df=df, 
              wd=wd,
              n_ea = n_ea,
              mean_sigma_n = 0.1,
              mean_sigma_v = 0.1,
              sd_sigma_n = 0.1,
              sd_sigma_v = 0.1,
              wd_out = paste0(wd, "/optimised_Timor/", "EA_102_sesp_rururb_log", "/"),
              n_it = 15000,
              model = "Model_TL_rururb.stan")

#Explore the outputs
df1 <- readRDS( paste0(wd, "/optimised_Timor/", "EA_102_sesp3_rururb_log", "/fit.rds"))
df2 <- readRDS( paste0(wd, "/optimised_Timor/", "EA_102_sesp_log", "/fit.rds"))

chain1 <- rstan::extract(df1)
chain2 <- rstan::extract(df2)

#associate village ID to village name
tag <- df %>%
  distinct(EA, village_id)


#Extract FOI
lambda1 <- lambda2  <- data.frame(matrix(NA, n_ea, 4))
for (l in 1: n_ea) lambda1[l,1:3] <- quantile(chain1$lambda_v[,l], c(0.5, 0.025, 0.975))
for (l in 1: n_ea) lambda2[l,1:3] <- quantile(chain2$lambda_v[,l], c(0.5, 0.025, 0.975))

lambda1[,4] <- 1:n_ea
colnames(lambda1) <- c("med", "ciL", "ciU", "village_id")
lambda1 <- left_join(lambda1, tag, by = "village_id")

lambda2[,4] <- 1:n_ea
colnames(lambda2) <- c("med", "ciL", "ciU", "village_id")
lambda2 <- left_join(lambda2, tag, by = "village_id")

lambda1$type <- "rural/urban"
lambda2$type <- "sesp"

lambda <- rbind(lambda1, lambda2)

#Plot FOI comparisons
ggplot(lambda, aes(x=as.factor(EA), y=med, ymin=ciL, ymax=ciU, col = type))+
  geom_point(position=position_dodge(width=0.5))+
  geom_errorbar(position=position_dodge(width=0.5), width =0.2)+
  theme_bw() + 
  xlab("EA") +ylab("FOI (95% CrI)")+ labs(col = "Model")+# White background theme
  theme(
    plot.title = element_text(hjust = 0.5),      # Center title
    axis.text.x = element_text(angle = 45,       # Rotate x-axis labels
                               hjust = 1)
  )

#extract test sensitivity and specificity
sesp1 <- data.frame(matrix(NA, 2, 4))
sesp1[1,1:3] <- quantile(chain1$se, c(0.5, 0.025, 0.975))
sesp1[1,1:3] <- quantile(chain1$sp, c(0.5, 0.025, 0.975))
sesp1[,4] <- c("se", "sp")
sesp1$tupe <- "rural/urban"

sesp2 <- data.frame(matrix(NA, 2, 4))
sesp2[1,1:3] <- quantile(chain2$se, c(0.5, 0.025, 0.975))
sesp2[1,1:3] <- quantile(chain2$sp, c(0.5, 0.025, 0.975))
sesp2[,4] <- c("se", "sp")
sesp2$tupe <- "sesp"

#extract urban effect from Model rural urban
beta <- data.frame(matrix(NA, 1, 4))
beta[1,1:3] <- quantile(chain1$beta, c(0.5, 0.025, 0.975))
beta



