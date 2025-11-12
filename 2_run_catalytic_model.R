###########################################################
#
#      Code to estimate the FOI in Timor-Leste 
#
###########################################################

#Load libraries
library(dplyr)
library(tidyverse)
library(rstan)
library(Hmisc)
library(gridExtra)
library("ggpubr")

#Load data
wd1 <- getwd()
wd2 <- paste0(wd1, "/data/")
df_multiple_tested_sing_pos <- read.csv(paste0(wd2, "Timor_single_pos.csv")) 
df_multiple_tested_multi_pos  <- read.csv(paste0(wd2, "Timor_multi_pos.csv")) 
df_multiple_tested_neg  <- read.csv(paste0(wd2, "Timor_multipleHHneg_singleHH.csv"))

#source functions to plot
source(paste0(wd1,"/R/functions.R"))

#create output directory
dir.create(file.path(wd1, "result"))
wd3 <- paste0(wd1, "/result/")

#Create a single dataframe with the entire dataset
complete_data <- rbind(df_multiple_tested_sing_pos, df_multiple_tested_multi_pos, df_multiple_tested_neg)

# create a list of EA 
EA_list <- unique(complete_data$EA)
n_EA <- length(EA_list)
n_HH <- length(unique(complete_data$household))

# Group by household and age, sum individuals
aggregated_complete_data <- complete_data %>%
  group_by(EA, age, tag) %>%
  summarise(total_positive = sum(positive_count, na.rm = TRUE),
            total_tested = sum(total_count, na.rm = TRUE)) %>%
  ungroup() %>%
  complete(EA, age, tag, fill = list(total_positive = 0, total_tested = 0))

aggregated_complete_data.tag1 <- filter(aggregated_complete_data, tag == 1) #assigned to Eq2 
aggregated_complete_data.tag0 <- filter(aggregated_complete_data, tag == 0) #assigned to Eq1

age <- unique(aggregated_complete_data$age)
AgeG <- length(age)

#Work with tag 1 --> #assigned to Eq2 (with bH)
Npos.tag1 <- aggregated_complete_data.tag1  %>%
select(EA, age, total_positive) %>%
  pivot_wider(
    names_from = age,          #age cols
    values_from = total_positive,  #tot pos in cells
    values_fill = list(total_positive = 0)   # fill with 0 if no value
  )

#Work with tested among the tag 1
N.tag1 <- aggregated_complete_data.tag1  %>%
  select(EA, age, total_tested) %>%
  pivot_wider(
    names_from = age,           #age cols
    values_from = total_tested,  #tot pos in cells
    values_fill = list(total_tested = 0)  # fill with 0 if no value
  )

#Format better
NPos_list1 <- Npos.tag1[,-1]
N_list1 <- N.tag1[,-1]

#Work with tag 2 --> #assigned to Eq1 (NO bH) 
Npos.tag0 <- aggregated_complete_data.tag0  %>%
  select(EA, age, total_positive) %>%
  pivot_wider(
    names_from = age,           
    values_from = total_positive,  
    values_fill = list(total_positive = 0)  
  )

#Work with tested among tag 2
N.tag0 <- aggregated_complete_data.tag0  %>%
  select(EA, age, total_tested) %>%
  pivot_wider(
    names_from = age,           
    values_from = total_tested, 
    values_fill = list(total_tested = 0) 
  )

#Format better
NPos_list0 <- Npos.tag0[,-1]
N_list0 <- N.tag0[,-1]


#Save all information into a list
data <- list(NPos1 = NPos_list1, #assigned to eq 2, positive
             N1 = N_list1, #assigned to eq 2, tested
             NPos0 = NPos_list0, #assigned to eq 1, positive
             N0 = N_list0, #assigned to eq 1, tested
             age = age, #age
             AgeG = AgeG, #number of age group
             EA = n_EA) #number of EA

#Load rstan model 
mod <- stan_model(paste0(wd1,"/R/model_TL.stan"))

#set initial conditions if necessary
# init_pars <- list()
num_lambdas <- 1
n_chains<-3

# init_pars <- list(
#   list(lambda = rep(0.1, n_EA)),  # Chain 1
#   list(lambda = rep(0.2, n_EA)),  # Chain 2
#   list(lambda = rep(0.3, n_EA))   # Chain 3
# )

# Fit the model using rstan::sampling()
stanfit1 <- rstan::sampling(
  object = mod,
  data = data,               # Data list
  chains = 3,                # Number of chains
  iter = 15000,              # Total iterations 
  warmup = 5000,             # Number of warmup iterations
  # init = init_pars,          # Initial values
  cores = 3,                 # Number of cores for parallel chains
  refresh = 50               # Frequency of progress reporting
)

saveRDS(stanfit1, file = paste0(wd3, "output_TL.rds"))
# stanfit1 <- readRDS(paste0(wd3, "output_TL.rds"))

##############################################################
#
#                     POSTERIOR CHECKS
#
##############################################################

# Extract chains (equivalent to cmdstanr::fit$draws())
chains <- rstan::extract(stanfit1)

#get EA list
EA_list <- unique(aggregated_complete_data$EA)

#Extract the clustering effect (b_H) and its stat
val_b_H <- c(chains[["b_H"]])
b_H <- as.data.frame(t(quantile(val_b_H , c(0.5, 0.025, 0.975))))

#Extract EA_specific FOI and its stat
val_lambda <- as.data.frame(matrix(0, n_EA, 30000))
lambda <- as.data.frame(matrix(0, n_EA, 3))
for (i in 1:n_EA) {
  lambda[i,] <- quantile(chains$lambda[,i], c(0.5, 0.025, 0.975))
}

lambda$EA <- unique(EA_list)
val_lambda <- as.data.frame(val_lambda)
for (i in 1:n_EA) { val_lambda[i,] <- chains$lambda[,i]}
val_lambda$EA <- EA_list

#Save outputs
write.csv(b_H, paste0(wd3, "/bH.csv"))

write.csv(lambda, paste0(wd3, "/lambda.csv"))

write.csv(val_lambda, paste0(wd3, "lambda_posterior.csv"))

write.csv(val_b_H, paste0(wd3, "bH_posterior.csv"))

#Plot traceplots
ggsave(
  paste0("traceplot_bh.png"),
  plot = traceplot(stanfit1, "b_H"),
  device = "png",
  path = wd3,
  scale = 1,
  width = 100,
  units = "mm",
  dpi = 300
)

ggsave(
  paste0("traceplot_lambda.png"),
  plot = traceplot(stanfit1, "lambda"),
  device = "png",
  path = wd3,
  scale = 5,
  height = 100,
  width = 100,
  units = "mm",
  dpi = 300
)


# PRIOR: truncated Normal(0, 1) with lower bound 0
set.seed(123)
n_samples <- 30000
prior_samples <- rnorm(n_samples, mean = 0.5, sd=0.5)
prior_samples <- prior_samples[prior_samples > 0]  # Troncamento a 0
prior_long <- data.frame(
  sample = 1:length(prior_samples),
  value = prior_samples,
  type = "prior"
)

# POSTERIOR
posterior_long <-data.frame(
  sample = seq(1, length(val_b_H)),
  value = val_b_H,
  type = "posterior"
)

combined_df <- bind_rows(prior_long, posterior_long)

#check prior posterior bH
plot <- ggplot(combined_df, aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("green", "red")) +
  labs(title = NULL,
       x = "Value",
       y = "Density",
       fill = "Distribution") +
  theme_minimal()

ggsave(
  paste0("multiplot_bH_posterior.png"),
  plot = plot,
  device = "png",
  path = wd3,
  scale = 1.5,
  height = 100,
  width = 150,
  units = "mm",
  dpi = 300
)

##############################################################
#
#                             PLOTS
#
##############################################################

### Plot goodness of fit
age_labels <- c("(0,10]", "(10,20]", "(20,30]", "(30,40]", "(40,50]", "(50,60]", "(60+)")
age_breaks <- c(0,9,19,29,39,49,59)

#organise data
complete_data_plot <- complete_data
complete_data_plot$age_group <- cut(complete_data_plot$age, 
                                    breaks = c(age_breaks, Inf),
                                    include.lowest = TRUE, 
                                    labels = age_labels)

#rename age groups
complete_data_plot <- complete_data_plot %>%
  group_by(age_group, EA) %>%
  dplyr::summarize(
    total_tested = sum(total_count),
    total_positive = sum(positive_count),
    age_min = case_when(age_group == "(0,10]" ~ 0,
                        
                        age_group == "(10,20]" ~ 10,
                        
                        age_group == "(20,30]" ~ 20,
                        
                        age_group == "(30,40]" ~ 30,
                        
                        age_group == "(40,50]" ~ 40,
                        
                        age_group == "(50,60]" ~ 50,
                        
                        age_group == "(60+)" ~ 60),
    
    age_max = case_when(age_group == "(0,10]" ~ 9,
                        
                        age_group == "(10,20]" ~ 19,
                        
                        age_group == "(20,30]" ~ 29,
                        
                        age_group == "(30,40]" ~ 39,
                        
                        age_group == "(40,50]" ~ 49,
                        
                        age_group == "(50,60]" ~ 59,
                        
                        age_group == "(60+)" ~ 102)) %>%
  distinct()

#calcualte mid age point
complete_data_plot$age <- (complete_data_plot$age_min + complete_data_plot$age_max) / 2

#number of household EA
nn_EA <- dim(val_lambda)[1]
#number og age groups
N_AgeG <- length(complete_data_plot$age)

true_IC_serop<- list()

#Calculate OBSERVED seroprevalence
for (k in 1: n_EA) {
  print(k)
  data <- filter(complete_data_plot, EA == EA_list[k])
  true_serop <- data.frame(matrix(NA, length(unique(data$age)), 5))
  true_serop[,1:3] <- binconf(x=data$total_positive, data$total_tested, alpha=.05)
  true_serop[,4] <- EA_list[k]
  true_serop[,5] <- unique(data$age)
  true_IC_serop[[k]] <- true_serop
}

#extract 1000 random samples from lambda/FOI and bH posterior and calculate ESTIMATED seroprevalence
nsamples <- 1000
val_lambda <- val_lambda[,-30001]

lambda <- as.data.frame(matrix(NA, nn_EA, nsamples))
nsamp <- sample(seq(1,dim(val_lambda)[2],1), nsamples)

b_h <- val_b_H[nsamp]
for (i in 1:nn_EA) lambda[i,] <- val_lambda[i,c(nsamp)]

age <- complete_data_plot$age
max_age <- 102
simulated_seroprevalence <- array(NA,dim = c(max_age , nn_EA, nsamples ))

for (k in 1:nn_EA){
  for (i in 1:nsamples){
    for (a in 1:max_age){ #each age group
      
      simulated_seroprevalence[a,k,i] <- 1-exp(-lambda[k,i] * (a + b_h[i])) # calculate seroprevalence with estimated lambda
      
    }}}


serop_list <- list()

for (k in 1:nn_EA){
  simulated_seroprevalence_stat <- matrix(NA,max_age , 5 )
  for(a in 1: max_age){
    simulated_seroprevalence_stat[a,1:3] <- quantile(simulated_seroprevalence[a,k,] , c(0.5, 0.025, 0.975))
    simulated_seroprevalence_stat[,4] <- EA_list[k]
    simulated_seroprevalence_stat[,5] <- seq(1,max_age ,1)
    serop_list[[k]] <- simulated_seroprevalence_stat
  }}

sim_serop <- as.data.frame(do.call(rbind,serop_list))
true_serop <- as.data.frame(do.call(rbind,true_IC_serop))
colnames(sim_serop) <- c("med", "ciL", "ciU", "EA", "age")
colnames(true_serop) <- c("med", "ciL", "ciU", "EA", "age")


plot_serop <- plot_f(sim_serop,true_serop, NULL) 
plot_serop

ggsave(
  paste0("multiplot_EA_serop.png"),
  plot = plot_serop,
  device = "png",
  path = wd3,
  scale = 3.7,
  height = 100,
  width = 100,
  units = "mm",
  dpi = 300
)

#Plot FOI by EA
lambda <- read.csv(paste0(wd3, "lambda.csv"))
colnames(lambda) <- c("X", "med", "ciL", "ciU", "EA")
complete_data_EA <- complete_data[, c("district", "EA")]
lambda <- left_join(lambda, complete_data_EA, by ="EA")
lambda$EA <- gsub("^", "EA_", lambda$EA)

plot_FOI <- ggplot(data = lambda, aes(x=EA, y= med, ymin=ciL, ymax=ciU, col = district)) +
  geom_point(
    size = 2)+
  geom_errorbar(width= 0.5)+
  ylim(0,1) +
  labs(x = element_blank(), y = element_blank())+
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab(NULL) + ylab("FOI")

ggsave(
  paste0("multiplot_EA_FOI.png"),
  plot = plot_FOI,
  device = "png",
  path = wd3,
  scale = 2,
  height = 80,
  width = 200,
  units = "mm",
  dpi = 300
)
