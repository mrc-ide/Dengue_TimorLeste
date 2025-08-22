
#LOAD LIBRARIES

library(dplyr)
library(tidyverse)
library(rstan)
library(ggplot2)
library("ggpubr")
library(binom)
library(Hmisc)
library(gridExtra)


#Set working directory
wd1 <- getwd()

dir.create(file.path(wd1, "data"))

#Load Original dataset
data_hous <- read.csv("Dataset_TL.csv") 

# Obtain the unique EA values for each district
unique_ea_list <- lapply(split(data_hous$EA, data_hous$district), unique)

# Create a dataset to store useful info
df <- data.frame(matrix(NA,  length(data_hous$EA),2))
df$EA <- data_hous$EA
df$age <- data_hous$age
df$result <- data_hous$result
df$household <- data_hous$household
df$district <- data_hous$district

#Grouping into age groups  of 10yrs
df$age_group <- cut(df$age, breaks = c(0, 10, 20, 30, 40, 50, 60, max(df$age)),
                    include.highest = FALSE)

#Select only columns you are interested in
df <- df[,c(3,4,5,6,7)]

#remove eventual NA (test results not given)
na_rows <- which(!complete.cases(df))
df <- df[-na_rows, ]

# Group by EA, age group, and result, and summarize the counts (n people positive =1 and n people negative (0))
summary_df <- df %>%
  group_by(EA, household, district, result, age)  %>%
  dplyr::summarize(
    count = n()
  )

# substitute 1 as Positive and 0 as Negative
for (i in 1:length(summary_df$result)) {
  if(summary_df$result[i] == "1") {
    summary_df$result[i] <- "Positive"
  } else if(summary_df$result[i] == "0") {
    summary_df$result[i] <- "Negative"
  } else {print(paste0("error_line",i))}
}

# Create a new column 'Household_Status' with initial value 'HH_Neg'
summary_df$Household_Status <- "HH_Neg"

# # Group the data by household
grouped_data <- summary_df %>%
  group_by(household)

# Check if any individual in the household (HH) has a positive result and if yes substitute with HH_Pos (also for single-individual household)
household_status <- grouped_data %>%
  mutate(Household_Status = case_when(
    sum(result == "Positive") > 1 ~ "HH_multi_Pos",   ##more than one pos individuals in a HH
    sum(result == "Positive") == 1 ~ "HH_single_Pos",  # only one pos individual in a HH
    TRUE ~ "HH_Neg"                                   # all negative
  ))


# Ungroup the data
ungrouped_data <- ungroup(household_status)

# Pivot the summary table to wide format
summary_wide <- pivot_wider(ungrouped_data, names_from = result, values_from = count, values_fill = 0)

# Calculate the total number of people tested and the number of positive results for each EA and age group
df <- summary_wide %>%
  group_by(EA, district, household, Household_Status, age) %>%
  dplyr::summarize(
    positive_count = sum(`Positive`),
    negative_count = sum(`Negative`),
    total_count = sum(`Positive`) + sum(`Negative`)
  )


#check tot number tests per EA 
df_summary <- df %>% 
  group_by(EA) %>% 
  dplyr::summarize(total_tested = sum(total_count))
# df_summary$EA <- paste("EA", df_summary$EA, sep = "_")
# df_summary <- filter(df_summary, total_tested > 30) # remove those with lower than 30 samples tested --> NOT done in the final analysis (commented)
# EA_accepted <- unique(df_summary$EA)
# df<- df[df$EA %in% EA_accepted,] 

# Dataframe with households where only a single individual was tested
household_count_list <-  group_by(df, household) %>% #individuate the number of individuals tested per house
  dplyr::summarize(tot_memb_tested = sum(`total_count`))

EA_single_list <- filter(household_count_list, tot_memb_tested <= 1) #407 household single tested

df_single_tested <- filter(df, household %in% EA_single_list$household) #keep only the household with single tested


# Dataframe with households where multiple individuals were tested and all resulted negative
df_multiple_tested <- df[!(df$household %in% df_single_tested$household), ] #remove single tested

df_multiple_tested_neg <- filter(df_multiple_tested, Household_Status == "HH_Neg")

df_multiple_tested_neg <- rbind(df_multiple_tested_neg, df_single_tested) #it contains also single tested

df_multiple_tested_sing_pos <- filter(df_multiple_tested, Household_Status == "HH_single_Pos")
df_multiple_tested_multi_pos <- filter(df_multiple_tested, Household_Status == "HH_multi_Pos")

#Add tag for pos individual in multi tested household with single positive
df_multiple_tested_sing_pos <- df_multiple_tested_sing_pos %>%
  mutate(tag = ifelse(positive_count > 0, 1, 0)) 
df_multiple_tested_neg$tag <- 0 #add tag 0 to neg individuals in neg HH or single tested HH
df_multiple_tested_multi_pos$tag <- 1 #add tag 1 to HH with multiple positive

#save outputs
wd2 <- paste0(wd1,"/data/")
write.csv(df_multiple_tested_sing_pos, paste0(wd2, "Timor_single_pos.csv"))
write.csv(df_multiple_tested_multi_pos, paste0(wd2, "Timor_multi_pos.csv"))
write.csv(df_multiple_tested_neg, paste0(wd2, "Timor_multipleHHneg_singleHH.csv"))

df_multiple_tested_pos <- rbind(df_multiple_tested_sing_pos, df_multiple_tested_multi_pos)
#total dataframe
tot.df <- rbind(df_multiple_tested_pos, df_multiple_tested_neg)

#household with multiple tested and at least one pos
a <-  group_by(df_multiple_tested_pos, household) %>% #individuate the number of individuals tested per house
  dplyr::summarize(tot_memb_pos = sum(`positive_count`))

#household with multiple tested and multiple positive
a.pos <- filter(a, tot_memb_pos > 1)
a.neg <- filter(a, tot_memb_pos == 1)

HH.pos <- unique(a.pos$household)
HH.neg <- c(unique(a.neg$household), unique(df_multiple_tested_neg$household))

tot.df.pos <- filter(tot.df, household %in% HH.pos)
tot.df.neg <- filter(tot.df, household %in% HH.neg)

tot.df.pos <- tot.df.pos %>%
  group_by(age) %>%
  summarise(
    total_positive = sum(positive_count, na.rm = TRUE),
    total_negative = sum(negative_count, na.rm = TRUE)
  ) %>%
  rowwise() %>%
  mutate(
    total = total_positive + total_negative,
    seroprevalence = total_positive / total * 100,
    ci_lower = binom.confint(total_positive, total, methods = "exact")$lower * 100,
    ci_upper = binom.confint(total_positive, total, methods = "exact")$upper * 100
  ) %>%
  ungroup()

tot.df.pos$tag <- "HH with other DENV IgG+ individuals"

tot.df.neg <- tot.df.neg %>%
  group_by(age) %>%
  summarise(
    total_positive = sum(positive_count, na.rm = TRUE),
    total_negative = sum(negative_count, na.rm = TRUE)
  ) %>%
  rowwise() %>%
  mutate(
    total = total_positive + total_negative,
    seroprevalence = total_positive / total * 100,
    ci_lower = binom.confint(total_positive, total, methods = "exact")$lower * 100,
    ci_upper = binom.confint(total_positive, total, methods = "exact")$upper * 100
  ) %>%
  ungroup()

tot.df.neg$tag <- "HH without other DENV IgG+ individuals"

tot.df.fin <- rbind(tot.df.neg, tot.df.pos)

# Plot with ggplot2
plot <- ggplot(tot.df.fin, aes(x = age, y = seroprevalence, col = tag)) +
  geom_point(size = 2) +                         # Dot for each age group
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) + # Error bars
  geom_line(group = 1) +                         # Line connecting the dots
  labs(
    title = NULL,
    x = "Age (years)",
    y = "Seroprevalence (%)"
  ) +
  theme_minimal() +                                              # Clean theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),           # Rotate x-axis labels
    text = element_text(size = 14),
    legend.position = "top"
  )

# ggsave(filename = "seroprevalence_stat_plot.png", plot = plot, width = 8, height = 6, dpi = 300)


#chi square tests
# Counts for single-positive households
single_positive_count <- sum(tot.df.neg$total_positive)
single_negative_count <- sum(tot.df.neg$total_negative)

# Counts for multiple-positive households
multiple_positive_count <- sum(tot.df.pos$total_positive)
multiple_negative_count <- sum(tot.df.pos$total_negative)

# Create a contingency table
contingency_table <- matrix(
  c(single_positive_count, single_negative_count,
    multiple_positive_count, multiple_negative_count),
  nrow = 2, byrow = TRUE
)
colnames(contingency_table) <- c("Positive", "Negative")
rownames(contingency_table) <- c("Single_Positive", "Multiple_Positive")
print(contingency_table)

# Perform chi-square test
chi_test <- chisq.test(contingency_table)

# View the test result
print(chi_test)
