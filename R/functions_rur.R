run_foi_model <- function(model, df, wd, n_ea = 20, n_it = 10000, n_wr = 5000, n_chain = 3,  mean_sigma_n = 0,
                          mean_sigma_v = 0,
                          sd_sigma_n = 0.5,
                          sd_sigma_v = 0.5,
                          wd_out) {
  dir.create(wd_out)
  data_list <- prepare_stan_data(df, n_ea)
  
  saveRDS(data_list, paste0(wd_out,"/data_stan_Timor.rds"))
  
  stan_model <- compile_stan_model(wd, model)
  stan_fit <- fit_stan_model(stan_model, data_list, n_it, n_wr, n_chain, wd_out)
  diagnostics_and_plots(stan_fit, df, wd,wd_out)
  plot_fit(stan_fit, df, wd, wd_out, data_list)
  return(stan_fit)
}

check_foi_model <- function(df, wd, n_ea = 20, n_it = 10000, n_wr = 5000, n_chain = 3,  mean_sigma_n = 0,
                            mean_sigma_v = 0,
                            sd_sigma_n = 0.5,
                            sd_sigma_v = 0.5,
                            wd_out) {
  stan_fit <- readRDS(paste0(wd_out, "fit.rds"))
  diagnostics_and_plots(stan_fit, df, wd,wd_out)
  plot_fit(stan_fit, df, wd, wd_out, data_list)
  return(stan_fit)
}



prepare_stan_data <- function(df, n_ea = 20,  mean_sigma_n = 0,
                              mean_sigma_v = 0,
                              sd_sigma_n = 0.5,
                              sd_sigma_v = 0.5) {
  N <- nrow(df)
  X_v <- df$class
  df$household_id <- as.integer(as.factor(df$household))
  H <- length(unique(df$household_id))
  
  df$village_id <- as.integer(as.factor(df$EA))
  V <- length(unique(df$village_id))
  
  hh_village <- df |>
    dplyr::select(household_id, village_id) |>
    dplyr::distinct() |>
    dplyr::arrange(household_id)
  
  df <- df %>% mutate(
    household_id = as.integer(as.factor(household)),
    village_id   = as.integer(as.factor(EA))
  )
  
  village_class <- df %>%
    group_by(village_id) %>%
    summarise(
      # check unique values in the village
      n_unique = n_distinct(class),
      # majority rule: mean(class) >= 0.5 -> urban (1), else rural (0)
      X_v = as.numeric(mean(class, na.rm = TRUE) >= 0.5),
      .groups = "drop"
    ) %>%
    arrange(village_id)  
  
  X_v_vec <- village_class$X_v   # length V, values 0 or 1
  
  hh <- df$household_id
  vv <- hh_village$village_id
  
  age <- df$age
  y <- df$result
  
  national_mu <- log (0.116114 / (1-0.116114)) #taken from map, logit of the value
  national_sigma <- 0.1
  
  data <- list(
    N = N,
    H = H,
    V = V,
    hh = hh,
    vv = vv,
    age = age,
    y = y,
    national_mu = national_mu,
    national_sigma = national_sigma,
    mean_sigma_n = mean_sigma_n,
    mean_sigma_v = mean_sigma_v,
    sd_sigma_n = sd_sigma_n,
    sd_sigma_v = sd_sigma_v,
    se = 0.77,
    sp = 0.93,
    X_v = X_v_vec
  )
  
  message("Data list prepared for Stan:")
  str(data)
  return(data)
}

compile_stan_model <- function(wd, model) {
  stan_file <- paste0(wd, "/R/", model)
  mod <- rstan::stan_model(stan_file)
  message("Stan model compiled successfully.")
  return(mod)
}

fit_stan_model <- function(mod, data, n_it = 10000, n_wr = 5000, n_chain = 3, wd_out=wd_out) {
  
  stanfit <- rstan::sampling(
    object = mod,
    data = data,
    chains = n_chain,
    iter = n_it,
    warmup = n_wr,
    cores = n_chain,
    refresh = 50
  )
  
  saveRDS(stanfit,paste0(wd_out, "fit.rds"))
  message("Model fitting complete.")
  return(stanfit)
}

diagnostics_and_plots <- function(fit, df, wd, wd_out, mean_sigma_n = 0,
                                  mean_sigma_v = 0,
                                  sd_sigma_n = 0.5,
                                  sd_sigma_v = 0.5) {
  t <- traceplot(fit, pars = c("se", "sp", "beta", "sigma_n", "sigma_v", "delta_v[1]", "delta_v[2]",
                               "delta_v[3]", "delta_h[1]", "delta_h[7]", "delta_h[30]"))
  
  ggsave(paste0(wd_out, "/traceplots.png"), plot = t, width = 10, height = 6)
  
  posterior <- rstan::extract(fit)
  n_chain <- length(fit@sim$samples)
  n_iter <- fit@sim$iter - fit@sim$warmup
  n_draws <- n_chain * n_iter
  V <- dim(posterior$mu_v)[2]
  H <- dim(posterior$lambda_h)[2]
  
  ## FOI
  FOI <- data.frame(matrix(NA,nrow= V, ncol=3))
  for (v in 1:V) FOI[v,1:3] <- quantile(posterior$lambda_v[,v], c(0.5,0.025,0.975))
  FOI$EA <- unique(df$EA) 
  
  f <- ggplot(FOI, aes(x=as.factor(EA), y = X1, ymin =X2, ymax=X3))+
    geom_point()+
    geom_errorbar()+
    geom_hline(yintercept = 1, color = "red", linetype = "dashed", size = 1)+
    xlab("EA") + ylab("FOI") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotazione delle etichette
    )
  
  ggsave(paste0(wd_out, "/FOI.png"), plot = f, width = 10, height = 6)
  message("FOI plots saved.")
  
}

plot_fit <- function(stan_fit, df, wd, wd_out, data_list){
  posterior <- rstan::extract(stan_fit)
  n_chain <- length(stan_fit@sim$samples)
  n_iter <- stan_fit@sim$iter - stan_fit@sim$warmup
  n_draws <- n_chain * n_iter
  V <- dim(posterior$mu_v)[2]
  H <- dim(posterior$lambda_h)[2]
  
  data <- data_list
  
  data_f <- data.frame(
    age = data_list$age,
    y = data_list$y,
    hh = data_list$hh
  )
  
  df <- df %>%
    distinct(household, EA, household_id, village_id)
  
  data_f$village <- data_list$vv[data_f$hh]  
  
  data_f <- data_f %>%
    left_join(
      df %>% select(household_id , village_id, EA, household),
      by = c("hh" = "household_id",
             "village" = "village_id")
    )
  
  
  data_f <- data_f %>%
    mutate(
      age_group = cut(age, breaks = seq(0, max(age) + 10, by = 10), right = FALSE),
      age_mid = as.numeric(sub("\\[(\\d+),.*", "\\1", age_group)) + 5
    )
  
  
  obs_df <- data_f %>%
    group_by(EA, age_group, age_mid) %>%
    summarise(
      obs_prev = mean(y),
      n = n(),
      se = sqrt(obs_prev * (1 - obs_prev) / n),
      lower = pmax(0, obs_prev - 1.96*se),
      upper = pmin(1, obs_prev + 1.96*se)
    ) %>% ungroup()
  
  post <- posterior
  
  # For each iteration, compute mean predicted seropositivity per age group and village
  age_village <- unique(data.frame(age = data$age, village = data$vv[data$hh]))
  
  pred_list <- lapply(1:nrow(post$p_obs), function(iter) {
    p_iter <- post$p_obs[iter, ]
    df_iter <- data.frame(
      village = data$vv[data$hh],
      age = data$age,
      p = p_iter
    )
    df_iter %>%
      group_by(village, age) %>%
      summarise(p_mean = mean(p), .groups = "drop")
  })
  
  pred_df <- bind_rows(pred_list, .id = "iteration") %>%
    group_by(village, age) %>%
    summarise(
      median = median(p_mean),
      lower = quantile(p_mean, 0.025),
      upper = quantile(p_mean, 0.975)
    ) %>% ungroup()
  
  
  pred_df <- pred_df %>%
    left_join(
      df %>% select( village_id, EA),
      by = c("village" = "village_id")
    )
  
  
  # --- 3. Plot by age, faceted by village
  pdf <- ggplot() +
    geom_point(data = obs_df, aes(x = age_mid, y = obs_prev), color = "red", size = 1) +
    geom_errorbar(data = obs_df, aes(x = age_mid, ymin = lower, ymax = upper), color = "red", width = 0.5) +
    geom_ribbon(data = pred_df, aes(x = age, ymin = lower, ymax = upper), fill = "green", alpha = 0.3) +
    geom_line(data = pred_df, aes(x = age, y = median), color = "green", size = 1, alpha = 0.3) +
    facet_wrap(~EA) +
    theme_bw()+
    labs(
      x = "Age",
      y = "Seroprevalence",
      title = NULL
    ) 
  
  ggsave(paste0(wd_out, "/Fit_final.png"), plot = pdf, width = 10, height = 10, scale = 1)
  
}
