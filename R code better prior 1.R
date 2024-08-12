# Load necessary libraries
library(MCMCpack)  # For Beta distribution functions
library(ggplot2)   # For plotting
library(dplyr)     # For data manipulation
library(tidyr)     # For data reshaping

# Function to simulate treatment outcomes
simulate_treatment_outcomes <- function(theta, n_treatments) {
  rbinom(n_treatments, 1, theta)
}

# Function to perform Bayesian updating
bayesian_update <- function(prior_alpha, prior_beta, outcome) {
  if (outcome == 1) {
    prior_alpha <- prior_alpha + 1
  } else {
    prior_beta <- prior_beta + 1
  }
  return(c(prior_alpha, prior_beta))
}

# Simulate datasets
set.seed(123)
n_patients <- 100
n_treatments <- 5
n_datasets <- 3  # Number of datasets to simulate

all_datasets <- list()
for (d in 1:n_datasets) {
  true_theta <- rbeta(n_patients, 2, 5)  # True success probabilities for each patient
  treatment_outcomes <- t(sapply(true_theta, simulate_treatment_outcomes, n_treatments = n_treatments))
  all_datasets[[d]] <- list(true_theta = true_theta, treatment_outcomes = treatment_outcomes)
}

# Initialize prior parameters for Beta distribution
alpha_0 <- 2
beta_0 <- 5

# Loop through each dataset and perform Bayesian updating
results <- list()
for (d in 1:n_datasets) {
  dataset <- all_datasets[[d]]
  true_theta <- dataset$true_theta
  treatment_outcomes <- dataset$treatment_outcomes

  alpha_post <- matrix(NA, nrow = n_patients, ncol = n_treatments + 1)
  beta_post <- matrix(NA, nrow = n_patients, ncol = n_treatments + 1)
  alpha_post[, 1] <- alpha_0
  beta_post[, 1] <- beta_0

  for (i in 1:n_patients) {
    for (t in 1:n_treatments) {
      outcome <- treatment_outcomes[i, t]
      alpha_beta <- bayesian_update(alpha_post[i, t], beta_post[i, t], outcome)
      alpha_post[i, t + 1] <- alpha_beta[1]
      beta_post[i, t + 1] <- alpha_beta[2]
    }
  }

  posterior_mean <- alpha_post / (alpha_post + beta_post)
  results[[d]] <- list(alpha_post = alpha_post, beta_post = beta_post, posterior_mean = posterior_mean, true_theta = true_theta, treatment_outcomes = treatment_outcomes)
}

# Function to save plots
save_plot <- function(plot, filename) {
  ggsave(filename, plot = plot, width = 10, height = 8)
}

# Plot and save results
for (d in 1:n_datasets) {
  dataset_results <- results[[d]]

  # Create data frame for posterior parameters
  df_posterior <- data.frame(
    Patient = rep(1:n_patients, each = n_treatments + 1),
    Treatment = rep(0:n_treatments, times = n_patients),
    Alpha = as.vector(t(dataset_results$alpha_post)),
    Beta = as.vector(t(dataset_results$beta_post)),
    Posterior_Mean = as.vector(t(dataset_results$posterior_mean)),
    True_Theta = rep(dataset_results$true_theta, each = n_treatments + 1)
  )

  # Plot raw data
  df_outcomes <- as.data.frame(dataset_results$treatment_outcomes)
  df_outcomes <- df_outcomes %>%
    mutate(Patient = row_number()) %>%
    pivot_longer(-Patient, names_to = "Treatment", values_to = "Outcome") %>%
    mutate(Treatment = as.numeric(gsub("V", "", Treatment)))

  plot_raw_data <- ggplot(df_outcomes, aes(x = Treatment, y = Outcome, group = Patient)) +
    geom_line(alpha = 0.1) +
    geom_point(alpha = 0.1) +
    labs(title = paste("Raw Treatment Outcomes for Each Patient (Dataset", d, ")"),
         x = "Treatment Number",
         y = "Outcome (0 = Failure, 1 = Success)") +
    theme_minimal()
  save_plot(plot_raw_data, paste0("raw_data_dataset_", d, ".png"))

  # Plot posterior mean and true theta for the first few patients
  patients_to_plot <- 1:10
  df_subset <- df_posterior %>% filter(Patient %in% patients_to_plot)

  plot_posterior <- ggplot(df_subset, aes(x = Treatment)) +
    geom_line(aes(y = Posterior_Mean, color = "Posterior Mean")) +
    geom_line(aes(y = True_Theta, color = "True Theta"), linetype = "dashed") +
    facet_wrap(~ Patient, scales = "free_y") +
    labs(title = paste("Posterior Mean and True Theta for Selected Patients (Dataset", d, ")"),
         y = "Value",
         color = "Legend") +
    theme_minimal()
  save_plot(plot_posterior, paste0("posterior_mean_true_theta_dataset_", d, ".png"))

  # Plot overall distribution of posterior means vs true theta for the last treatment
  df_final <- df_posterior %>% filter(Treatment == n_treatments)

  plot_comparison <- ggplot(df_final, aes(x = True_Theta, y = Posterior_Mean)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = paste("Comparison of Posterior Mean and True Theta (Dataset", d, ")"),
         x = "True Theta",
         y = "Posterior Mean") +
    theme_minimal()
  save_plot(plot_comparison, paste0("comparison_posterior_true_theta_dataset_", d, ".png"))

  # Create and save summary table for the first few patients
  summary_table <- df_posterior %>%
    filter(Patient %in% patients_to_plot, Treatment == n_treatments) %>%
    select(Patient, Alpha, Beta, True_Theta, Posterior_Mean)
  write.csv(summary_table, paste0("summary_table_dataset_", d, ".csv"), row.names = FALSE)
}

# Reliability Checking: Calculate mean squared error (MSE) between posterior mean and true theta
mse_results <- data.frame(Dataset = integer(), MSE = numeric())
for (d in 1:n_datasets) {
  dataset_results <- results[[d]]
  df_final <- dataset_results$posterior_mean[, n_treatments + 1]
  true_theta <- dataset_results$true_theta
  mse <- mean((df_final - true_theta)^2)
  mse_results <- rbind(mse_results, data.frame(Dataset = d, MSE = mse))
}

print(mse_results)
write.csv(mse_results, "mse_results.csv", row.names = FALSE)

