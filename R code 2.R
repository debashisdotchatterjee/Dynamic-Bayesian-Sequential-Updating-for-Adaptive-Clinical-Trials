# Load necessary libraries
library(MCMCpack)  # For Beta distribution functions
library(ggplot2)   # For plotting
library(dplyr)     # For data manipulation
library(tidyr)     # For data reshaping

# Create directories for saving plots if they don't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Function to save plot
save_plot <- function(plot, filename) {
  ggsave(filename = filename, plot = plot, path = "plots", width = 8, height = 6)
}

# Simulate multiple datasets
set.seed(123)
n_datasets <- 5
datasets <- list()

for (i in 1:n_datasets) {
  n_patients <- 100
  n_treatments <- 5
  true_theta <- rbeta(n_patients, 2, 5)  # True success probabilities for each patient

  simulate_treatment_outcomes <- function(theta, n_treatments) {
    rbinom(n_treatments, 1, theta)
  }

  treatment_outcomes <- t(sapply(true_theta, simulate_treatment_outcomes, n_treatments = n_treatments))

  df_outcomes <- as.data.frame(treatment_outcomes)
  df_outcomes$Patient <- 1:n_patients
  df_outcomes <- gather(df_outcomes, key = "Treatment", value = "Outcome", -Patient)
  df_outcomes$Treatment <- as.numeric(gsub("V", "", df_outcomes$Treatment))
  df_outcomes$True_Theta <- rep(true_theta, each = n_treatments)

  datasets[[i]] <- df_outcomes
}

# Plot raw data for all datasets
for (i in 1:n_datasets) {
  df_outcomes <- datasets[[i]]
  p <- ggplot(df_outcomes, aes(x = Treatment, y = Outcome, group = Patient)) +
    geom_line(alpha = 0.1) +
    geom_point(alpha = 0.1) +
    labs(title = paste("Raw Treatment Outcomes for Each Patient - Dataset", i),
         x = "Treatment Number",
         y = "Outcome (0 = Failure, 1 = Success)") +
    theme_minimal()

  save_plot(p, paste("raw_data_dataset_", i, ".png", sep = ""))
}

# Bayesian updating and plotting for each dataset
for (i in 1:n_datasets) {
  df_outcomes <- datasets[[i]]

  # Initialize prior parameters for Beta distribution
  alpha_0 <- 1
  beta_0 <- 1

  # Bayesian updating function
  bayesian_update <- function(prior_alpha, prior_beta, outcome) {
    if (outcome == 1) {
      prior_alpha <- prior_alpha + 1
    } else {
      prior_beta <- prior_beta + 1
    }
    return(c(prior_alpha, prior_beta))
  }

  # Store the updated parameters
  n_patients <- length(unique(df_outcomes$Patient))
  n_treatments <- length(unique(df_outcomes$Treatment))

  alpha_post <- matrix(NA, nrow = n_patients, ncol = n_treatments + 1)
  beta_post <- matrix(NA, nrow = n_patients, ncol = n_treatments + 1)
  alpha_post[, 1] <- alpha_0
  beta_post[, 1] <- beta_0

  # Sequential treatment and Bayesian updating
  for (p in 1:n_patients) {
    for (t in 1:n_treatments) {
      outcome <- df_outcomes %>% filter(Patient == p & Treatment == t) %>% select(Outcome) %>% as.numeric()
      alpha_beta <- bayesian_update(alpha_post[p, t], beta_post[p, t], outcome)
      alpha_post[p, t + 1] <- alpha_beta[1]
      beta_post[p, t + 1] <- alpha_beta[2]
    }
  }

  # Create a data frame for posterior parameters
  df_posterior <- data.frame(
    Patient = rep(1:n_patients, each = n_treatments + 1),
    Treatment = rep(0:n_treatments, times = n_patients),
    Alpha = as.vector(t(alpha_post)),
    Beta = as.vector(t(beta_post)),
    True_Theta = rep(unique(df_outcomes$True_Theta), each = n_treatments + 1)
  )

  # Calculate the posterior mean (expected success probability)
  df_posterior <- df_posterior %>%
    mutate(Posterior_Mean = Alpha / (Alpha + Beta))

  # Plotting the posterior mean and true theta for the first few patients
  patients_to_plot <- 1:10
  df_subset <- df_posterior %>% filter(Patient %in% patients_to_plot)

  p <- ggplot(df_subset, aes(x = Treatment)) +
    geom_line(aes(y = Posterior_Mean, color = "Posterior Mean")) +
    geom_line(aes(y = True_Theta, color = "True Theta"), linetype = "dashed") +
    facet_wrap(~ Patient, scales = "free_y") +
    labs(title = paste("Posterior Mean and True Theta for Selected Patients - Dataset", i),
         y = "Value",
         color = "Legend") +
    theme_minimal()

  save_plot(p, paste("posterior_true_theta_dataset_", i, ".png", sep = ""))

  # Plot overall distribution of posterior means vs true theta for the last treatment
  df_final <- df_posterior %>% filter(Treatment == n_treatments)

  p <- ggplot(df_final, aes(x = True_Theta, y = Posterior_Mean)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = paste("Comparison of Posterior Mean and True Theta - Dataset", i),
         x = "True Theta",
         y = "Posterior Mean") +
    theme_minimal()

  save_plot(p, paste("posterior_vs_true_theta_dataset_", i, ".png", sep = ""))

  # Create a summary table for the first few patients
  summary_table <- df_posterior %>%
    filter(Patient %in% patients_to_plot, Treatment == n_treatments) %>%
    select(Patient, Alpha, Beta, True_Theta, Posterior_Mean)

  write.csv(summary_table, file = paste("summary_table_dataset_", i, ".csv", sep = ""), row.names = FALSE)

  # Plotting the results for the first patient
  patient_id <- 1
  df_patient <- df_posterior %>% filter(Patient == patient_id)

  p <- ggplot(df_patient, aes(x = Treatment)) +
    geom_line(aes(y = Alpha, color = "Alpha")) +
    geom_line(aes(y = Beta, color = "Beta")) +
    labs(title = paste("Posterior Parameters for Patient", patient_id, "- Dataset", i),
         y = "Parameter Value",
         color = "Parameter") +
    theme_minimal()

  save_plot(p, paste("posterior_parameters_patient_", patient_id, "_dataset_", i, ".png", sep = ""))

  # Calculate the posterior mean (expected success probability) for the first patient
  posterior_mean <- df_patient$Posterior_Mean
  print(posterior_mean)

  # Decision-making: Select the treatment with the highest expected success probability for the next patient
  # Expected success probability after each treatment
  expected_success <- numeric(n_treatments)
  for (t in 1:n_treatments) {
    expected_success[t] <- alpha_post[patient_id, t] / (alpha_post[patient_id, t] + beta_post[patient_id, t])
  }

  print(expected_success)

  # The treatment with the highest expected success probability
  best_treatment <- which.max(expected_success)
  print(best_treatment)
}

