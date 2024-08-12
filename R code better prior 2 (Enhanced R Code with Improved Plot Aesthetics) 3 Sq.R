
# Load necessary libraries
library(MCMCpack)  # For Beta distribution functions
library(ggplot2)   # For plotting
library(dplyr)     # For data manipulation
library(tidyr)     # For data reshaping

# Function to simulate treatment outcomes
simulate_treatment_outcomes <- function(theta, n_treatments) {
  rbinom(n_treatments, 1, theta)
}

# Function to perform Bayesian updating sequentially
bayesian_update_sequential <- function(success, n_treatments, alpha_prior, beta_prior) {
  alpha_post <- numeric(n_treatments + 1)
  beta_post <- numeric(n_treatments + 1)
  posterior_mean <- numeric(n_treatments + 1)

  alpha_post[1] <- alpha_prior
  beta_post[1] <- beta_prior
  posterior_mean[1] <- alpha_post[1] / (alpha_post[1] + beta_post[1])

  for (t in 1:n_treatments) {
    alpha_post[t + 1] <- alpha_post[t] + success[t]
    beta_post[t + 1] <- beta_post[t] + 1 - success[t]
    posterior_mean[t + 1] <- alpha_post[t + 1] / (alpha_post[t + 1] + beta_post[t + 1])
  }

  return(data.frame(Treatment = 0:n_treatments, Alpha = alpha_post, Beta = beta_post, Posterior_Mean = posterior_mean))
}

# Simulate a dataset
set.seed(123)
n_patients <- 50
n_treatments <- 20
true_theta <- rbeta(n_patients, 2, 5)  # True success probabilities for each patient
treatment_outcomes <- t(sapply(true_theta, simulate_treatment_outcomes, n_treatments = n_treatments))

# Initialize prior parameters for Beta distribution
alpha_0 <- 1
beta_0 <- 1

# Apply Bayesian updating sequentially for each patient
sequential_posterior_data <- lapply(1:n_patients, function(i) {
  bayesian_update_sequential(treatment_outcomes[i, ], n_treatments, alpha_0, beta_0) %>%
    mutate(Patient = i, True_Theta = true_theta[i])
})

# Combine into a single data frame
df_posterior <- bind_rows(sequential_posterior_data)

# Plot the posterior means over time for a subset of patients
patients_to_plot <- sample(unique(df_posterior$Patient), 5)  # Randomly select 5 patients for visualization
plot_data <- df_posterior %>%
  filter(Patient %in% patients_to_plot)

plot_sequential <- ggplot(plot_data, aes(x = Treatment, y = Posterior_Mean, color = factor(Patient))) +
  geom_line() +
  geom_hline(aes(yintercept = True_Theta), linetype = "dashed") +
  labs(title = "Sequential Posterior Means Over Time for Selected Patients",
       x = "Treatment Number",
       y = "Posterior Mean Success Probability",
       color = "Patient") +
  theme_minimal(base_size = 15)

ggsave("sequential_posterior_means_example.png", plot = plot_sequential, width = 10, height = 8)

# Summary statistics for the final posterior means
summary_fit <- df_posterior %>%
  filter(Treatment == n_treatments) %>%
  summarize(
    Mean_True_Theta = mean(True_Theta),
    Mean_Posterior_Mean = mean(Posterior_Mean),
    SD_Posterior_Mean = sd(Posterior_Mean)
  )

write.csv(summary_fit, "summary_fit_sequential_example.csv", row.names = FALSE)
print(summary_fit)

# Create data frame for posterior parameters
df_outcomes <- as.data.frame(treatment_outcomes)
df_outcomes <- df_outcomes %>%
  mutate(Patient = row_number()) %>%
  pivot_longer(-Patient, names_to = "Treatment", values_to = "Outcome") %>%
  mutate(Treatment = as.numeric(gsub("V", "", Treatment)))

# Calculate observed proportions
df_outcomes <- df_outcomes %>%
  group_by(Patient) %>%
  mutate(Cumulative_Success = cumsum(Outcome),
         Observed_Proportion = Cumulative_Success / Treatment)

# Plot actual dataset outcomes with proportions
plot_actual_data <- ggplot(df_outcomes, aes(x = Treatment, y = Outcome, group = Patient, color = factor(Patient))) +
  geom_point(alpha = 0.7) +
  geom_line(aes(y = Observed_Proportion), linetype = "dotted", size = 1) +
  labs(title = "Actual Treatment Outcomes and Observed Proportions for Each Patient",
       x = "Treatment Number",
       y = "Outcome (0 = Failure, 1 = Success)",
       color = "Patient") +
  theme_minimal(base_size = 15)
ggsave("actual_data_example.png", plot = plot_actual_data, width = 10, height = 8)

# Plot posterior mean and true theta for the first few patients
df_subset <- df_posterior %>% filter(Patient %in% patients_to_plot)

plot_posterior <- ggplot(df_subset, aes(x = Treatment)) +
  geom_line(aes(y = Posterior_Mean, color = "Posterior Mean"), size = 1.2) +
  geom_line(aes(y = True_Theta, color = "True Theta"), linetype = "dashed", size = 1.2) +
  facet_wrap(~ Patient, scales = "free_y") +
  labs(title = "Posterior Mean and True Theta for Selected Patients",
       y = "Value",
       color = "Legend") +
  theme_minimal(base_size = 15)
ggsave("posterior_mean_true_theta_example.png", plot = plot_posterior, width = 10, height = 8)

# Plot overall distribution of posterior means vs true theta for the last treatment
df_final <- df_posterior %>% filter(Treatment == n_treatments)

plot_comparison <- ggplot(df_final, aes(x = True_Theta, y = Posterior_Mean)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Comparison of Posterior Mean and True Theta",
       x = "True Theta",
       y = "Posterior Mean") +
  theme_minimal(base_size = 15)
ggsave("comparison_posterior_true_theta_example.png", plot = plot_comparison, width = 10, height = 8)

# Create and save summary table for the first few patients
summary_table <- df_posterior %>%
  filter(Patient %in% patients_to_plot, Treatment == n_treatments) %>%
  select(Patient, Alpha, Beta, True_Theta, Posterior_Mean)
write.csv(summary_table, "summary_table_example.csv", row.names = FALSE)

# Reliability Checking: Calculate mean squared error (MSE) between posterior mean and true theta
mse_results <- data.frame(Dataset = 1, MSE = mean((df_final$Posterior_Mean - df_final$True_Theta)^2))

print(mse_results)
write.csv(mse_results, "mse_results_example.csv", row.names = FALSE)

