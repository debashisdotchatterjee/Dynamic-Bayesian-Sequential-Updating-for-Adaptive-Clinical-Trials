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

# Simulate a dataset
set.seed(123)
n_patients <- 100
n_treatments <- 100  # Increased number of treatments to simulate more pegs
true_theta <- rbeta(n_patients, 2, 5)  # True success probabilities for each patient
treatment_outcomes <- t(sapply(true_theta, simulate_treatment_outcomes, n_treatments = n_treatments))

# Initialize prior parameters for Beta distribution
alpha_0 <- 1
beta_0 <- 6

# Perform Bayesian updating
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
results <- list(alpha_post = alpha_post, beta_post = beta_post, posterior_mean = posterior_mean, true_theta = true_theta, treatment_outcomes = treatment_outcomes)

# Create data frame for posterior parameters
df_posterior <- data.frame(
  Patient = rep(1:n_patients, each = n_treatments + 1),
  Treatment = rep(0:n_treatments, times = n_patients),
  Alpha = as.vector(t(results$alpha_post)),
  Beta = as.vector(t(results$beta_post)),
  Posterior_Mean = as.vector(t(results$posterior_mean)),
  True_Theta = rep(results$true_theta, each = n_treatments + 1)
)

# Plot actual dataset outcomes with proportions
df_outcomes <- as.data.frame(results$treatment_outcomes)
df_outcomes <- df_outcomes %>%
  mutate(Patient = row_number()) %>%
  pivot_longer(-Patient, names_to = "Treatment", values_to = "Outcome") %>%
  mutate(Treatment = as.numeric(gsub("V", "", Treatment)))

# Calculate observed proportions
df_outcomes <- df_outcomes %>%
  group_by(Patient) %>%
  mutate(Cumulative_Success = cumsum(Outcome),
         Observed_Proportion = Cumulative_Success / Treatment)

plot_actual_data <- ggplot(df_outcomes, aes(x = Treatment, y = Outcome, group = Patient, color = factor(Patient))) +
 # geom_line(alpha = 0.7) +
  geom_point(alpha = 0.7) +
  geom_line(aes(y = Observed_Proportion), linetype = "dotted", size = 1) +
  labs(title = "Actual Treatment Outcomes and Observed Proportions for Each Patient",
       x = "Treatment Number",
       y = "Outcome (0 = Failure, 1 = Success)",
       color = "Patient") +
  theme_minimal(base_size = 15)
ggsave("actual_data.png", plot = plot_actual_data, width = 10, height = 8)

# Plot posterior mean and true theta for the first few patients
patients_to_plot <- 1:12
df_subset <- df_posterior %>% filter(Patient %in% patients_to_plot)

plot_posterior <- ggplot(df_subset, aes(x = Treatment)) +
  geom_line(aes(y = Posterior_Mean, color = "Posterior Mean"), size = 1.2) +
  geom_line(aes(y = True_Theta, color = "True Theta"), linetype = "dashed", size = 1.2) +
  facet_wrap(~ Patient, scales = "free_y") +
  labs(title = "Posterior Mean and True Theta for Selected Patients",
       y = "Value",
       color = "Legend") +
  theme_minimal(base_size = 15)
ggsave("posterior_mean_true_theta.png", plot = plot_posterior, width = 10, height = 8)

# Plot overall distribution of posterior means vs true theta for the last treatment
df_final <- df_posterior %>% filter(Treatment == n_treatments)

plot_comparison <- ggplot(df_final, aes(x = True_Theta, y = Posterior_Mean)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Comparison of Posterior Mean and True Theta",
       x = "True Theta",
       y = "Posterior Mean") +
  theme_minimal(base_size = 15)
ggsave("comparison_posterior_true_theta.png", plot = plot_comparison, width = 10, height = 8)

# Create and save summary table for the first few patients
summary_table <- df_posterior %>%
  filter(Patient %in% patients_to_plot, Treatment == n_treatments) %>%
  select(Patient, Alpha, Beta, True_Theta, Posterior_Mean)
write.csv(summary_table, "summary_table.csv", row.names = FALSE)

# Reliability Checking: Calculate mean squared error (MSE) between posterior mean and true theta
mse_results <- data.frame(Dataset = 1, MSE = mean((df_final$Posterior_Mean - df_final$True_Theta)^2))

print(mse_results)
write.csv(mse_results, "mse_results.csv", row.names = FALSE)


##############################

# Accuracy Checking: Plot observed vs. predicted success probabilities

# Calculate observed success proportions
observed_success <- df_outcomes %>%
  group_by(Patient) %>%
  summarize(Observed_Proportion = sum(Outcome) / n_treatments)

# Combine with posterior means for the last treatment
df_accuracy <- df_posterior %>%
  filter(Treatment == n_treatments) %>%
  select(Patient, Posterior_Mean, True_Theta) %>%
  left_join(observed_success, by = "Patient")

# Plot observed vs. predicted success probabilities
plot_accuracy <- ggplot(df_accuracy, aes(x = Observed_Proportion, y = Posterior_Mean)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Observed vs. Predicted Success Probabilities",
       x = "Observed Success Proportion",
       y = "Predicted Success Proportion (Posterior Mean)") +
  theme_minimal(base_size = 15)
ggsave("accuracy_plot.png", plot = plot_accuracy, width = 10, height = 8)

# Create summary table for model fit
summary_fit <- df_accuracy %>%
  summarize(Mean_Observed_Proportion = mean(Observed_Proportion),
            Mean_Predicted_Proportion = mean(Posterior_Mean),
            Correlation = cor(Observed_Proportion, Posterior_Mean))

write.csv(summary_fit, "summary_fit.csv", row.names = FALSE)

########################

# Load necessary libraries if not already loaded
library(ggplot2)
library(dplyr)
library(tidyr)

# Calculate observed success proportions
observed_success <- df_outcomes %>%
  group_by(Patient) %>%
  summarize(Observed_Proportion = sum(Outcome) / n_treatments, .groups = 'drop')

# Combine with posterior means for the last treatment
df_accuracy <- df_posterior %>%
  filter(Treatment == n_treatments) %>%
  select(Patient, Posterior_Mean, True_Theta) %>%
  left_join(observed_success, by = "Patient")

# Plot observed vs. predicted success probabilities
plot_accuracy <- ggplot(df_accuracy, aes(x = Observed_Proportion, y = Posterior_Mean)) +
  geom_point(alpha = 0.8, color = "forestgreen", size = 3, show.legend = TRUE, aes(shape = "Observed")) + # Increase point size and add shape legend
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1.2, show.legend = TRUE) + # Adjust line width and add to legend
  labs(title = "Observed vs. Predicted Success Probabilities",
       subtitle = "Comparison of each patient's observed success proportion against the model's predictions.",
       x = "Observed Success Proportion",
       y = "Predicted Success Proportion (Posterior Mean)",
       caption = "Data points show individual patient outcomes; red dashed line indicates perfect prediction.") +
  theme_minimal(base_size = 15) +
  theme(legend.title = element_blank(), # Removing the legend title
        legend.position = "bottom", # Move legend to bottom
        plot.title = element_text(size = 16, face = "bold"), # Enhance title
        plot.subtitle = element_text(size = 14), # Enhance subtitle
        plot.caption = element_text(size = 12)) +
  scale_shape_manual(values = 16) # Adjust shape to solid circles

# Save the plot
ggsave("accuracy_plot.png", plot = plot_accuracy, width = 10, height = 8, dpi = 300)
