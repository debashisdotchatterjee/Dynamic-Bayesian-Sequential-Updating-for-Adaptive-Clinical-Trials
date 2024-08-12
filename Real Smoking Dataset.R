# Install and load required packages
if (!requireNamespace("HSAUR3", quietly = TRUE)) {
  install.packages("HSAUR3")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
if (!requireNamespace("MCMCpack", quietly = TRUE)) {
  install.packages("MCMCpack")
}

library(HSAUR3)
library(dplyr)
library(ggplot2)
library(tidyr)
library(MCMCpack)

# Load the data
data("smoking", package = "HSAUR3")

# Display the structure and first few rows of the dataset
print(head(smoking))

# Prior parameters for Beta distribution (assuming non-informative priors for simplicity)
alpha_prior <- 1
beta_prior <- 1

# Function to perform Bayesian updating
bayesian_update <- function(success, total, alpha_prior, beta_prior) {
  alpha_post <- alpha_prior + success
  beta_post <- beta_prior + total - success
  return(c(alpha_post, beta_post))
}

# Apply Bayesian updating for each row in the dataset
smoking <- smoking %>%
  rowwise() %>%
  mutate(
    alpha_post_treatment = bayesian_update(qt, tt, alpha_prior, beta_prior)[1],
    beta_post_treatment = bayesian_update(qt, tt, alpha_prior, beta_prior)[2],
    alpha_post_control = bayesian_update(qc, tc, alpha_prior, beta_prior)[1],
    beta_post_control = bayesian_update(qc, tc, alpha_prior, beta_prior)[2]
  )

# Calculate posterior means
smoking <- smoking %>%
  mutate(
    posterior_mean_treatment = alpha_post_treatment / (alpha_post_treatment + beta_post_treatment),
    posterior_mean_control = alpha_post_control / (alpha_post_control + beta_post_control),
    observed_proportion_treatment = qt / tt,
    observed_proportion_control = qc / tc
  )

# Posterior distributions using Beta distribution
posterior_samples_treatment <- rbeta(10000, mean(smoking$alpha_post_treatment), mean(smoking$beta_post_treatment))
posterior_samples_control <- rbeta(10000, mean(smoking$alpha_post_control), mean(smoking$beta_post_control))

# Plot posterior distributions
plot_data <- data.frame(
  Treatment = posterior_samples_treatment,
  Control = posterior_samples_control
)

df_melted <- gather(plot_data, key = "Group", value = "Posterior_Sample")

plot_posterior <- ggplot(df_melted, aes(x = Posterior_Sample, fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Posterior Distributions of Success Probabilities",
       x = "Success Probability",
       y = "Density",
       fill = "Group") +
  theme_minimal(base_size = 15)

ggsave("posterior_distributions_smoke.png", plot = plot_posterior, width = 10, height = 8)
plot_posterior
# Summary statistics for model fit
summary_fit <- smoking %>%
  summarize(
    Mean_Observed_Proportion_Treatment = mean(observed_proportion_treatment),
    Mean_Observed_Proportion_Control = mean(observed_proportion_control),
    Mean_Posterior_Mean_Treatment = mean(posterior_mean_treatment),
    Mean_Posterior_Mean_Control = mean(posterior_mean_control)
  )

write.csv(summary_fit, "summary_fit.csv", row.names = FALSE)
print(summary_fit)

# Histograms of posterior means for treatment and control groups
hist_treatment <- ggplot(smoking, aes(x = posterior_mean_treatment)) +
  geom_histogram(binwidth = 0.05, fill = "blue", alpha = 0.7) +
  labs(title = "Histogram of Posterior Means for Treatment Group",
       x = "Posterior Mean Success Probability",
       y = "Frequency") +
  theme_minimal(base_size = 15)

ggsave("histogram_treatment.png", plot = hist_treatment, width = 10, height = 8)

hist_control <- ggplot(smoking, aes(x = posterior_mean_control)) +
  geom_histogram(binwidth = 0.05, fill = "red", alpha = 0.7) +
  labs(title = "Histogram of Posterior Means for Control Group",
       x = "Posterior Mean Success Probability",
       y = "Frequency") +
  theme_minimal(base_size = 15)

ggsave("histogram_control.png", plot = hist_control, width = 10, height = 8)
