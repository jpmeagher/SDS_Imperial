# Set seed for Reproducability.
set.seed(8)

# Create the samples of calls balanced by species and individual.
N <- 1000
resampled_calls <- balanced_samples(rescale_smooth_densities(), n_samples = N)

# Proportion of variance accounted for by each component.
get_proportion_variance <- function(df_col){
  single_sample <- matrix(unlist(df_col), nrow = length(df_col), ncol = length(df_col[[1]]), byrow = T)
  pca <- prcomp(single_sample, scale. = T)
  proportion_variance <- pca$sdev^2 /(sum(pca$sdev^2))
  return(proportion_variance)
}

sample_proportion_variance <- apply(resampled_calls, 2, get_proportion_variance)
matplot(sample_proportion_variance[1:20, 1:100])

mean_proportion_variance <- apply(sample_proportion_variance, 1, mean)
plot(mean_proportion_variance)
mean_proportion_variance[1:10]
plot(cumsum(mean_proportion_variance))

sum(sample_proportion_variance[9,] > 0.03)
proportion_mask <- (sample_proportion_variance > 0.03)
sum(proportion_mask)
# The centering / mean curves for the samples.
get_centring_curve <- function(df_col){
  single_sample <- matrix(unlist(df_col), nrow = length(df_col), ncol = length(df_col[[1]]), byrow = T)
  pca <- prcomp(single_sample, scale. = T)

  return(pca$center)
}

sample_centring_curve <- apply(resampled_calls, 2, get_centring_curve)
matplot(sample_centring_curve[, 1:100], type = 'l')
mean_centring_curve <- apply(sample_centring_curve, 1, mean)
plot(mean_centring_curve, type = 'l')

# Get 10 most descriptive components.
get_components <-function(df_col, n_components = 10){
  single_sample <- matrix(unlist(df_col), nrow = length(df_col), ncol = length(df_col[[1]]), byrow = T)
  pca <- prcomp(single_sample, scale. = T)

  return(pca$rotation[, 1:n_components])
}

examined_components <- 10
sample_components <- apply(resampled_calls, 2, get_components, n_components = examined components)
sample_components <- array(sample_components, dim = c(208, 10, N))
matplot(sample_components[,6,], type = 'l')
proportion_mask <- proportion_mask[1:10, ]
sum(proportion_mask)

# Identify component peak frequencies
peak_frequency <- apply(abs(sample_components), c(2,3), which.max)
boxplot_components <- boxplot(t(peak_frequency), xlab = 'Principal Components', ylab = 'Frequency', main = 'Boxplot of Max Asolute value Position')
boxplot_components$stats
matplot(abs(sample_components[,6,]), type = 'l')


# Select components based on peak frequency
# First Component
absmax <- function(x) { x[which.max( abs(x) )]}

sample_first_component <- apply(sample_components[, 1, ], 2, function(x) sign(absmax(x))*x)


matplot(sample_first_component, type = 'l')
lines(mean_first_component, lwd = 3)


# Second Component
matplot(sample_components[,3,], type = 'l')

sample_second_component <- apply(sample_components[, 2, ], 2, function(x) sign(absmax(x[50:100]))*x)


matplot(sample_second_component, type = 'l')
lines(mean_second_component, lwd = 3)

which.min(sample_components[, 2, 1])

# Identify curves where components have switched order and switch them back.
mask <- apply(sample_first_component, 2, function(x) absmax(x[150:200]) < 0)
matplot(sample_second_component[, mask], type = 'l')

final_first_component <- cbind(sample_first_component[, !mask], sample_second_component[, mask] )
final_second_component <- cbind(sample_second_component[, !mask], sample_first_component[, mask] )

mean_second_component <- apply(final_second_component, 1, mean)
mean_first_component <- apply(final_first_component, 1, mean)
)

# Identify next 5 components
# Peaks under 50
res_proportion_mask <- proportion_mask[3:10,]
res_proportion_mask <- array(rep(res_proportion_mask, 208), c(8, 1000, 208))
identical(res_proportion_mask[1,,], res_proportion_mask[2,,])
res_proportion_mask <- aperm(res_proportion_mask, c(3,1,2))
res_components <- sample_components[, 3:10,]
res_components <- res_components[res_proportion_mask]
plot(res_components[4*208+(1:208)])
dim(res_components)
dim(res_proportion_mask)
# Peak frequency at 42 to 45
peak_mask <- (peak_frequency == 42 | peak_frequency == 43 | peak_frequency == 44)
sum(peak_mask)
matplot(sample_components[, ])


