# Set seed for Reproducability.
set.seed(8)

# Create the samples of calls balanced by species and individual.
N <- 1000
resampled_calls <- balanced_samples(smooth_and_regularise_call_densities(), n_samples = N)

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
cumsum(mean_proportion_variance)[1:10]

sum(sample_proportion_variance[7,] > 0.03)
proportion_mask <- (sample_proportion_variance > 0.03)
sum(proportion_mask)
# The centering / mean curves for the samples.
get_centring_curve <- function(df_col){
  single_sample <- matrix(unlist(df_col), nrow = length(df_col), ncol = length(df_col[[1]]), byrow = T)
  pca <- prcomp(single_sample, scale. = T)

  return(pca$center)
}

sample_centring_curve <- apply(resampled_calls, 2, get_centring_curve)
matplot(sample_centring_curve, type = 'l')
mean_centring_curve <- apply(sample_centring_curve, 1, mean)
lines(mean_centring_curve, lwd = 3)

# Get 10 most descriptive components.
get_components <-function(df_col, n_components = 10){
  single_sample <- matrix(unlist(df_col), nrow = length(df_col), ncol = length(df_col[[1]]), byrow = T)
  pca <- prcomp(single_sample, scale. = T)

  return(pca$rotation[, 1:n_components])
}

examined_components <- 7
sample_components <- apply(resampled_calls, 2, get_components, n_components = examined_components)
sample_components <- array(sample_components, dim = c(208, examined_components, N))

proportion_mask <- proportion_mask[1:examined_components, ]
sum(proportion_mask)

# Identify component peak frequencies
peak_frequency <- apply(abs(sample_components), c(2,3), which.max)
boxplot_components <- boxplot(t(peak_frequency), xlab = 'Principal Components', ylab = 'Frequency', main = 'Boxplot of Max Asolute value Position')
boxplot_components$stats




# Orientate components so absolute maximum is positive
absmax <- function(x) { x[which.max( abs(x) )]}

# first component
matplot(sample_components[,1,], type = 'l')
first_sample <-  apply(sample_components[, 1, ], 2, function(x) sign(absmax(x))*x)
matplot(first_sample, type = 'l')

first_mean <- apply(first_sample, 1, mean)
lines(first_mean, lwd = 3)

# second component
matplot(sample_components[,2,], type = 'l')
second_sample <- apply(sample_components[, 2, ], 2, function(x) sign(absmax(x[50:100]))*x)

matplot(second_sample, type = 'l')

second_mean <- apply(second_sample, 1, mean)
lines(second_mean, lwd = 3)

# subsequent components
identify <- function(components = sample_components[, -c(1,2),], max_range = c(41, 45)){
  peaks <- apply(abs(components), c(2,3), which.max)

  mask <- peaks >= max_range[1] & peaks <= max_range[2]

  mask <- array(rep(mask, dim(components)[1]), dim = c(dim(components)[2:3], dim(components)[1]))
  mask <- aperm(mask, c(3,1,2))

  curves <- apply(components, c(2,3), function(x) sign(absmax(x))*x)
  curves <- curves[mask]
  curves <- matrix(curves, ncol = dim(components)[1], byrow = T)

  return(t(curves))
}

# third
third_sample <- identify(max_range = c(41,45))
matplot((third_sample), type = 'l')

third_mean <- apply(third_sample, 1, mean)
lines(third_mean, lwd = 3)

# fourth
fourth_sample <- identify(max_range = c(30,34))
matplot((fourth_sample), type = 'l')

fourth_mean <- apply(fourth_sample, 1, mean)
lines(fourth_mean, lwd = 3)

# fifth
fifth_sample <- identify(max_range = c(15,19))
matplot((fifth_sample), type = 'l')

fifth_mean <- apply(fifth_sample, 1, mean)
lines(fifth_mean, lwd = 3)

# sixth
sixth_sample <- identify(max_range = c(24, 28))
matplot((sixth_sample), type = 'l')

sixth_mean <- apply(sixth_sample, 1, mean)
lines(sixth_mean, lwd = 3)

components <- cbind(first_mean, second_mean, third_mean, fourth_mean, fifth_mean, sixth_mean)

principal_components <- apply(components, 1, '/', sqrt(diag(t(components)%*%(components))))
image(t(principal_components)%*%principal_components)
image(((principal_components)%*%t(principal_components)))
matplot(t(principal_components), type = 'l', lty = 1)

getwd()

spectral_fpca <- t(principal_components)
spectral_mean <- mean_centring_curve


devtools::use_data(spectral_fpca)
devtools::use_data(spectral_mean)
