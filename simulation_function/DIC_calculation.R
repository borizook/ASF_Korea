##################################################
# François, O., & Laval, G. (2011). Deviance information criteria for model selection in approximate Bayesian computation. Statistical Applications in Genetics and Molecular Biology, 10(1).

# Define the Deviance Information Criteria (DIC) function for use in Approximate Bayesian Computation (ABC)
# This function evaluates the fit of statistical models by computing the DIC using various kernel density estimations.
#
# Arguments:
#   - target: A numeric vector or list of observed summary statistics.
#   - postsumstat_from_posterior: A matrix or data frame of simulated summary statistics from the posterior distribution.
#   - postsumstat_from_point: A matrix, data frame, or vector of simulated summary statistics from a point estimate.
#   - kernel: A character string specifying the kernel density estimation method. Default options include "gaussian",
#             "epanechnikov", "rectangular", "triangular", "biweight", "cosine".
#
# Returns:
#   A list containing the expected deviance, pD, DIC, and adjusted distances based on the specified kernel.

DIC_ABC <-function (target, postsumstat_from_posterior, postsumstat_from_point ,kernel = c("gaussian", "epanechnikov", "rectangular", 
                                                                                                     "triangular", "biweight", "cosine")) 
{
  # Function to normalize data
  normalise <- function(x, y) {
    if (mad(y) == 0)
      return(x)
    else
      return(x / mad(y))
  }
  
  # Match the kernel argument with the available choices
  kernel <- match.arg(arg = kernel)
  
  # Convert data.frames to matrices
  if (is.data.frame(postsumstat_from_posterior)) {
    postsumstat_from_posterior <- as.matrix(postsumstat_from_posterior)
  }
  
  if (is.data.frame(postsumstat_from_point)) {
    postsumstat_from_point <- as.matrix(postsumstat_from_point)
    v_postsumstat_from_point <- apply(postsumstat_from_point, MARGIN = 2, FUN = median)
  } else {
    v_postsumstat_from_point <- unlist(postsumstat_from_point)
  }
  
  # Unlist target if it is a list
  if (is.list(target)) {
    target <- unlist(target)
  }
  
  # Ensure the summary statistics are matrices
  if (is.vector(postsumstat_from_posterior)) {
    postsumstat_from_posterior <- matrix(postsumstat_from_posterior, ncol = 1)
  }
  
  if (is.vector(postsumstat_from_point)) {
    postsumstat_from_point <- matrix(postsumstat_from_point, ncol = 1)
  }
  
  # Prepare for normalization
  nss <- length(postsumstat_from_posterior[1, ])
  gwt <- !is.na(postsumstat_from_posterior[, 1])
  scaled.postsumstat_posterior <- postsumstat_from_posterior
  scaled.postsumstat_point <- postsumstat_from_point
  
  # Normalize the summary statistics
  for (j in 1:nss) {
    scaled.postsumstat_posterior[, j] <- normalise(postsumstat_from_posterior[, j], postsumstat_from_posterior[, j][gwt])
    target[j] <- normalise(target[j], postsumstat_from_posterior[, j][gwt])
    scaled.postsumstat_point[j] <- normalise(v_postsumstat_from_point[j], postsumstat_from_posterior[, j][gwt])
  }
  
  sum1 <- rowSums((scaled.postsumstat_posterior - target)^2)
  sum2 <- abs(scaled.postsumstat_point - target)
  
  # Adjust distances based on the kernel function
  adjust_distance <- function(dist, kernel) {
    epsilon <- max(dist)
    dist <- dist / epsilon
    if (kernel == "gaussian") {
      dist <- exp(-0.5 * dist^2) / sqrt(2 * pi)
    } else if (kernel == "epanechnikov") {
      dist <- 1 - dist^2
    } else if (kernel == "triangular") {
      dist <- 1 - abs(dist)
    } else if (kernel == "biweight") {
      dist <- (1 - dist^2)^2
    } else if (kernel == "cosine") {
      dist <- cos(pi / 2 * dist)
    }
    dist <- dist / epsilon
    dist[dist == 0] <- NA
    return(log(dist))
  }
  
  
  
  dist1 <- adjust_distance(sqrt(sum1), kernel)
  dist2 <- adjust_distance(sqrt(sum2), kernel)
  
  
  dev <- -2 * mean(dist1, na.rm = TRUE)
  pD <- -2 * dist2
  DIC <- dev + pD

  return(list(expected.deviance = dev, pD = pD, DIC = DIC, dist = dist1))
}
