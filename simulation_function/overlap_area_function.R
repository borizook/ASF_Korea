overlap_area <- function(distri_prior, posterior, priors_names){
  
  # Step 1: Determine the range for plotting the prior and posterior
  # xmin and xmax are set to the minimum and maximum of the prior distribution.
  xmin <- min(distri_prior)
  xmax <- max(distri_prior)
  
  # Step 2: Calculate bandwidth for posterior density estimation using Silverman's rule of thumb
  # n1 is the number of posterior samples, and h1 is the bandwidth (smoothing factor) for density estimation.
  n1 <- length(posterior)
  h1 <- 1.06 * min(sd(posterior), IQR(posterior) / 1.34) * n1^(- 1 / 5)
  
  # Step 3: Calculate the kernel density estimate of the posterior
  # The 'density' function estimates the posterior density across the range defined by xmin and xmax.
  d1 <- density(posterior,
                kernel = "gaussian",
                bw = h1,
                from = xmin,
                to = xmax)
  
  # Step 4: Estimate uniform prior density (since prior is uniform, it's a constant)
  p_range = distri_prior[2] - distri_prior[1]  # Step between prior points
  1 / p_range -> p_dens  # Uniform density, since prior is uniformly distributed
  
  # Step 5: Create a data frame to store x (shared between posterior and prior) and corresponding y values for the prior
  # Here, y for the prior is set to p_dens (constant, as it's a uniform distribution).
  d2 = as.data.frame(matrix(ncol = 2, nrow = length(d1$x)))
  colnames(d2) = c("x", "y")
  d2[[1]] = d1$x  # X values
  d2[[2]] = p_dens  # Prior's density values (constant)
  
  # Step 6: Visualizing the overlap by taking the minimum y-value at each x between the prior and posterior distributions
  # y_overlap gives the points where prior and posterior overlap by taking the minimum y value at each x.
  y_overlap <- pmin(d1$y, d2$y)
  
  # Step 7: Create a data frame to hold both the prior, posterior, and overlap curves for plotting
  # 'ddf' will store all the x-values, y-values for each curve, and the type of distribution (posterior, prior, or overlap).
  ddf <- data.frame(x = rep(d1$x, 3),  # X is repeated three times for posterior, prior, and overlap
                    y = c(d1$y, d2$y, y_overlap),  # Combine y-values for posterior, prior, and overlap
                    distri = rep(c("posterior", "distri_prior", "overlap"), 
                                 each = length(d1$x)))
  
  # Step 8: Define a function to calculate the area of the overlap between the prior and posterior
  # 'integrand' takes x values and returns the minimum y values from the posterior and prior, which is then integrated.
  integrand <- function(x) {
    return(pmin(approx(d1$x, d1$y, xout = x)$y,
                approx(d2$x, d2$y, xout = x)$y))
  }
  
  # Step 9: Use numerical integration to calculate the area of overlap between the prior and posterior
  overlap <- integrate(integrand, 
                       lower = xmin,  # Lower bound of integration
                       upper = xmax,  # Upper bound of integration
                       subdivisions = 10000,  # Number of subintervals for integration
                       stop.on.error = FALSE)$value
  
  # Step 10: Plot the prior, posterior, and overlap using ggplot2
  font <- "Helvetica"
  
  gg <- ggplot(ddf, aes(x, y, fill = distri)) + 
    geom_area(position = position_identity(), color = "black", alpha = 0.4) +  # Area plot for posterior, prior, and overlap
    scale_fill_brewer(palette = "YlOrRd") +  # Color palette
    annotate(geom = "text", 
             x = pars_d[,priors_names],  # Position of annotation
             y = 0.2,  # Y-position of the annotation
             label = paste0(round(overlap, 2) * 100, "%"),  # Label the overlap percentage
             size = 5,  # Text size
             color = "black", 
             fontface = "bold") +
    ggplot2::theme(plot.title = ggplot2::element_text(family = font, 
                                                      size = 28, face = "bold", color = "#222222"), 
                   plot.subtitle = ggplot2::element_text(family = font, 
                                                         size = 22, margin = ggplot2::margin(9, 0, 9, 0)), 
                   plot.caption = ggplot2::element_blank(), legend.position = "top", 
                   legend.text.align = 0, legend.background = ggplot2::element_blank(), 
                   legend.title = ggplot2::element_blank(), legend.key = ggplot2::element_blank(), 
                   legend.text = ggplot2::element_text(family = font, size = 18, 
                                                       color = "#222222"), axis.title = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_text(family = font, size = 18, 
                                                     color = "#222222"), axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, 
                                                                                                                                      b = 10)), axis.ticks = ggplot2::element_blank(), 
                   axis.line = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
                   panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"), 
                   panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                   strip.background = ggplot2::element_rect(fill = "white"), 
                   strip.text = ggplot2::element_text(size = 22, hjust = 0))
  
  # Step 11: Return the overlap value and the ggplot object for visualization
  return(list(overlap = overlap, gg = gg))
}