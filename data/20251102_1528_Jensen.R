########################################################################
# Overview:
# Simulation model for assessing how metabolic inequality affects 
# ecosystem processes. We assume that microbial respiration is influenced
# by cellular enzyme levels according to a non-linear saturating function
# (i.e., Michaelis-Menten). With this model, we assess how, not only variance, 
# but also skew influences aggregate function, by comparing a Gaussian 
# distribution of enzyme levels to lognormal distribution with different 
# levels of sigma (standard deviation in log space)
########################################################################


########################################################################
# Part 1. Set up the conditions for simulating enzymes and respiration

# Michaelis-Menten dynamics with different half-saturation constants (Km)
########################################################################

set.seed(42)

# ----------------------
# Parameters & functions
# ----------------------

n <- 10^6 # number of cells
mean_norm <- 1e5 # mean of distribution
sd_norm   <- 0.25 * mean_norm # sd of for Gaussian distribution
Vmax <- 1 # standardized maximum respiration
sigma_vals <- seq(0.01, 3.5, by = 0.1) # log-space sd for skew in lognormal sims
Km_vals <- c(1e2, 1e3, 1e4, 1e5) # Km values
colors <- c("red","darkorange","darkgreen","blue")[1:length(Km_vals)]
f_mm <- function(x, Km, Vmax = 1) (Vmax * x) / (Km + x) # Michaelis-Menten function

# --------------------------------------------------------------
# Part 2: Simulations of respiration with Gaussian and lognormal
# --------------------------------------------------------------

# Pull cells from a Gaussian distribution
activity_norm <- rnorm(n, mean = mean_norm, sd = sd_norm)

# Replace any negative values with very small positive ones
activity_norm_pos <- pmax(activity_norm, 1e-12)

# Apply Michaelis-Menten function to Gaussian-sampled cells
# Calculates expected (mean) community respiration 
baseline_resp <- sapply(Km_vals, function(Km) mean(f_mm(activity_norm_pos, Km = Km, Vmax = Vmax)))

# Create matrices for output of Km and sigma combinations
resp_ratio_matrix <- matrix(NA, nrow = length(sigma_vals), ncol = length(Km_vals)) # R ratios
Ef_lognorm_matrix <- matrix(NA, nrow = length(sigma_vals), ncol = length(Km_vals)) # mean R for lognormal
for (k in seq_along(Km_vals)) {
  Km <- Km_vals[k]
  for (i in seq_along(sigma_vals)) {
    sigma <- sigma_vals[i]  # divided by the Gaussian baseline for Km = Km_vals[k]
    mu <- log(mean_norm) - 0.5 * sigma^2
    # Sets lognormal’s meanlog so arithmetic (unlogged) mean of the lognormal is approximately mean_norm. 
    # For a lognormal with parameters (meanlog = μ, sdlog = σ), the arithmetic mean is:exp(μ + σ^2/2). 
    # Solving for μ gives μ = log(mean) − σ^2/2. 
    # Lets you vary skew (σ) while holding arithmetic mean fixed.
    # Preserve arithmetic mean across σ values
    activity_lognorm <- rlnorm(n, meanlog = mu, sdlog = sigma)
    # Draw n samples from the lognormal with parameters mu and sigma
    # Models per-cell activity under multiplicative variability (skewed), controlling skew via sigma
    activity_lognorm <- pmax(activity_lognorm, 1e-12)
    # Mirrors the Gaussian clipping behavior for numerical consistency
    mean_resp_lognorm <- mean(f_mm(activity_lognorm, Km = Km, Vmax = Vmax))
    # Calculate average respiration across n sampled cells under the lognormal activity for current Km
    resp_ratio_matrix[i, k] <- mean_resp_lognorm / baseline_resp[k]
    # Store ratio of lognormal mean respiration to the Gaussian-baseline mean respiration for this (sigma, Km)
    # Values >1 indicate the lognormal yields higher predicted community respiration than the Gaussian baseline (with same arithmetic mean)
    # Values <1 indicate the lognormal yields lower predicted community respiration than the Gaussian baseline (with same arithmetic mean)
    Ef_lognorm_matrix[i, k] <- mean_resp_lognorm
    # Store raw mean respiration (useful for plotting in absolute units or checking numeric stability).
  }
}

# --------------------
# Part 3: Make figure

# Respiration ratio as function of skew
# Different colored lines for different Km values
# --------------------

cutoff <- 3.2
idx <- sigma_vals <= cutoff

op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)

par(mar = c(7, 9, 4, 2) + 0.1, mgp = c(4, 1, 0), lwd = 1.6)
par(pty = "s")

# Main plot
plot(sigma_vals[idx], resp_ratio_matrix[idx, 1],
     type = "l", lwd = 2.2, col = colors[1],
     ylim = range(resp_ratio_matrix, na.rm = TRUE),
     xlab = "", ylab = "",
     xaxt = "n", yaxt = "n",
     xlim = c(0, 3.7),
     cex.lab = 1.8)

# Bottom axis
axis(1, at = c(0,1,2,3), labels = c("0.0","1.0","2.0","3.0"),
     cex.axis = 1.6, lwd = 1.8, lwd.ticks = 1.8)

# Left axis
y_at <- pretty(range(resp_ratio_matrix, na.rm = TRUE), n = 5)
axis(2, at = y_at, labels = format(y_at, trim = TRUE, digits = 3),
     cex.axis = 1.6, lwd = 1.8, lwd.ticks = 1.8, las = 1)

# Add ticks: top (x) and right (y) axes
axis(3, at = c(0, 1, 2, 3), labels = FALSE, tck = -0.02,
     lwd = 1.8, lwd.ticks = 1.8)
axis(4, at = y_at, labels = FALSE, tck = -0.02,
     lwd = 1.8, lwd.ticks = 1.8)

# Add axis labels
mtext(side = 1, text = expression("Standard deviation (" * sigma * ")"),
      line = 4.5, cex = 1.95)
mtext(side = 2, text = "Respiration ratio", line = 4.5, cex = 1.95)

# Add colored text next to lines for Km values
for (k in seq_along(Km_vals)) {
  lines(sigma_vals[idx], resp_ratio_matrix[idx, k], col = colors[k], lwd = 2)
  y_end <- tail(resp_ratio_matrix[idx, k], 1)
  x_end <- tail(sigma_vals[idx], 1)
  exp_val <- as.integer(round(log10(Km_vals[k])))
  
  # Improved label format: K[m] = 10^x
  label_expr <- as.expression(bquote(K[m] == 10^.(exp_val)))
  text(x_end + 0.02, y_end, labels = label_expr,
       col = colors[k], cex = 1.25, pos = 4)
}

abline(h = 1, col = "gray50", lty = 2)
box(lwd = 1.8)

