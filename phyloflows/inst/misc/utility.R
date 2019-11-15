plot_gp_fit_quantiles <- function(pars, data, pointdata ,title, id,xcol,cred,i) {
  # probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  # cred <- sapply(1:nrow(data$x),
  #                function(n) quantile(pars[,n], probs=probs))
  limits <- c(min(cred[1,id],pointdata[id]),
              max(cred[9,id],pointdata[id]))
  plot(1, type="n", main=paste0(title ,i+15, ' to ', i+16),
       xlab="x", ylab="y", xlim=c(15, 50), ylim=limits)
  polygon(c(data$x[id,xcol], rev(data$x[id,xcol])), c(cred[1,id], rev(cred[9,id])),
          col = c_light, border = NA)
  polygon(c(data$x[id,xcol], rev(data$x[id,xcol])), c(cred[2,id], rev(cred[8,id])),
          col = c_light_highlight, border = NA)
  polygon(c(data$x[id,xcol], rev(data$x[id,xcol])), c(cred[3,id], rev(cred[7,id])),
          col = c_mid, border = NA)
  polygon(c(data$x[id,xcol], rev(data$x[id,xcol])), c(cred[4,id], rev(cred[6,id])),
          col = c_mid_highlight, border = NA)
  lines(data$x[id,xcol], cred[5,id], col=c_dark, lwd=2)
  smoothingSpline = smooth.spline(data$x[id,1], pointdata[id],
                                  spar=0.5, cv=TRUE)
  points(data$x[id,1], pointdata[id], col='black', pch=16, cex=0.8,ylim=c(0,2))	
  lines(smoothingSpline,col='black')
}


# Check transitions that ended with a divergence
check_div <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  
  print(sprintf('%s of %s iterations ended with a divergence (%s%%)',
                n, N, 100 * n / N))
  if (n > 0)
    print('  Try running with larger adapt_delta to remove the divergences')
}

# Check transitions that ended prematurely due to maximum tree depth limit
check_treedepth <- function(fit, max_depth = 10) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  
  print(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                n, N, max_depth, 100 * n / N))
  if (n > 0)
    print('  Run again with max_depth set to a larger value to avoid saturation')
}

# Checks the energy Bayesian fraction of missing information (E-BFMI)
check_energy <- function(fit) {
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  no_warning <- TRUE
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      print(sprintf('Chain %s: E-BFMI = %s', n, numer / denom))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('E-BFMI indicated no pathological behavior')
  else
    print('  E-BFMI below 0.2 indicates you may need to reparameterize your model')
}

# Checks the effective sample size per iteration
check_n_eff <- function(fit) {
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  iter <- dim(extract(fit)[[1]])[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    ratio <- fit_summary[,5][n] / iter
    if (ratio < 0.001) {
      print(sprintf('n_eff / iter for parameter %s is %s!',
                    rownames(fit_summary)[n], ratio))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('n_eff / iter looks reasonable for all parameters')
  else
    print('  n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated')
}

# Checks the potential scale reduction factors
check_rhat <- function(fit) {
  fit_summary <- summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    rhat <- fit_summary[,6][n]
    if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
      print(sprintf('Rhat for parameter %s is %s!',
                    rownames(fit_summary)[n], rhat))
      no_warning <- FALSE
    }
  }
  if (no_warning)
    print('Rhat looks reasonable for all parameters')
  else
    print('  Rhat above 1.1 indicates that the chains very likely have not mixed')
}

check_all_diagnostics <- function(fit) {
  check_n_eff(fit)
  check_rhat(fit)
  check_div(fit)
  check_treedepth(fit)
  check_energy(fit)
}

# Returns parameter arrays separated into divergent and non-divergent transitions
partition_div <- function(fit) {
  nom_params <- extract(fit, permuted=FALSE)
  n_chains <- dim(nom_params)[2]
  params <- as.data.frame(do.call(rbind, lapply(1:n_chains, function(n) nom_params[,n,])))
  
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  params$divergent <- divergent
  
  div_params <- params[params$divergent == 1,]
  nondiv_params <- params[params$divergent == 0,]
  
  return(list(div_params, nondiv_params))
}

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_light_highlight_trans <- c("#C7999980")
c_mid_trans <- c("#B97C7C80")
c_mid_highlight_trans <- c("#A2505080")
c_dark_trans <- c("#8F272780")
c_dark_highlight_trans <- c("#7C000080")

c_light_teal <- c("#6B8E8E")
c_mid_teal <- c("#487575")
c_dark_teal <- c("#1D4F4F")
