plot_age_source_recipient <- function(data, title, lab, outdir){
  
  data <- data[!is.na(age_infection.SOURCE) & !is.na(age_infection.RECIPIENT)]
  data[, `Cohort round recipient` := cohort_round.RECIPIENT]
  data[, `Cohort round source` := cohort_round.SOURCE]
  data[, `Community recipient` := comm.RECIPIENT]
  data[, `Community source` := comm.SOURCE]
  
  # all pairs
  p <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point() + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) 
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_AllPairs_', lab, '.png')), w = 4, h = 4)
  
  # by cohort round
  p <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point() + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    facet_grid(.~`Cohort round source`, label = 'label_both') 
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_CohortSource_', lab, '.png')), w = 7, h = 7)
  
  p <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point() + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    facet_grid(.~`Cohort round recipient`, label = 'label_both') 
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_CohortRecipient_', lab, '.png')), w = 7, h = 7)
  
  # by community
  p <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point() + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    facet_grid(.~`Community source`, label = 'label_both') 
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_CommunitySource_', lab, '.png')), w = 7, h = 7)
  
  p <- ggplot(data, aes(x = age_infection.SOURCE, y = age_infection.RECIPIENT)) + 
    geom_point() + 
    labs(x = 'Age at infection source', y = 'Age at infection recipient') +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed', col = 'grey50') + 
    theme_bw() + 
    coord_fixed() +
    scale_x_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT)))+
    scale_y_continuous(limits = range(c(data$age_infection.SOURCE, data$age_infection.RECIPIENT))) +
    ggtitle(paste0(title, ' - ', paste0(nrow(data), ' pairs'))) + 
    facet_grid(.~`Community recipient`, label = 'label_both') 
  ggsave(p, filename = file.path(outdir, paste0('AgeInfection_CommunityRecipient_', lab, '.png')), w = 7, h = 7)
  
}

plot_hist_age_infection <- function(pairs, outdir){
  
  pairs[, Sex := sex.SOURCE]
  p1 <- ggplot(pairs, aes(x = age_infection.SOURCE)) + 
    geom_histogram(bins = 30) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection source') +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  pairs[, Sex := sex.SOURCE]
  p2 <- ggplot(pairs, aes(x = age_infection.RECIPIENT)) + 
    geom_histogram(bins = 30) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient')  +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  p <- ggarrange(p1, p2, ncol = 2)
  
  file = file.path(outdir, paste0('hist_age_infection_source_', lab, '.png'))
  cat('saving', file)
  ggsave(p, file = file, w = 6, h = 6)
  
  return(p)
}

plot_hist_age_infection_diff_threshold <- function(pairs, outdir){
  
  chain <- keep.likely.transmission.pairs(as.data.table(dchain), 0.5)
  pairs <- pairs.get.meta.data(chain, meta_data)
  pairs$threshold = '0.5'
  
  chain <- keep.likely.transmission.pairs(as.data.table(dchain), 0.6)
  pairs2 <- pairs.get.meta.data(chain, meta_data)
  pairs2$threshold = '0.6'
  
  pairs = rbind(pairs, pairs2)
  
  if(!include.mrc){
    cat('Keep only pairs in RCCS')
    pairs <- pairs[cohort.RECIPIENT == 'RCCS' & cohort.SOURCE == 'RCCS']
  }
  if(include.only.heterosexual.pairs){
    cat('Keep only heterosexual pairs')
    pairs <- pairs[(sex.RECIPIENT == 'M' & sex.SOURCE == 'F') | (sex.RECIPIENT == 'F' & sex.SOURCE == 'M')]
  }
  
  pairs[, Sex := sex.SOURCE]
  p1 <- ggplot(pairs, aes(x = age_infection.SOURCE)) + 
    geom_density(aes( group = threshold, fill = threshold), alpha = 0.5) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection source') +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  pairs[, Sex := sex.SOURCE]
  p2 <- ggplot(pairs, aes(x = age_infection.RECIPIENT)) + 
    geom_density(aes( group = threshold, fill = threshold), alpha = 0.5) + 
    facet_wrap(~Sex, ncol = 1, label = 'label_both') + 
    theme_bw() + 
    labs(x = 'Age at infection recipient')  +
    scale_x_continuous(limits = range(c(pairs$age_infection.SOURCE, pairs$age_infection.RECIPIENT)))
  
  p <- ggarrange(p1, p2, ncol = 2, common.legend = T, legend = 'bottom')
  
  file = file.path(outdir, paste0('hist_age_infection_source_', gsub('(.+)_threshold.*', '\\1', lab), '.png'))
  cat('saving', file)
  ggsave(p, file = file, w = 6, h = 6)
  
  return(p)
}
