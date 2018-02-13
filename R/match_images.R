match_images <- function(source, target, nbins = 64) {
  
  source <- source / max(source) * 255
  target <- target / max(target) * 255
  
  h_source <- hist(x = as.vector(source), plot = FALSE, breaks = nbins)
  h_target <- hist(x = as.vector(target), plot = FALSE, breaks = nbins)
  
  counts_source <- h_source$counts / sum(h_source$counts)
  counts_target <- h_target$counts / sum(h_target$counts)
  
  counts_source <- cumsum(counts_source)
  counts_target <- cumsum(counts_target)
  
  require(dtw)
  
  t <- dtw::dtw(x = counts_source, y = counts_target, step = symmetric1)
  
  tv <- 1
  sv <- 1
  
  input_interval <- min(h_source$breaks)
  reference_interval <- min(h_target$breaks)
  
  
  while (tv <= max(t$index2) && sv <= max(t$index1)) {
    
    which_t <- which(t$index2 == tv)
    which_s <- which(t$index1 == sv)
    
    if (length(which_s) > 1) {
      
      t_values <- t$index2[which_s]
      input_interval <- c(input_interval, h_source$breaks[sv + 1])
      reference_interval <- c(reference_interval, h_target$breaks[max(t_values) + 1])
      
      sv <- sv + 1
      tv <- max(t_values) + 1
      
      next
    } 
    
    if (length(which_t) > 1) {
      
      s_values <- t$index1[which_t]
      input_interval <- c(input_interval, h_source$breaks[max(s_values) + 1])
      reference_interval <- c(reference_interval, h_target$breaks[tv + 1])
      
      tv <- tv + 1
      sv <- max(s_values) + 1
      
      next
      
    }
    
    if ((length(which_s) == 1) && (length(which_t) == 1)) {
      
      input_interval <- c(input_interval, h_source$breaks[sv + 1])
      reference_interval <- c(reference_interval, h_target$breaks[tv + 1])
      
      sv <- sv + 1
      tv <- tv + 1
      
    }
    
  }
  
  new_s <- map_intervals(input = source, 
                         input_intervals = input_interval, 
                         reference_intervals = reference_interval)
  
  return(new_s)
  
}