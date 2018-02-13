histogram_mode <- function(data, nbins = NA) {
  
  mu <- mean(data, na.rm = TRUE)
  sigma <- sd(data, na.rm = TRUE)
  
  y <- data[data >= mu - 2 * sigma & data <= mu + 2 * sigma]
  
  if (is.na(nbins)) {
    
    h <- hist(y, plot = FALSE)
    
  } else {
    
    h <- hist(y, breaks = nbins, plot = FALSE)
    
  }
  
  # Find peaks
  h_pos <- findpeaks(h$density)
  mode <- h$mids[h_pos$pos_x_max[which.max(h_pos$pos_y_max)]]
  
  return(mode)
  
  
}

noise_estimation_aja <- function(img, nbins = 64) {
  
  img <- (img - min(img))/(max(img) - min(img))
  
  sigma_orig <- histogram_mode(img, nbins)
  sigma <- histogram_mode(img[img < 2 * sigma_orig], nbins)
  
  return(sigma)
  
}

findpeaks <- function(vec, bw = 1, x.coo = c(1:length(vec))) {
  
  pos.x.max <- NULL
  pos.y.max <- NULL
  pos.x.min <- NULL
  pos.y.min <- NULL 	
  
  for (i in 1:(length(vec) - 1)) 	{
    
    if ((i + 1 + bw) > length(vec)) {
      
      sup.stop <- length(vec)
    } else {
      
      sup.stop <- i + 1 + bw
      
    }
    
    if ((i - bw) < 1) { 
      
      inf.stop <- 1
      
    } else {
        
      inf.stop <- i - bw
      
    }
    
    subset.sup <- vec[(i + 1):sup.stop]
    subset.inf <- vec[inf.stop:(i - 1)]
    
    is.max   <- sum(subset.inf > vec[i]) == 0
    is.nomin <- sum(subset.sup > vec[i]) == 0
    
    no.max   <- sum(subset.inf > vec[i]) == length(subset.inf)
    no.nomin <- sum(subset.sup > vec[i]) == length(subset.sup)
    
    if (is.max & is.nomin) {
      
      pos.x.max <- c(pos.x.max, x.coo[i])
      pos.y.max <- c(pos.y.max, vec[i])
      
    }
    
    if (no.max & no.nomin) {
      
      pos.x.min <- c(pos.x.min, x.coo[i])
      pos.y.min <- c(pos.y.min, vec[i])
      
    }
    
  }
  
  return(list(pos_x_max = pos.x.max, 
              pos_y_max = pos.y.max, 
              pos_x_min = pos.x.min, 
              pos_y_min = pos.y.min))
  
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  
  H <- 1.5 * IQR(x, na.rm = na.rm)
  
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  
  return(y[!is.na(y)])
  
}

map_histograms <- function(reference_img, input_img) {
  
  # Compute m1', m2', LIR, HIR for reference image
  foreground <- reference_img[reference_img > 0]
  vec <- remove_outliers(foreground)
  h <- hist(vec, plot = FALSE)
  ref_min <- min(foreground) # m1 = 87
  ref_max <- max(foreground) # m2 = 2929.06958007812
  ref_mean <- mean(foreground) # mu = 841.908577812341
  
  q <- quantile(vec, probs = c(0.1, 0.9), na.rm = TRUE) # q is LIR, HIR (low intensity region, high intensity region)
  LIR <- q[1] # LIR = 410
  HIR <- q[2] # HIR = 1243
  
  # For input image, m1, m2, smin, smax
  foreground <- input_img[input_img > 0]
  vec <- remove_outliers(foreground)
  h <- hist(vec, plot = FALSE)
  input_min <- min(foreground) # m1
  input_max <- max(foreground) # m2
  input_mean <- mean(foreground)
  # img_min_scaled <- (img_min - q[1]) / (q[2] - q[1])
  # img_max_scaled <- (img_max - q[1]) / (q[2] - q[1])
  # foreground_scaled <- (foreground - q[1]) / (q[2] - q[1])
  
  smin <- min(h$mids)
  smax <- max(h$mids)
  
  # Actual transform/mapping
  res <- foreground
  # interval 1: [input_min, smin] to [ref_min, LIR]
  idx1 <- which(foreground >= input_min & foreground <= smin)
  res[idx1] <- (res[idx1] - input_min) / (smin - input_min) * (LIR - ref_min) + ref_min
  
  # interval 2: [smin, input_mean] to [LIR, ref_mean]
  idx2 <- which(foreground >= smin & foreground <= input_mean)
  res[idx2] <- (res[idx2] - smin) / (input_mean - smin) * (ref_mean - LIR) + LIR
  
  # interval 3: [input_mean, smax] to [ref_mean, HIR]
  idx3 <- which(foreground >= input_mean & foreground <= smax)
  res[idx3] <- (res[idx3] - input_mean) / (smax - input_mean) * (HIR - ref_mean) + ref_mean
  
  # interval 4: [smax, input_max] to [HIR, ref_max]
  idx4 <- which(foreground >= smax & foreground <= input_max)
  res[idx4] <- (res[idx4] - smax) / (input_max - smax) * (ref_max - HIR) + HIR
  
  output <- input_img
  output[output > 0] <- res
  
  return(output)
  
}

map_histogram_to_oasis <- function(input_img) {
  
  # Compute m1', m2', LIR, HIR for reference image
  ref_min <- 87
  ref_max <- 2929.06958007812
  ref_mean <- 841.908577812341
  
  # q <- quantile(vec, probs = c(0.1, 0.9), na.rm = TRUE) # q is LIR, HIR (low intensity region, high intensity region)
  LIR <- 410
  HIR <- 1243
  
  # For input image, m1, m2, smin, smax
  foreground <- input_img[input_img > 0]
  vec <- remove_outliers(foreground)
  h <- hist(vec, plot = FALSE)
  input_min <- min(foreground) # m1
  input_max <- max(foreground) # m2
  input_mean <- mean(foreground)

  smin <- min(h$mids)
  smax <- max(h$mids)
  
  # Actual transform/mapping
  res <- foreground
  # interval 1: [input_min, smin] to [ref_min, LIR]
  idx1 <- which(foreground >= input_min & foreground <= smin)
  res[idx1] <- (res[idx1] - input_min) / (smin - input_min) * (LIR - ref_min) + ref_min
  
  # interval 2: [smin, input_mean] to [LIR, ref_mean]
  idx2 <- which(foreground >= smin & foreground <= input_mean)
  res[idx2] <- (res[idx2] - smin) / (input_mean - smin) * (ref_mean - LIR) + LIR
  
  # interval 3: [input_mean, smax] to [ref_mean, HIR]
  idx3 <- which(foreground >= input_mean & foreground <= smax)
  res[idx3] <- (res[idx3] - input_mean) / (smax - input_mean) * (HIR - ref_mean) + ref_mean
  
  # interval 4: [smax, input_max] to [HIR, ref_max]
  idx4 <- which(foreground >= smax & foreground <= input_max)
  res[idx4] <- (res[idx4] - smax) / (input_max - smax) * (ref_max - HIR) + HIR
  
  output <- input_img
  output[output > 0] <- res
  
  return(output)
  
}
