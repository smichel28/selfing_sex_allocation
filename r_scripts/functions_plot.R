
boxplot_pheno <- function(alpha, 
                          value, 
                          param, 
                          delta,
                          color,
                          xlabel = '\u03B1', 
                          ylabel = 'Mean value',
                          jitter.level = 1, 
                          point.size = 1) {
  
  deltas <- sort(unique(delta))
  par(mfrow = c(1,length(deltas)), 
      bty = 'n', 
      pch = 16, 
      cex.axis = 1.7,
      cex.main = 2,
      mar = c(2, 2, 2, 2),
      oma = c(4, 4, 1, 0),
      mgp = c(3.5, 1, 0))
  first <- TRUE
  for (d in deltas) {
    selection <- delta == d
    a <- alpha[selection]
    p <- param[selection]
    y <- value[selection]
    
    lab <- levels(a)
    x <- jitter(as.numeric(a), jitter.level)
    
    if (first) {
      plot(x, y,
           col = color[p],
           xlab = '',
           ylab = '',
           xaxt = 'n',
           yaxt = 'n',
           main = paste('\u03B4 =', d))
      axis(2, pos = 0.5)
      legend(1, 0.8, legend = c('b', 'a'), pch = 16, col = color)
    } else {
      plot(x, y,
           col = color[p],
           xaxt = 'n',
           yaxt = 'n',
           xlab = '',
           ylab = '',
           main = paste('\u03B4 =', d))
    }
    first <- FALSE
    axis(1, at = seq_along(lab), labels = lab)
    
    abline(h = 0.5, lty = 2, col = color[1])
    abline(h = 0, lty = 2, col = color[2])
    
    mtext(xlabel, side = 1, outer = TRUE, line = 2, cex = 1.5)
    mtext(ylabel, side = 2, outer = TRUE, line = 2, cex = 1.5)
    
    for (i in levels(p)) {
      
      select_p <- p==i
      x_box <- as.numeric(a)[select_p]

      boxplot(y[select_p] ~ x_box, add = TRUE,
              boxwex = 0.3, 
              border = color[which(i == levels(p))],
              lwd = 1.5,
              col = NA,
              outline = FALSE,
              axes = FALSE)
    }
    
  }
  
}

