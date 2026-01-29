
boxplot_pheno <- function(alpha, 
                          value, 
                          param, 
                          delta,
                          color,
                          xlabel = '\u03B1', 
                          ylabel = 'Mean value',
                          jitter.level = 1, 
                          point.size = 1) {
  
  # creates plot structure and layout
  deltas <- sort(unique(delta))
  par(mfrow = c(1,length(deltas)), 
      bty = 'n', 
      pch = 16, 
      cex.axis = 1.3,
      cex.main = 2,
      mar = c(2, 2, 2, 2),
      oma = c(4, 4, 1, 0),
      mgp = c(3.5, 1, 0))
  
  # creates plots
  first <- TRUE
  for (d in deltas) {
    selection <- delta == d
    a <- alpha[selection]
    p <- param[selection]
    y <- value[selection]
    
    lab <- levels(a)
    x <- jitter(as.numeric(a), jitter.level)
    
    # creates points
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
    
    # adds ablines to mark starting value
    abline(h = 0.5, lty = 2, col = color[1])
    abline(h = 0, lty = 2, col = color[2])
    
    # adds boxplots for each parameter (slope & intercept)
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
  
  # adds axis labels
  mtext(xlabel, side = 1, outer = TRUE, line = 2, cex = 1.5)
  mtext(ylabel, side = 2, outer = TRUE, line = 2, cex = 1.5)
  
}

sample_through_time <- function(gen,
                                value,
                                param,
                                delta,
                                color = c('blue', 'red'),
                                xlabel = expression("Generation [x" * 10^3 * "]"),
                                ylabel = 'Value') {
  
  # creates plot structure and layout
  par(mfrow=c(1,1), 
      bty = 'l',
      pch = 16, 
      cex.axis = 1.3,
      cex.lab = 1.5,
      cex.main = 2,
      mar = c(5, 5, 2, 2))
  
  # creates scatter plot
  plot(gen/1000, value, 
       col = color[param],
       xlab = xlabel,
       ylab = ylabel)
  # adds legend
  legend(0.5, 1.15, legend = levels(param), col = color, pch = 16, cex = 1.5)
}

## To test ##
library(dplyr)
library(tidyr)

data <- read.table("/home/samuel/Desktop/data/test_6_data_sampled_hapl.tsv", header = TRUE) %>%
  filter(seed == 7197913992516104192)

sample_through_time(data$generation, data$Param_value, param = factor(data$param), delta = data$delta)

data1 <- read.table("/home/samuel/Desktop/data/test_6_data_sampled_hapl.tsv", header = TRUE) %>%
  filter(seed == 7197913992516104192) %>%
  pivot_wider(names_from = param, values_from = Param_value)

