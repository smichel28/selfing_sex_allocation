
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
                                alpha,
                                color = c('blue', 'red'),
                                xlabel = expression("Generation [x" * 10^3 * "]"),
                                ylabel = 'Value') {
  
  param <- as.factor(param)
  
  # creates plot structure and layout
  par(bty = 'l',
      pch = 16, 
      cex.axis = 1.3,
      cex.lab = 1.5,
      cex.main = 2,
      mar = c(5, 5, 2, 2))
  
  for (a in unique(alpha)) {
    sel_alpha <- alpha == a
    g <- gen[sel_alpha]
    # creates scatter plot
    plot(g/1000, 
         value[sel_alpha], 
         col = color[param[sel_alpha]],
         xlab = xlabel,
         ylab = ylabel,
         main = paste0("\u03B1 = ", a))
    # adds legend
    legend(0.5, max(value[sel_alpha]), legend = levels(param), col = color, pch = 16, cex = 1)
  }
}

mean_through_time <- function(gen,
                              value,
                              param,
                              alpha,
                              delta,
                              migration = NULL,
                              sim,
                              color = c('orange', 'darkblue'),
                              xlabel = 'Generation',
                              ylabel = 'Mean value',
                              xlimit = c(0,100000),
                              ylimit = c(0,1)
                              ) {
  
  al <- unique(alpha)
  del <- unique(delta)
  
  migration_is_null <- is.null(migration)
  
  if (!migration_is_null) mig <- unique(migration)
  else mig <- 1
    
    
  for (m in mig) {
    for (a in al) {
      for (d in del) {
        
        # plotting the evolution of the mean is useless when delta > 0.5 due to branching
        if (d > 0.5 | d == "0.6") next
        
        if (migration_is_null) {
          seeds <- unique(sim[alpha == a & delta == d])
          title <- paste0("\u03B1 = ", a, ", \u03B4 = ", d)
        }
        else {
          seeds <- unique(sim[alpha == a & delta == d & migration == m])
          title <- paste0("\u03B1 = ", a, ", \u03B4 = ", d, ", m = ", m)
        }
        
        for (i in 1:length(seeds)) {
          if (i == 1) {
            plot(gen[sim == seeds[1] & param == 'intercept'], 
                 value[sim == seeds[1] & param == 'intercept'],
                 type = 'l',
                 col = color[1],
                 xlim = xlimit,
                 ylim = ylimit,
                 xlab = xlabel,
                 ylab = ylabel,
                 main = title,
                 bty = 'l')
            lines(gen[sim == seeds[1] & param == 'slope'], 
                  value[sim == seeds[1] & param == 'slope'],
                  col = color[2])
          } else {
            lines(gen[sim == seeds[i] & param == 'intercept'], 
                  value[sim == seeds[i] & param == 'intercept'],
                  col = color[1])
            lines(gen[sim == seeds[i] & param == 'slope'], 
                  value[sim == seeds[i] & param == 'slope'],
                  col = color[2])
          }
        }
          
      }
    }
  }
  
}

## To test ##
# 
# library(dplyr)
# library(tidyr)
# 
# data <- read.table("/home/samuel/Desktop/data/test_6_data_sampled_hapl.tsv", header = TRUE) %>%
#   filter(seed == 7197913992516104192)
# 
# sample_through_time(data$generation, data$Param_value, param = factor(data$param), delta = data$delta)
# 
# data1 <- read.table("/home/samuel/Desktop/data/test_6_data_sampled_hapl.tsv", header = TRUE) %>%
#   filter(seed == 7197913992516104192) %>%
#   pivot_wider(names_from = param, values_from = Param_value)
