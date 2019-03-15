#!/usr/bin/env RScript
# Title: ggplot_n_samples.R
# Authors: Alexandre Pellan Cheng
# Brief description: Convenience function to plot ns at the bottom of a boxplot-like figure.

n_samples<-function(dat, x, y, y_nudge, fontsize){
  labs=c(paste("n=", table(x), sep=""))
  xs=levels(x)
  ys=y
  df_x_y<-data.frame(x,y)
  heights=aggregate(.~x, df_x_y, max)$y
  return(geom_text(data=data.frame(), aes(x=xs, y=heights, label=labs), nudge_y = y_nudge, family="Helvetica", size=fontsize/ggplot2:::.pt))
}