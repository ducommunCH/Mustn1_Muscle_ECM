### general functions ###
library(egg)
library(tidyverse)

default_theme = function() {
  theme_classic() + 
    theme(
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      # axis options
      axis.title = element_text(size = 12, vjust = 0.5),
      axis.text = element_text(size = 11, colour="black")
      
    )
}

# Defines the size of the graphs based on the axis size
ggsave_fixed = function(file, plot = ggplot2::last_plot(), 
                        units = "cm",
                        margin = 1, 
                        plot_width = 10,
                        plot_height = 5, 
                        width = round(dev.size()[1], digits = 1), 
                        height = round(dev.size()[1], digits = 1)) {
  pf = egg::set_panel_size(p = plot,
                           file = NULL, 
                           margin = unit(margin, units),
                           width = unit(plot_width, units), 
                           height = unit(plot_height, units))
  ggplot2::ggsave(file, plot = pf, units = units, width = width, height = height, dpi = 300)
}

