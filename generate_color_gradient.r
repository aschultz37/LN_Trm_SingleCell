library(RColorBrewer)

# change the colors at ends of gradient
start_color <- '#1B9E77'
end_color   <- '#D95F02'

# set the number of colors in the palette/gradient
num_colors  <- 8

# generate the palette and display it
paletteFunc <- colorRampPalette(c(start_color, end_color))
palette     <- paletteFunc(num_colors)
barplot(1:num_colors, col=palette)