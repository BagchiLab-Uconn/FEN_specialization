library(tidyverse)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(png) # Needed for compass rose
library(grid) # Needed for compass rose
library(tidyterra)

# State of connecticut
ct.shp <- st_read("../Fragment-Spatial/data/shapefiles/ct2010_state_wgs84.shp")

# State of connecticut, with counties
ct.counties.shp <- st_read("../Fragment-Spatial/data/shapefiles/ct2010_counties_wgs84.shp")

# FEN site polygons and metadata
fen.shp <- st_read("../Fragment-Spatial/data/shapefiles/ct2015_FEN_core_wgs84.shp")

# Core forest for the whole state of connecticut
ct.coreforest.shp <- st_read("../Fragment-Spatial/data/shapefiles/ct2015_core_wgs84.shp")

# Import compass rose image and convert to a raster-based grob for ggplot
compass <- readPNG("../Fragment-Spatial/data/compass_rose.png") %>% rasterGrob()

sf_use_s2(FALSE)
# Edit the map of core forest in connecticut
ct.coreforest.shp <- ct.coreforest.shp %>% 
    
    # Crop to the pixels inside CT
    st_intersection(ct.counties.shp)

#### NaturalEarth Data
# Get sf data for the USA map using RNaturalEarth
USA <- ne_countries(scale = 'medium', returnclass = 'sf') %>% 
    filter(admin == "United States of America")

# Get sf data for CT, matching the USA map
CT <- ne_states(country = 'United States of America', returnclass = "sf") %>% 
    filter(name == "Connecticut")

usa.map <- ggplot() + 
    
    # Add the united states
    geom_sf(data = USA, fill = NA) + 
    
    # Add the state of connecticut
    geom_sf(data = CT, fill = "blue") +
    
    # Add a point for connecticut
    annotate("point", x = 2220000, y = -72.593930, color = "red", size = 3) + 
    
    # Add a rectangle around the map
    annotate("rect", xmin = -2080000, xmax = 2550000, 
             ymin = -2200000, ymax = 800000, fill = NA, color = "black") + 
    
    # Crop the map 
    coord_sf(crs = st_crs(2163), 
             xlim = c(-1880000, 2350000), ylim = c(-2060000, 670000)) + 
    
    # Add custom map theme
    theme(
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_line(colour = "transparent"), 
        panel.background = element_rect(fill = "transparent", colour = NA)
    )

# Base map of connecticut
ct.map <- ggplot() + 
    
    # Add the outline of CT, with counties
    geom_sf(data = ct.counties.shp, color = "black", fill = NA) +
    
    # Add the core forest background across the state of CT
       geom_sf(data = ct.coreforest.shp, color = "#0c851a", fill = "#0c851a50",
               size = 0.05) +
    # # Add a scale bar
       annotation_scale(location = "bl", pad_x = unit(0.4, "native"),
                        width_hint = 0.15)+
    # # # Add compass rose
    # annotation_custom(compass, xmin = -73, xmax = -72.45, ymin = 40.94, 
    #                    ymax = 41.19) +
    # Add USA inset
    annotation_custom(ggplotGrob(usa.map), xmin = -72.39, xmax = -71.73, 
                      ymin = 40.71, ymax = 41.38) + 
    ggthemes::theme_tufte() # # Add the map theme
 

fen.shp$BlockID <- factor(fen.shp$BlockID, levels = 
                              sample(unique(fen.shp$BlockID)))

# Construct the custom map for Riley's Frontiers paper
plot_map <- ct.map + 
    # # Add the outline of CT
  #   geom_sf(data = ct.shp, color = "black", fill = NA) + 
   
    # Add points for all of the FEN sites
    geom_point(aes(x = Longitude, y = Latitude, color = BlockID, 
                   shape = SizeClass), data = fen.shp, size = 3, alpha = 0.9) +
    guides(colour = "none") +
    scale_shape_manual(values = c(15, 16, 17), name = "Fragment size class") +
    theme(legend.position = "top") + 
    guides(shape = guide_legend(title.position = "top", title.hjust = 0.5))




# Find the vertices for an equilateral triangle with centroid [0,0]
# Top vertex is 2/3 of the length of the line bisecting it with the bottom side
# Side specifies the length of a side in meters
vertices <- function(side = 25, start = 1, type = "point", offset = 2, 
                     center = c(0,0)) {
    
    # Make a data frame to hold data on the three vertices of a triangle
    points <- data.frame(
        
        # Add the point type
        type = type,
        
        # Number the points
        point = seq(start, start + 2), 
        
        # Calculate the 3 x-coordinates, starting from lower left
        x = c(side/2*-1 + center[1], center[1], side/2 + center[1]), 
        
        # Calculate the 3 y-coordinates, starting from lower left
        y = c((sqrt((side^2) - ((side/2)^2))/-3) + center[2], 
              (sqrt((side^2) - ((side/2)^2))/3*2) + center[2], 
              (sqrt((side^2) - ((side/2)^2))/-3) + center[2]))
    
    # Calculate offsets for bi-directional arrows along each triangle side: 
    # (left-top, right-top, bottom).
    # Calculated w/trig, using the offset (dist from vertex) as hypotenuse
    x1 <- points$x + c(offset * cos(pi / 3), offset * sin(pi / 6), -offset)
    y1 <- points$y + c(offset * sin(pi / 3), -offset * cos(pi / 6), 0)
    x2 <- points$x[c(2,3,1)] + c(-offset * sin(pi / 6), -offset * cos(pi / 3), offset)
    y2 <- points$y[c(2,3,1)] + c(-offset * cos(pi / 6), offset * sin(pi / 3), 0)
    
    # Combine the vertices with the arrow offsets
    points <- cbind(points, x1, y1, x2, y2)
    
    # Return the 
    return(points)
}


##### Triangle Vertices & Deer Points #####

# Calculate 25m vertices, numbering starting with 2, arrow offset of 3m
beat.pts <- vertices(25, 2, "beat", offset = 3)

# Calculate 125m vertices, numbering starting with 5, arrow offset of 10m
burlap.pts <- vertices(125, 5, "burlap", offset = 10)

# Calculate all the vegetation plot vertices
veg.pts <-  vertices(25, 2, "veg", offset = 4) 


# Calculate 50 m grid of deer plot points
deer.pts <- data.frame(type = "deer", point = paste0("DP", seq(1,25)),
                       x = rep(seq(-100, 100, by = 50), each = 5), 
                       y = rep(seq(100, -100, by = -50), times = 5),
                       x1 = NA, y1 = NA, x2 = NA, y2 = NA)

# Combine all of the points
fen.pts <- rbind(beat.pts,  veg.pts, deer.pts)


# Make a dataframe to hold vertices of the corners of caterpillar beat plots
# First, we calculate one corner of CP1
beat.plots <- beat.pts %>% 
    
    # Only keep the point #, and coordinates columns
    select(point, x, y) %>% 
    
    # Remove Vegplot point #3, because it has a different arrangement
    filter(point != 3) %>% 
    
    # Calculate the position of 2_CP1, SE corner, and 4_CP1, SE corner
    mutate(x = x + c(-2.5, 7.5), plot = paste0(point, "_CP1"), corner = "SE")

# Next, we add the same corner (SE), to CP2, CP3, and CP4
beat.plots <- beat.plots %>% 
    
    # Calculate the SE corner for each caterpillar beat plot
    rbind(beat.plots %>% mutate(x = x + 7.5, y = y - 7.5, 
                                plot = paste0(point, "_CP2")), 
          beat.plots %>% mutate(y = y - 15, 
                                plot = paste0(point, "_CP3")), 
          beat.plots %>% mutate(x = x - 7.5, y = y - 7.5, 
                                plot = paste0(point, "_CP4"))) 


# Finally, we add the three other corners to each caterpillar beat plot.
beat.plots <- beat.plots %>% 
    
    # Calculate the other three corners, using the SE corner as a reference
    rbind(beat.plots %>% mutate(y = y + 5, corner = "NE"), 
          beat.plots %>% mutate(x = x - 5, y = y + 5, corner = "NW"), 
          beat.plots %>% mutate(x = x - 5, corner = "SW")
    ) %>%
    
    # Arrange the data by Vegplot and Caterpillar beat plot
    arrange(plot)

beat.plots

# Vegplot 3 oriented differently; its plots need different treatmenty
plot3 <- beat.pts %>%     
    
    # Only keep the point #, and coordinates columns
    select(point, x, y) %>% 
    
    # Use only vegetation plot 3
    filter(point == 3) %>% 
    
    # Calculate the E corner of 3_CP3
    mutate(x = x - sqrt((2.5^2)/2), y = y + sqrt((2.5^2)/2), 
           plot = paste0(point, "_CP3"), corner = "E")

# Next, we add the same corner (E), to CP1, CP2, and CP4
plot3 <- plot3 %>% 
    
    # Calculate the SE corner for each caterpillar beat plot
    rbind(plot3 %>% mutate(x = x + sqrt((2.5^2)/2)*2 + sqrt(5^2*2), 
                           plot = paste0(point, "_CP2")), 
          plot3 %>% mutate(x = x + sqrt((2.5^2)/2)*2 + sqrt(5^2*2), 
                           y = y + sqrt((2.5^2)/2)*2 + sqrt(5^2*2), 
                           plot = paste0(point, "_CP1")), 
          plot3 %>% mutate(y = y + sqrt((2.5^2)/2)*2 + sqrt(5^2*2), 
                           plot = paste0(point, "_CP4"))
    )


# Finally, we add the three other corners to each caterpillar beat plot.
plot3 <- plot3 %>%
    
    # Calculate the other three corners, using the E corner as a reference
    rbind(plot3 %>% mutate(x = x - sqrt(5^2*2)/2, y = y + sqrt(5^2*2)/2, 
                           corner = "N"), 
          plot3 %>% mutate(x = x - sqrt(5^2*2), 
                           corner = "W"), 
          plot3 %>% mutate(x = x - sqrt(5^2*2)/2, y = y - sqrt(5^2*2)/2, 
                           corner = "S")) %>%
    
    # Arrange the data by Caterpillar beat plot
    arrange(plot)


# Add the caterpillar plots from vegplot 3 to the rest of the caterpillar plots
beat.plots <- beat.plots %>% rbind(plot3) %>% arrange(plot)

##### Vegetation Plots #####


# Make a dataframe to hold vertices of the corners of vegetation plots
# First, we calculate one corner of of the western vegplot
veg.plots <- veg.pts %>% 
    
    # Filter out the beat plots, they need separate treatment
    filter(point > 4) %>%
    
    # Only keep the point #, and coordinates columns
    select(point, x, y) %>% 
    
    # Calculate the NE corner of the western vegplot at each vertex
    mutate(x = x - 12.5, plot = paste0(point, "_1"), corner = "NE")

# Add the eastern vegplots and vegplots 2 & 4
veg.plots <- veg.plots %>% 
    
    # Add NE corner of additional vegetation plots
    rbind(
        
        # Add the second (eastern) vegplot at each vertex
        veg.plots %>% mutate(x = x + 35, plot = paste0(point, "_2")),
        
        # Add Vegplot #2, NE corner
        veg.pts %>% filter(point == 2) %>% select(point, x, y) %>% 
            mutate(plot = 2, corner = "NE"),
        
        # Add Vegplot #4, NE corner
        veg.pts %>% filter(point == 4) %>% select(point, x, y) %>% 
            mutate(x = x + 10, plot = 4, corner = "NE")
        
    )

# Finally, we add the three other corners to each vegetation plot.
veg.plots <- veg.plots %>%
    
    # Calculate the other three corners, using the NE corner as a reference
    rbind(veg.plots %>% mutate(y = y - 10, corner = "SE"), 
          veg.plots %>% mutate(x = x - 10, y = y - 10, corner = "SW"), 
          veg.plots %>% mutate(x = x - 10, corner = "NW")) %>%
    
    # Arrange the data by vegetation plot
    arrange(plot)


# Vegplot 3 is oriented differently; needs special treatmenty
vegplot3 <- veg.pts %>% 
    
    # Only keep the point #, and coordinates columns
    select(point, x, y) %>% 
    
    # Use only vegetation plot 3
    filter(point == 3) %>% 
    
    # Make the vertex the S corner
    mutate(plot = 3, corner = "S")


# Finally, we add the three other corners vegplot 3
vegplot3 <- vegplot3 %>% 
    
    # Calculate the other three corners, using the S corner as a reference
    rbind(
        vegplot3 %>% mutate(x = x + sqrt((10^2)*2)/2, y = y + sqrt((10^2)*2)/2, 
                            corner = "E"), 
        vegplot3 %>% mutate(y = y + sqrt((10^2)*2), corner = "N"), 
        vegplot3 %>% mutate(x = x - sqrt((10^2)*2)/2, y = y + sqrt((10^2)*2)/2, 
                            corner = "W")
    )

# Combine all the vegetation plots
veg.plots <- veg.plots %>% rbind(vegplot3)

arrow_locs <-  tibble(
    x = c(-12, 15, -22), y = c(-7.5, -1, -12.5), 
    xend = c(0, 20, -13), yend = c(14, -1, -12.5), 
    dist = c("25 m", "5 m", "10 m"),
    xlab = c(-9, 17.5, -17.5), ylab = c(2.5, 2, -9.5), 
    ang = c(60, 0, 0))
# Create a plot combining burlap trees and caterpillar beat plots
plot_layout <- ggplot() + 
    # # Add the deer survey plots
    # geom_point(aes(x = x, y = y, color = "Deer Plots"), data = deer.pts,
    #            size = 1) +
    # Add the caterpillar beat plots
    geom_polygon(aes(x = x, y = y, group = plot, fill = "Beat Plots"), 
                 data = beat.plots) +
    
    # Add the vegetation
    geom_polygon(aes(x = x, y = y, group = plot, fill = "Vegetation Plots"), 
                 data = veg.plots) + 
   
    
    # Add legend and colors
    scale_color_manual(name = "", values = c("blue", "forestgreen"), 
                       aesthetics = c("color", "color", "color", "fill", "fill")) + 
    
    # Add arrows showing distance between triangle vertices
    # Only showing a subset of the arrows using filter()
    geom_segment(data = arrow_locs, 
                 aes(x = x, y = y, xend = xend, yend = yend), 
                 size = 0.6,
                 arrow = arrow(ends = "both", type = "closed", 
                               length = unit(0.15, "cm"))) +
     
    geom_text(data = arrow_locs, aes(label = dist, x=xlab, y = ylab, angle= ang)) +
    # Add compass rose
    annotation_custom(compass, xmin = -30, xmax = -10, ymin = 10, 
                      ymax = 30) +
    coord_equal() +    ggthemes::theme_tufte() +
    theme(
        # Turn off axes
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank()
    )

library(patchwork)
sitemap <- (plot_map + plot_layout) +
    plot_annotation(tag_levels = "A") +
    plot_layout(widths = c(1.5, 1))


ggsave(sitemap, file = "figures/sitemap.pdf", 
       height = 4, width = 7)