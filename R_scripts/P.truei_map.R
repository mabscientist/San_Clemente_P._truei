
# P.truei_map.R -----------------------------------------------------------
#
# MAB 10/29/25
# Code for Figure 1 map for the following manuscript:
# Becker MA, Bell KC, Lopez-Ortiz S, Maldonado JE, Edwards CW, Castañeda-Rico S. (2026) 
# Genomic detection of a novel Peromyscus truei invasion on San Clemente Island, California. 
# Therya 17:277-294. doi: 10.12933/therya.2026.6265
#
# load_packages -----------------------------------------------------------

library(sf)
library(tidyverse)
library(cowplot)
library(ggspatial)
library(ggpubr)


# Channel_Islands_map_for_inset -------------------------------------------

# read in CA county map incl. CI

counties <-
  st_read('../ca_counties/California_County_Boundaries.geojson') %>% 
  set_names(
    names(.) %>% 
      tolower()
  )

# Change generic "Channel Islands" to specific localities

counties$island[59:69] <- c("Sutil", "SBI", "SRI", "SMI", "SCZ", "SNI",
                            "AI-M", "AI-W", "AI-E", "SCL", "CAT")

counties$county_name[59:69] <- c("Sutil", "SBI", "SRI", "SMI", "SCZ", "SNI",
                                 "AI-M", "AI-W", "AI-E", "SCL", "CAT")

my_counties <- c("Los Angeles", "Ventura", "Santa Barbara",
                 "Orange", "San Diego", "San Bernardino", "Riverside",
                 "Sutil", "SBI", "SRI", "SMI", "SCZ", "SNI",
                 "AI-M", "AI-W", "AI-E", "SCL", "CAT")

islands_map <-
  counties %>% 
  filter(county_name %in% my_counties) %>% 
  ggplot() +
  geom_sf() +
  coord_sf(
    xlim = c(-120.7,-117.5),
    ylim = c(32.6,34.5),
    expand = FALSE) +
  #annotate(geom = 'text', x = -120, y = 33.5, label = "California Channel Islands", size = 6, fontface = "bold") +
  annotate(
    geom = 'text', x = -118.4, y = 33.205, 
    label = "Santa Catalina", size = 3.2) +
  
  annotate(
    geom = 'text', x = -118.45, y = 32.72, 
    label = "San Clemente", size = 3.2) +
  
  annotate(
    geom = 'text', x = -119.04, y = 33.38, 
    label = "Santa Barbara", size = 3.2) +
  
  annotate(
    geom = 'text', x = -119.46, y = 33.12, 
    label = "San Nicolas", size = 3.2) +
  
  annotate(
    geom = 'text', x = -120.355, y = 34.17, 
    label = "San Miguel", size = 3.2) +
  
  annotate(
    geom = 'text', x = -120.1, y = 33.8, 
    label = "Santa Rosa", size = 3.2) +
  
  annotate(
    geom = 'text', x = -119.7, y = 34.16, 
    label = "Santa Cruz", size = 3.2) +
  
  annotate(
    geom = 'text', x = -119.3, y = 33.9, 
    label = "Anacapa", size = 3.2) +
  
  annotate(
    geom = 'text', x = -118.244, y = 34.1, 
    label = "Los Angeles", size = 3.2, fontface = "bold") +
  #annotate(
   # geom = 'text', x = -120.4, y = 34.41, 
    #label = "Point Conception", size = 3.5, fontface = "bold") +
  
  theme_void() +
  xlab(NULL) + 
  ylab(NULL) +
  theme(panel.border = element_rect(colour = "black", fill="white", linewidth=1))
 
  #+ annotation_scale()


# csv_processing ----------------------------------------------------------

P_localites <- read.csv("SCL_localities_tibble.csv") %>% 
  st_as_sf(coords = c("x", "y"), crs = 4326)

# SCL map -----------------------------------------------------------------

SCL_map <- counties %>% 
  filter(county_name == "SCL") %>% 
  ggplot() +
  geom_sf() +
  
  geom_sf(data = P_localites, shape = 21, aes(fill = genetic_species, size = n)) +
  scale_size(range = c(1, 8)) +
  scale_fill_manual(labels = c(~italic("P. gambelii"),~italic("P. truei"),"Undetermined"), 
                    values = c("#4477AA","#EE6677","dark gray"),
                    name = "Genetic identification") +
  
  coord_sf(
    xlim = c(-118.7,-118.1),
    ylim = c(32.75,33.15),
    expand = FALSE) +
  
  annotate(
    geom = 'text', x = -118.522	, y = 33.015, 
    label = "Wilson's Cove", size = 3.5) +
  
  annotate(
    geom = 'text', x = -118.488, y = 32.9, 
    label = "Middle Ranch", size = 3.5) +
  
  annotate(
    geom = 'text', x = -118.430875, y = 32.79, 
    label = "China Point", size = 3.5) +
  
  annotate(
    geom = 'text', x = -118.365, y = 32.81, 
    label = "Pyramid Cove", size = 3.5) +
  
  # annotate(
  #   geom = 'text', x = -118.53, y = 32.946, 
  #   label = "1934\nAirfield", size = 3.5, hjust = 0.5) +
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, color = "black"),
        axis.text.y= element_text(color = "black"),
        legend.position = c(0.868,0.098),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  scale_size(guide="none") +
  
  annotation_scale() 


# both_maps ---------------------------------------------------------------

ggdraw() +
  draw_plot(SCL_map) +
  draw_plot(islands_map, x = 0.484, y = 0.514, width = 0.5, height = 0.5)
