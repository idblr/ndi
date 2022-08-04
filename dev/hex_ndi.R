# ------------------------------------------------------------------------------ #
# Hexsticker for the GitHub Repository idblr/ndi
# ------------------------------------------------------------------------------ #
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: July 23, 2022
#
# Recently modified by: @idblr
# Recently modified on: August 04, 2022
#
# Notes:
# A) Uses the "hexSticker" package
# B) Subplot from an example computation of tract-level NDI (Messer) for Washington, D.C. (2020)
# C) Hexsticker for the GitHub Repository https://github.com/idblr/ndi
# ------------------------------------------------------------------------------ #

############
# PACKAGES #
############

loadedPackages <- c("hexSticker", "ndi")
suppressMessages(invisible(lapply(loadedPackages, library, character.only = TRUE)))

############
# SETTINGS #
############

## Access Key for census data download
### Obtain one at http://api.census.gov/data/key_signup.html
tidycensus::census_api_key("...") # INSERT YOUR OWN KEY FROM U.S. CENSUS API

######################
# SUBPLOT GENERATION #
######################

# NDI 2020
messer2020DC <- ndi::messer(state = "DC", year = 2020, imp = TRUE)

# Tracts 2020
tract2020DC <- tigris::tracts(state = "DC", year = 2020, cb = TRUE)

# Join
DC2020messer <- merge(tract2020DC, messer2020DC$ndi, by = "GEOID")

# Plot of tract-level NDI (Messer) for Washington, D.C. (2020)
dcp <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = DC2020messer, 
                   ggplot2::aes(fill = NDI), 
                   color = NA,
                   show.legend = FALSE) +
  ggplot2::theme_void() + 
  ggplot2::theme(axis.text = ggplot2::element_blank()) +
  ggplot2::scale_fill_viridis_c() +
  ggplot2::labs(fill = "", 
                caption = "")+
  ggplot2::ggtitle("", subtitle = "")

#####################
# CREATE HEXSTICKER #
#####################

s <- hexSticker::sticker(subplot = dcp,
                         package = "ndi",
                         p_size = 75, p_x = 0.55, p_y = 0.75, p_color = "#FDE724", # title
                         s_x = 1.15, s_y = 1.05, s_width = 2.1, s_height = 2.1, # symbol
                         h_fill = "#695488", # inside
                         h_color = "#440C54", # outline
                         dpi = 1000, # resolution
                         filename = "man/figures/ndi.png",
                         white_around_sticker = F)

# -------------------------------- END OF CODE --------------------------------- #
