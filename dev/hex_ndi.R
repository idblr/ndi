# ----------------------------------------------------------------------------------------------- #
# Hexagon sticker for the GitHub Repository idblr/ndi
# ----------------------------------------------------------------------------------------------- #
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: 2022-07-23
#
# Recently modified by: @idblr
# Recently modified on: 2024-07-06
#
# Notes:
# A) Uses the 'hexSticker' package
# B) Subplot from an example computation of tract-level NDI (Messer) for Washington, D.C. (2020)
# C) Hexsticker for the GitHub Repository https://github.com/idblr/ndi
# ----------------------------------------------------------------------------------------------- #

# -------- #
# PACKAGES #
# -------- #

loadedPackages <- c('ggplot2', 'hexSticker', 'ndi', 'tidycensus', 'tigris')
suppressMessages(invisible(lapply(loadedPackages, library, character.only = TRUE)))

# -------- #
# SETTINGS #
# -------- #

## Access Key for census data download
### Obtain one at http://api.census.gov/data/key_signup.html
census_api_key('...') # INSERT YOUR OWN KEY FROM U.S. CENSUS API

# ------------------ #
# SUBPLOT GENERATION #
# ------------------ #

# NDI 2020
messer2020DC <- messer(state = 'DC', year = 2020, imp = TRUE)

# Tracts 2020
tract2020DC <- tracts(state = 'DC', year = 2020, cb = TRUE)

# Join
DC2020messer <- merge(tract2020DC, messer2020DC$ndi, by = 'GEOID')

# Plot of tract-level NDI (Messer) for Washington, D.C. (2020)
dcp <- ggplot() + 
  geom_sf(data = DC2020messer, aes(fill = NDI), color = NA, show.legend = FALSE) +
  theme_void() + 
  theme(axis.text = element_blank()) +
  scale_fill_viridis_c() +
  labs(fill = '', caption = '')+
  ggtitle('', subtitle = '')

# ---------------------- #
# CREATE HEXAGON STICKER #
# ---------------------- #

s <-sticker(
  subplot = dcp,
  package = 'ndi',
  p_size = 75, p_x = 0.55, p_y = 0.75, p_color = '#FDE724', # title
  s_x = 1.15, s_y = 1.05, s_width = 2.1, s_height = 2.1, # symbol
  h_fill = '#695488', # inside
  h_color = '#440C54', # outline
  dpi = 1000, # resolution
  filename = file.path('man', 'figures', 'ndi.png'),
  white_around_sticker = FALSE
)

# ----------------------------------------- END OF CODE ----------------------------------------- #
