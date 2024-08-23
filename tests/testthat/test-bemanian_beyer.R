context('bemanian_beyer')

# ----------------------- #
# bemanian_beyer testthat #
# ----------------------- #

test_that('bemanian_beyer throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    bemanian_beyer(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  expect_error(
    bemanian_beyer(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(
    bemanian_beyer(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    bemanian_beyer(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  expect_error(
    bemanian_beyer(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'terran',
      quiet = TRUE
    )
  )
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(
    bemanian_beyer(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  
})

test_that('bemanian_beyer works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_warning(bemanian_beyer(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
    subgroup_ixn = c('NHoLW', 'HoLW')
  ))
  
  expect_warning(
    bemanian_beyer(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  
  expect_warning(bemanian_beyer(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
    subgroup_ixn = c('NHoLW', 'HoLW'),
    quiet = TRUE
  ))
  
})
