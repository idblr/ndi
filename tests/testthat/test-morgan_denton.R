context('morgan_denton')

# ---------------------- #
# morgan_denton testthat #
# ---------------------- #

test_that('morgan_denton throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    morgan_denton(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  expect_error(
    morgan_denton(
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
    morgan_denton(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    morgan_denton(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  expect_error(
    morgan_denton(
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
    morgan_denton(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  
})

test_that('morgan_denton works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(
    morgan_denton(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
      subgroup_ixn = c('NHoLW', 'HoLW')
    )
  )
  
  expect_silent(
    morgan_denton(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  
  expect_silent(
    morgan_denton(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
      subgroup_ixn = c('NHoLW', 'HoLW'),
      quiet = TRUE
    )
  )
  
})
