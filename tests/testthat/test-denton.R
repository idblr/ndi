context('denton')

# --------------- #
# denton testthat #
# --------------- #

test_that('denton throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    denton(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  expect_error(
    denton(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(
    denton(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    denton(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  expect_error(
    denton(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ref = 'terran',
      quiet = TRUE
    )
  )
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(
    denton(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  
})

test_that('denton works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(
    denton(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
      subgroup_ref = c('NHoLW', 'HoLW')
    )
  )
  
  expect_silent(
    denton(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  
  expect_silent(
    denton(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
      subgroup_ref = c('NHoLW', 'HoLW'),
      quiet = TRUE
    )
  )
  
})
