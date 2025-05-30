context('bell')

# ------------- #
# bell testthat #
# ------------- #

test_that('bell throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    bell(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  expect_error(
    bell(
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
    bell(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    bell(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  expect_error(
    bell(
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
    bell(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  
})

test_that('bell works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(
    bell(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
      subgroup_ixn = c('NHoLW', 'HoLW')
    )
  )
  
  expect_silent(
    bell(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ixn = 'NHoLW',
      quiet = TRUE
    )
  )
  
  expect_silent(
    bell(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
      subgroup_ixn = c('NHoLW', 'HoLW'),
      quiet = TRUE
    )
  )
  
})
