context('duncan_duncan')

# ---------------------- #
# duncan_duncan testthat #
# ---------------------- #

test_that('duncan_duncan throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    duncan_duncan(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  expect_error(
    duncan_duncan(
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
    duncan_duncan(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    duncan_duncan(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  expect_error(
    duncan_duncan(
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
    duncan_duncan(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  
})

test_that('duncan_duncan works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(duncan_duncan(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
    subgroup_ref = c('NHoLW', 'HoLW')
  ))
  
  expect_silent(
    duncan_duncan(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      subgroup_ref = 'NHoLW',
      quiet = TRUE
    )
  )
  
  expect_silent(duncan_duncan(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
    subgroup_ref = c('NHoLW', 'HoLW'),
    quiet = TRUE
  ))
  
})
