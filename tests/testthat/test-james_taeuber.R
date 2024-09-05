context('james_taeuber')

# ---------------------- #
# james_taeuber testthat #
# ---------------------- #

test_that('james_taeuber throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    james_taeuber(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  expect_error(
    james_taeuber(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(
    james_taeuber(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    james_taeuber(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      quiet = TRUE
    )
  )
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(
    james_taeuber(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
})

test_that('james_taeuber works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(
    james_taeuber(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
    )
  )
  
  expect_silent(
    james_taeuber(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  expect_silent(
    james_taeuber(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
      quiet = TRUE
    )
  )
  
})
