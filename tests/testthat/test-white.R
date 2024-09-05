context('white')

# -------------- #
# white testthat #
# -------------- #

test_that('white throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    white(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  expect_error(
    white(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(
    white(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    white(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      quiet = TRUE
    )
  )
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(
    white(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
})

test_that('white works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(
    white(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB')
    )
  )
  
  expect_silent(
    white(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  expect_silent(
    white(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
      quiet = TRUE
    )
  )
  
})
