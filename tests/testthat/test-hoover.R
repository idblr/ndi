context('hoover')

# --------------- #
# hoover testthat #
# --------------- #

test_that('hoover throws error with invalid arguments', {
  # Unavailable geography
  expect_error(hoover(
    geo_small = 'zcta',
    state = 'DC',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  expect_error(
    hoover(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(hoover(
    state = 'DC',
    year = 2005,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
  # Unavailable subgroup
  expect_error(hoover(
    state = 'DC',
    year = 2020,
    subgroup = 'terran',
    quiet = TRUE
  ))
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(hoover(
    state = 'AB',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
})

test_that('hoover works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(hoover(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB')
  ))
  
  expect_silent(hoover(
    state = 'DC',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
  expect_silent(hoover(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
    quiet = TRUE
  ))
  
})
