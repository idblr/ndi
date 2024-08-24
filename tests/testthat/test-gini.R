context('gini')

# ------------- #
# gini testthat #
# ------------- #

test_that('gini throws error with invalid arguments', {
  # Unavailable geography
  expect_error(gini(
    geo_small = 'zcta',
    state = 'DC',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  expect_error(
    gini(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(gini(
    state = 'DC',
    year = 2005,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
  # Unavailable subgroup
  expect_error(gini(
    state = 'DC',
    year = 2020,
    subgroup = 'terran',
    quiet = TRUE
  ))
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(gini(
    state = 'AB',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
})

test_that('gini works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_message(gini(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB')
  ))
  
  expect_silent(gini(
    state = 'DC',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
  expect_silent(gini(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
    quiet = TRUE
  ))
  
})
