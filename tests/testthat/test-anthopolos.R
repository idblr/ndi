context('anthopolos')

# ------------------- #
# anthopolos testthat #
# ------------------- #

test_that('anthopolos throws error with invalid arguments', {
  # Unavailable geography
  expect_error(anthopolos(
    geo = 'zcta',
    state = 'DC',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
  # Unavailable year
  expect_error(anthopolos(
    state = 'DC',
    year = 2005,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
  # Unavailable subgroup
  expect_error(anthopolos(
    state = 'DC',
    year = 2020,
    subgroup = 'terran',
    quiet = TRUE
  ))
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(anthopolos(
    state = 'AB',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
})

test_that('anthopolos works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_output(anthopolos(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB')
  ))
  
  expect_silent(anthopolos(
    state = 'DC',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
  expect_silent(anthopolos(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
    quiet = TRUE
  ))
  
})
