context('morgan_massey')

# ---------------------- #
# morgan_massey testthat #
# ---------------------- #

test_that('morgan_massey throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    morgan_massey(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  expect_error(
    morgan_massey(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(morgan_massey(
    state = 'DC',
    year = 2005,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
  # Unavailable subgroup
  expect_error(morgan_massey(
    state = 'DC',
    year = 2020,
    subgroup = 'terran',
    quiet = TRUE
  ))
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(morgan_massey(
    state = 'AB',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
})

test_that('morgan_massey works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(morgan_massey(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB')
  ))
  
  expect_silent(morgan_massey(
    state = 'DC',
    year = 2020,
    subgroup = 'NHoLB',
    quiet = TRUE
  ))
  
  expect_silent(morgan_massey(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
    quiet = TRUE
  ))
  
})
