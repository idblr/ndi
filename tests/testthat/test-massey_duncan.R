context('massey_duncan')

# ---------------------- #
# massey_duncan testthat #
# ---------------------- #

test_that('massey_duncan throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    massey_duncan(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  expect_error(
    massey_duncan(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(
    massey_duncan(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    massey_duncan(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      quiet = TRUE
    )
  )
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(
    massey_duncan(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
})

test_that('massey_duncan works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(massey_duncan(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
  ))
  
  expect_silent(
    massey_duncan(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  expect_silent(massey_duncan(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
    quiet = TRUE
  ))
  
})
