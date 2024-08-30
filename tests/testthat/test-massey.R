context('massey')

# --------------- #
# massey testthat #
# --------------- #

test_that('massey throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    massey(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  expect_error(
    massey(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(
    massey(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    massey(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      quiet = TRUE
    )
  )
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(
    massey(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
})

test_that('massey works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(massey(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
  ))
  
  expect_silent(
    massey(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  expect_silent(massey(
    state = 'DC',
    year = 2020,
    subgroup = c('NHoLB', 'HoLB'),
    quiet = TRUE
  ))
  
})
