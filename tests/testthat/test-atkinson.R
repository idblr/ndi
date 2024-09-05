context('atkinson')

# ----------------- #
# atkinson testthat #
# ----------------- #

test_that('atkinson throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    atkinson(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  expect_error(
    atkinson(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(
    atkinson(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    atkinson(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      quiet = TRUE
    )
  )
  
  # Incorrect epsilon
  expect_error(
    atkinson(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      epsilon = 2,
      quiet = TRUE
    )
  )
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(
    atkinson(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
})

test_that('atkinson works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(
    atkinson(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB')
    )
  )
  
  expect_silent(
    atkinson(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
      holder = TRUE
    )
  )
  
  expect_silent(
    atkinson(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  expect_silent(
    atkinson(
      state = 'DC',
      year = 2020,
      subgroup = 'MedHHInc',
      quiet = TRUE
    )
  )
  
})
