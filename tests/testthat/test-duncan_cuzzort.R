context('duncan_cuzzort')

# ----------------------- #
# duncan_cuzzort testthat #
# ----------------------- #

test_that('duncan_cuzzort throws error with invalid arguments', {
  # Unavailable geography
  expect_error(
    duncan_cuzzort(
      geo_small = 'zcta',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  expect_error(
    duncan_cuzzort(
      geo_large = 'block group',
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable year
  expect_error(
    duncan_cuzzort(
      state = 'DC',
      year = 2005,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  # Unavailable subgroup
  expect_error(
    duncan_cuzzort(
      state = 'DC',
      year = 2020,
      subgroup = 'terran',
      quiet = TRUE
    )
  )
  
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  # Incorrect state
  expect_error(
    duncan_cuzzort(
      state = 'AB',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
})

test_that('duncan_cuzzort works', {
  skip_if(Sys.getenv('CENSUS_API_KEY') == '')
  
  expect_silent(
    duncan_cuzzort(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
    )
  )
  
  expect_silent(
    duncan_cuzzort(
      state = 'DC',
      year = 2020,
      subgroup = 'NHoLB',
      quiet = TRUE
    )
  )
  
  expect_silent(
    duncan_cuzzort(
      state = 'DC',
      year = 2020,
      subgroup = c('NHoLB', 'HoLB'),
      quiet = TRUE
    )
  )
  
})
