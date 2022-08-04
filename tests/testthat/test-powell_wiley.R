context("powell_wiley")

####################
# powell_wiley testthat #
####################

test_that("powell_wiley throws error with invalid arguments", {
  
  # Unavailable geography
  expect_error(powell_wiley(geo = "zcta", state = "DC", year = 2009, quiet = TRUE))
  
  # Unavailable year
  expect_error(powell_wiley(state = "DC", year = 2009, quiet = TRUE))
  
  skip_if(Sys.getenv("CENSUS_API_KEY") == "")
  
  # Incorrect state
  expect_error(powell_wiley(state = "AB", year = 2020, quiet = TRUE))
  
  # Unavailable geography for DC (only 1 'county' in DC so, alone, NDI cannot be computed)
  expect_error(powell_wiley(geo = "county", state = "DC", year = 2009, quiet = TRUE))

}
)   

test_that("powell_wiley works", {  
  
  skip_if(Sys.getenv("CENSUS_API_KEY") == "")
  
  expect_message(powell_wiley(state = "DC", year = 2020))
  
  expect_message(powell_wiley(state = "DC", year = 2020, imp = TRUE))
  
  expect_silent(powell_wiley(state = "DC", year = 2020, quiet = TRUE))
  
  expect_silent(powell_wiley(state = "DC", year = 2020, imp = TRUE, quiet = TRUE))
  
}
)  
