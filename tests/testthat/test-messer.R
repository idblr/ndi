context("messer")

###################
# messer testthat #
###################

test_that("messer throws error with invalid arguments", {
  
  # Unavailable geography
  expect_error(messer(geo = "zcta", state = "DC", year = 2020, quiet = TRUE))
  
  # Unavailable year
  expect_error(messer(state = "DC", year = 2005, quiet = TRUE))
  
  skip_if(Sys.getenv("CENSUS_API_KEY") == "")
  
  # Incorrect state
  expect_error(messer(state = "AB", year = 2020))
  
  # Unavailable geography for DC (only 1 'county' in DC so, alone, NDI cannot be computed)
  expect_error(messer(geo = "county", state = "DC", year = 2009, quiet = TRUE))
  
}
)   

test_that("messer works", {  
  
  skip_if(Sys.getenv("CENSUS_API_KEY") == "")
  
  expect_message(messer(state = "DC", year = 2020)) 
  
  expect_message(messer(state = "DC", year = 2020, imp = TRUE))
  
  expect_silent(messer(state = "DC", year = 2020, quiet = TRUE))
  
  expect_silent(messer(state = "DC", year = 2020, imp = TRUE, quiet = TRUE))
  
}
)  
