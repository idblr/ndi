context("krieger")

####################
# krieger testthat #
####################

test_that(" throws error with invalid arguments", {
  
  # Unavailable geography
  expect_error(krieger(geo = "zcta", state = "DC", year = 2020, quiet = TRUE))
  
  # Unavailable year
  expect_error(krieger(state = "DC", year = 2005, quiet = TRUE))
  
  skip_if(Sys.getenv("CENSUS_API_KEY") == "")
  
  # Incorrect state
  expect_error(krieger(state = "AB", year = 2020))
  
}
)   

test_that("krieger works", {  
  
  skip_if(Sys.getenv("CENSUS_API_KEY") == "")
  
  expect_message(krieger(state = "DC", year = 2020))
  
  expect_silent(krieger(state = "DC", year = 2020, quiet = TRUE))
  
  expect_silent(krieger(state = "DC", year = 2020, quiet = TRUE)) 
  
}
)  
