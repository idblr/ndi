context("bravo")

##################
# bravo testthat #
##################

test_that("bravo throws error with invalid arguments", {
  
  # Unavailable geography
  expect_error(bravo(geo = "zcta", state = "DC", year = 2020, subgroup = "LtHS", quiet = TRUE))
  
  # Unavailable year
  expect_error(bravo(state = "DC", year = 2005, subgroup = "LtHS", quiet = TRUE))
  
  # Unavailable subgroup
  expect_error(bravo(state = "DC", year = 2020, subgroup = "terran", quiet = TRUE))
  
  skip_if(Sys.getenv("CENSUS_API_KEY") == "")
  
  # Incorrect state
  expect_error(bravo(state = "AB", year = 2020, subgroup = "LtHS", quiet = TRUE))
  
}
)   

test_that("bravo works", {  
  
  skip_if(Sys.getenv("CENSUS_API_KEY") == "")
  
  expect_silent(bravo(state = "DC", year = 2009, subgroup = c("LtHS", "HSGiE"))) 
  
  expect_silent(bravo(state = "DC", year = 2020, subgroup = c("LtHS", "HSGiE"))) 
  
  expect_silent(bravo(state = "DC", year = 2020, subgroup = "LtHS", quiet = TRUE))
  
  expect_silent(bravo(state = "DC", year = 2020, subgroup = c("LtHS", "HSGiE"), quiet = TRUE)) 
  
}
)  
