# ----------------------------------------------------------------------------------------------- #
# Code to prepare `DCtracts2020`
# ----------------------------------------------------------------------------------------------- #
#
# Created by: Ian Buller, Ph.D., M.A. (GitHub: @idblr)
# Created on: 2022-07-23
#
# Recently modified by: @idblr
# Recently modified on: 2024-07-06
#
# Notes:
# A) 2024-07-06 (@idblr): Re-formatted
# ----------------------------------------------------------------------------------------------- #

# ------------------ #
# NECESSARY PACKAGES #
# ------------------ #

loadedPackages <- c('dplyr', 'tidycensus', 'usethis')
suppressMessages(invisible(lapply(loadedPackages, library, character.only = TRUE)))

# -------- #
# SETTINGS #
# -------- #

## Access Key for census data download
### Obtain one at http://api.census.gov/data/key_signup.html
census_api_key('...') # INSERT YOUR OWN KEY FROM U.S. CENSUS API

# ---------------- #
# DATA PREPARATION #
# ---------------- #

# U.S. Census Bureau American Community Survey (ACS) 5-year variables
## For NDI (Messer)
### ACS-5 variables
messer_vars <- c(
  PctMenMgmtBusScArti_num1 = 'C24030_018',
  PctMenMgmtBusScArti_num2 = 'C24030_019',
  PctMenMgmtBusScArti_den = 'C24030_002',
  PctCrwdHH_num1 = 'B25014_005',
  PctCrwdHH_num2 = 'B25014_006',
  PctCrwdHH_num3 = 'B25014_007',
  PctCrwdHH_num4 = 'B25014_011',
  PctCrwdHH_num5 = 'B25014_012',
  PctCrwdHH_num6 = 'B25014_013',
  PctCrwdHH_den = 'B25014_001',
  PctHHPov_num = 'B17017_002',
  PctHHPov_den = 'B17017_001',
  PctFemHeadKids_num1 = 'B25115_012',
  PctFemHeadKids_num2 = 'B25115_025',
  PctFemHeadKids_den = 'B25115_001',
  PctPubAsst_num = 'B19058_002',
  PctPubAsst_den = 'B19058_001',
  PctHHUnder30K_num1 = 'B19001_002',
  PctHHUnder30K_num2 = 'B19001_003',
  PctHHUnder30K_num3 = 'B19001_004',
  PctHHUnder30K_num4 = 'B19001_005',
  PctHHUnder30K_num5 = 'B19001_006',
  PctHHUnder30K_den = 'B19001_001',
  PctEducLessThanHS_num = 'B06009_002',
  PctEducLessThanHS_den = 'B06009_001',
  PctUnemp_num = 'B23025_005',
  PctUnemp_den = 'B23025_003'
)

### Obtain ACS-5 data for DC tracts in 2020
DCtracts2020messer <- get_acs(
  geography = 'tract',
  year = 2020,
  output = 'wide',
  variables = messer_vars,
  state = 'DC'
)
### Format ACS-5 data for NDI (Messer) of DC tracts in 2020
DCtracts2020messer <- DCtracts2020messer[ , -2] # omit NAME feature (column)
DCtracts2020messer <- DCtracts2020messer %>%
  mutate(
    OCC = (PctMenMgmtBusScArti_num1E + PctMenMgmtBusScArti_num2E) / PctMenMgmtBusScArti_denE,
    CWD = (
      PctCrwdHH_num1E + PctCrwdHH_num2E + PctCrwdHH_num3E +
        PctCrwdHH_num4E + PctCrwdHH_num5E + PctCrwdHH_num6E
    ) / PctCrwdHH_denE,
    POV = PctHHPov_numE / PctHHPov_denE,
    FHH = (PctFemHeadKids_num1E + PctFemHeadKids_num2E) / PctFemHeadKids_denE,
    PUB = PctPubAsst_numE / PctPubAsst_denE,
    U30 = (
      PctHHUnder30K_num1E + PctHHUnder30K_num2E + PctHHUnder30K_num3E +
        PctHHUnder30K_num4E + PctHHUnder30K_num5E
    ) / PctHHUnder30K_denE,
    EDU = PctEducLessThanHS_numE / PctEducLessThanHS_denE,
    EMP = PctUnemp_numE / PctUnemp_denE
  )

### Clean-up and format
DCtracts2020messer <- DCtracts2020messer %>%
  select(GEOID, OCC, CWD, POV, FHH, PUB, U30, EDU, EMP)

## For NDI (Powell-Wiley)
### ACS-5 variables
powell_wiley_vars <- c(
  MedHHInc = 'B19013_001',
  PctRecvIDR_num = 'B19054_002',
  PctRecvIDR_den = 'B19054_001',
  PctPubAsst_num = 'B19058_002',
  PctPubAsst_den = 'B19058_001',
  MedHomeVal = 'B25077_001',
  PctMgmtBusScArti_num = 'C24060_002',
  PctMgmtBusScArti_den = 'C24060_001',
  PctFemHeadKids_num1 = 'B11005_007',
  PctFemHeadKids_num2 = 'B11005_010',
  PctFemHeadKids_den = 'B11005_001',
  PctOwnerOcc = 'DP04_0046P',
  PctNoPhone = 'DP04_0075P',
  PctNComPlmb = 'DP04_0073P',
  PctEduc_num25upHS = 'S1501_C01_009',
  PctEduc_num25upSC = 'S1501_C01_010',
  PctEduc_num25upAD = 'S1501_C01_011',
  PctEduc_num25upBD = 'S1501_C01_012',
  PctEduc_num25upGD = 'S1501_C01_013',
  PctEduc_den25up = 'S1501_C01_006',
  PctFamBelowPov = 'S1702_C02_001',
  PctUnempl = 'S2301_C04_001',
  TotalPopulation = 'B01001_001'
)

### Obtain ACS-5 data for DC tracts in 2020
DCtracts2020pw <- get_acs(
  geography = 'tract',
  year = 2020,
  output = 'wide',
  variables = powell_wiley_vars,
  state = 'DC'
)

### Format ACS-5 data for NDI (Powell-Wiley) of DC tracts in 2020
DCtracts2020pw <- DCtracts2020pw[,-2] # omit NAME feature (column)
DCtracts2020pw <- DCtracts2020pw %>%
  mutate(
    MedHHInc = MedHHIncE,
    PctRecvIDR = PctRecvIDR_numE / PctRecvIDR_denE * 100,
    PctPubAsst = PctPubAsst_numE / PctPubAsst_denE * 100,
    MedHomeVal = MedHomeValE,
    PctMgmtBusScArti = PctMgmtBusScArti_numE / PctMgmtBusScArti_denE * 100,
    PctFemHeadKids = (PctFemHeadKids_num1E + PctFemHeadKids_num2E) / PctFemHeadKids_denE * 100,
    PctOwnerOcc = PctOwnerOccE,
    PctNoPhone = PctNoPhoneE,
    PctNComPlmb = PctNComPlmbE,
    PctEducHSPlus = (
      PctEduc_num25upHSE + PctEduc_num25upSCE + PctEduc_num25upADE +
        PctEduc_num25upBDE + PctEduc_num25upGDE
    ) / PctEduc_den25upE * 100,
    PctEducBchPlus = (PctEduc_num25upBDE + PctEduc_num25upGDE) / PctEduc_den25upE * 100,
    PctFamBelowPov = PctFamBelowPovE,
    PctUnempl = PctUnemplE,
    TotalPop = TotalPopulationE
  ) %>%
  # Log transform median household income and median home value
  # Reverse code percentages so that higher values represent more deprivation
  # Round percentages to 1 decimal place
  mutate(
    logMedHHInc = log(MedHHInc),
    logMedHomeVal = log(MedHomeVal),
    PctNoIDR = 100 - PctRecvIDR,
    PctWorkClass = 100 - PctMgmtBusScArti,
    PctNotOwnerOcc = 100 - PctOwnerOcc,
    PctEducLTHS = 100 - PctEducHSPlus,
    PctEducLTBch = 100 - PctEducBchPlus
  ) %>%
  # Z-standardize the percentages
  mutate(
    PctNoIDRZ = scale(PctNoIDR),
    PctPubAsstZ = scale(PctPubAsst),
    PctWorkClassZ = scale(PctWorkClass),
    PctFemHeadKidsZ = scale(PctFemHeadKids),
    PctNotOwnerOccZ = scale(PctNotOwnerOcc),
    PctNoPhoneZ = scale(PctNoPhone),
    PctNComPlmbZ = scale(PctNComPlmb),
    PctEducLTHSZ = scale(PctEducLTHS),
    PctEducLTBchZ = scale(PctEducLTBch),
    PctFamBelowPovZ = scale(PctFamBelowPov),
    PctUnemplZ = scale(PctUnempl)
  )

### Clean-up and format
DCtracts2020pw <- DCtracts2020pw %>%
  select(
    GEOID,
    TotalPop,
    logMedHHInc,
    PctNoIDRZ,
    PctPubAsstZ,
    logMedHomeVal,
    PctWorkClassZ,
    PctFemHeadKidsZ,
    PctNotOwnerOccZ,
    PctNoPhoneZ,
    PctNComPlmbZ,
    PctEducLTHSZ,
    PctEducLTBchZ,
    PctFamBelowPovZ,
    PctUnemplZ
  )

# Combine
DCtracts2020 <- left_join(DCtracts2020messer, DCtracts2020pw, by = 'GEOID')
# reorder so TotalPop is second feature (column)
DCtracts2020 <- DCtracts2020[, c(1, 10, 2:9, 11:ncol(DCtracts2020))]

# ---------------- #
# DATA EXPORTATION #
# ---------------- #

use_data(DCtracts2020, overwrite = TRUE)

# ----------------------------------------- END OF CODE ----------------------------------------- #
