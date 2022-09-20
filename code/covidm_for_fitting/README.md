## Folder structure ##
- `~/code/`
  - `0_LoadAll.R`
    - call on everything that starts with `0_`
  - `0_1_util_functions.R`:
    - Load functions used cross the board
  - `0_2_StringencyIndex.R`:
    - Load NPI stringency data from oxcgrt
  - `0_3_Mobility.R`:
    - Functions that define the relationship between mobility and percentage 
    changes in the corresponding panels in the contact matrices
    - Load Google community mobility report for Thailand
    - General additive model that uses future NPI stringency to project for 
    possible changes in mobility and thus contacts
    - Definitions of school terms
  - `0_4_Vaccinations.R`
    - Load daily COVID-19 vaccination data from *Our World in Data*
    - Exploratory analysis on `people_fully_vaccinated` indicates that the 
    transition between primary course vaccination and booster campaign occurred
    between 2021-09-30 and 2021-10-17
    - Exploratory analysis on `daily_vaccinations` indicates that the daily 
    number of vaccinations may be okay to be set to 190000/day for starters
    - Definition of vaccine efficacies
  - `0_5_EpiParams.R`
    - Definitions of `cf` (clinical fraction) and `sus` (susceptibility)
    - Load contact matrices from [Prem et al. (2021)](https://github.com/kieshaprem/synthetic-contact-matrices/tree/master/output/syntheticcontactmatrices2020/overall)
  - `0_6_HealthCareSystem.R`
    - Definition of health care system processes: severe, critical, death
    - Dependent on the put from `0_4_Vaccinations.R` and thus must be ran after
    
