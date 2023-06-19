#--------------------------------------------------------------------------
# SET UP                          
#--------------------------------------------------------------------------
# Description -------------------------------------------------------------

# Derive ethnicity using CPRD Aurum & HES APC 
# ONS England & Wales 2011 Census - 18 categories
# Example code for deriving maternal ethnicity
# Suggested edits in CAPS if interested in a different age/sex

# Set working directory ---------------------------------------------------

setwd('X:/')

# Load packages -----------------------------------------------------------

pacman::p_load(tidyverse)

# Disable scientific notation ---------------------------------------------

options(scipen = 999)

#--------------------------------------------------------------------------
# DERIVE MATERNAL ETHNICITY IN AURUM                         
#--------------------------------------------------------------------------
# Identify ethnicity codes in observation files ---------------------------

# load ethnicity codelist 
ethnicity_codelist <- read_tsv('~/ethnicity_codes_publication.txt', col_names = TRUE)

# load clean Aurum observations file
load(file="~/obs_M_clean.Rdata")

# simplify codelist into medcode and 2011 census 18 categories (EDIT IF USING CENSUS 2001 OR 2021 CATEGORIES)
ethnicity_codelist <- ethnicity_codelist %>% 
  select(medcodeid, eth18_2011) %>% 
  mutate(eth18_2011 = as.factor(eth18_2011)) 

# identify obs with ethnicity codes
obs_M_ethnicity <- obs_M_clean %>% 
  inner_join(ethnicity_codelist, by = "medcodeid") %>% 
  select(M_patid, obsdate, eth18_2011)

rm(obs_M_clean, ethnicity_codelist)

# Drop and recode records based on dates ---------------------------------------------

# load patid list and aurum file with demographic data
load(file = "~/patids.Rdata")
load(file = "~/aurum_M_clean.Rdata")

# create ethnicity file which includes only mothers with at least 1 ethnicity record
# create date variables for child's fifth birthday (fifth birthday of their most recent child) and the date a year before mother's DOB
# EDIT IF DERIVING ETHNICITY FOR CHILDREN OR OTHER ADULTS
dates <- patids %>% 
  mutate(fifth_bday = deldate %m+% years(5)) %>% 
  select(-c(C_patid, deldate)) %>% 
  group_by(M_patid) %>% 
  slice(which.max(fifth_bday)) %>% 
  ungroup()

dates_2 <- aurum_M_clean %>% 
  select(M_patid, M_dob_aurum) %>% 
  mutate(year_before = M_dob_aurum %m-% years(1)) %>% 
  select(-M_dob_aurum)
  
obs_M_ethnicity <- obs_M_ethnicity %>% 
  left_join(dates, by = c("M_patid" = "M_patid")) %>% 
  left_join(dates_2, by = "M_patid")

# get denominators
n_ethnicity_records <- nrow(obs_M_ethnicity) 

# count how ethnicity records had implausible dates (>1 year before mum's dob given mum's parents could have nominated their ethnicity during pregnancy)
n_implausible_records <- sum(obs_M_ethnicity$obsdate < obs_M_ethnicity$year_before, na.rm = T) 
percent_implausible_records <- round(n_implausible_records/n_ethnicity_records*100, 2) 

# count how many ethnicity records had no obsdate
n_no_obsdate_records <- sum(is.na(obs_M_ethnicity$obsdate)) 
percent_no_obsdate_records <- round(n_no_obsdate_records/n_ethnicity_records*100, 2) 

# count how many ethnicity records had dates after their most recent child's fifth birthday
n_after_fifthbday_records <- sum(obs_M_ethnicity$obsdate > obs_M_ethnicity$fifth_bday, na.rm = T) 
percent_after_fifthbday_records <- round(n_after_fifthbday_records/n_ethnicity_records*100, 2) 

# recode implausible dates (>1 year before dob) to NA
# drop ethnicity records that were recorded after child's fifth birthday 
# don't exclude no obsdate. It gets processed later on via the final step of the algorithm
obs_M_ethnicity <- obs_M_ethnicity %>% 
  mutate(obsdate = replace(obsdate, (obsdate < year_before), NA)) %>% 
  filter(obsdate <= fifth_bday | is.na(obsdate)) %>%  
  select(-c(year_before, fifth_bday))

# count new n_ethnicity_records and missing obsdate after the above
n_ethnicity_records_excl <- nrow(obs_M_ethnicity) 
n_no_obsdate_records_excl <- sum(is.na(obs_M_ethnicity$obsdate)) 
percent_no_obsdate_records_excl <- round(n_no_obsdate_records_excl/n_ethnicity_records_excl*100, 2) 

# save
save(obs_M_ethnicity, file="~/obs_M_ethnicity.Rdata")

# Apply algorithm ---------------------------------------------------------

# SEE SUPPLEMENTAL FIGURE 2 in [citation]

# Remove duplicates with same ethnicity on the same date
obs_M_ethnicity_dedup <- obs_M_ethnicity %>% 
  distinct(eth18_2011, eth18_2011, eth19_2021, obsdate, M_patid, .keep_all = TRUE)

# remove 'not stated' records but keep them in a separate dataframe to add back in at the end
not_stated <- obs_M_ethnicity_dedup %>%
  filter(eth18_2011 == "Not stated")
obs_M_ethnicity_dedup <- obs_M_ethnicity_dedup %>%
  filter(eth18_2011 != "Not stated")

# select the most frequently occurring 
obs_M_ethnicity_dedup_freq <- obs_M_ethnicity_dedup %>% 
  group_by(M_patid, eth16_2001, eth18_2011, eth19_2021) %>%
  count() %>% 
  filter(n == max(n))

# find those who have achieved 1 unique record and separate them into a new dataframe
obs_M_ethnicity_unique_freq <- obs_M_ethnicity_dedup_freq %>% 
  group_by(M_patid) %>%
  count() %>% 
  filter(n == 1) %>% 
  select(-n) %>% 
  left_join(obs_M_ethnicity_dedup_freq, by = "M_patid") %>% 
  select(-n)

# check this new dataframe is all 1 row per patient
n_distinct(obs_M_ethnicity_unique_freq$M_patid) == nrow(obs_M_ethnicity_unique_freq)

# find patients who have not achieved 1 unique record (will expand out number of rows to attach obsdate again) and also separate them into a new dataframe
obs_M_ethnicity_nonunique <- obs_M_ethnicity_dedup_freq %>% 
  anti_join(obs_M_ethnicity_unique_freq, by = "M_patid") %>% 
  select(-n) %>% 
  left_join(obs_M_ethnicity_dedup, by = c("M_patid","eth18_2011"))

# check the 2 new dataframes still add up to the same number of patients as the original
n_distinct(obs_M_ethnicity_nonunique$M_patid) + n_distinct(obs_M_ethnicity_unique_freq$M_patid) == n_distinct(obs_M_ethnicity_dedup_freq$M_patid)

# for patients who have not achieved 1 unique record, select their most recent ethnicity record 
obs_M_ethnicity_unique_recent <- obs_M_ethnicity_nonunique %>%
  group_by(M_patid) %>%
  slice(which.max(obsdate))

# check whether all patients have been processed at this stage of the algorithm and whether there are any NA obsdates to deal with
sum(is.na(obs_M_ethnicity_unique_recent$obsdate)) # 0
nrow(obs_M_ethnicity_unique_recent) == n_distinct(obs_M_ethnicity_nonunique$M_patid) 
n_distinct(obs_M_ethnicity_nonunique$M_patid) - nrow(obs_M_ethnicity_unique_recent) 

# identify missing patients from obs_M_ethnicity_unique_recent, i.e. those who are proceeding to the next stage of the algorithm
obs_M_ethnicity_nonunique2 <- obs_M_ethnicity_nonunique %>% 
  anti_join(obs_M_ethnicity_unique_recent, by = "M_patid") 

# select category most frequently occurring in the general population - use 2011 as per CHIME but refine to subgroup 15-49 women in England as per representativeness analysis
# 2011 census: https://www.ethnicity-facts-figures.service.gov.uk/uk-population-by-ethnicity/national-and-regional-populations/population-of-england-and-wales/latest
# EDIT IF USING A DIFFERENT AGE GROUP/SEXX
genpop_ethnic_freq <- read_csv('~/e_ethnicity_2011_women_15to49y_freq_order.csv', col_names = TRUE)
genpop_ethnic_freq <- genpop_ethnic_freq %>% 
  mutate(eth18 = as.factor(eth18))
obs_M_ethnicity_unique_genpop <- obs_M_ethnicity_nonunique2 %>% 
  left_join(genpop_ethnic_freq, by=c(eth18_2011 = "eth18")) %>% 
  group_by(M_patid) %>% 
  slice(which.min(order)) %>% 
  select(-c(obsdate, order))

# check whether unqiue number of patients after algorithm applied is the same as before it was applied
n_distinct(obs_M_ethnicity_unique_recent$M_patid) + n_distinct(obs_M_ethnicity_unique_freq$M_patid) + n_distinct(obs_M_ethnicity_unique_genpop$M_patid) == n_distinct(obs_M_ethnicity_dedup$M_patid)

# join 
obs_M_ethnicity_unique <- obs_M_ethnicity_unique_recent %>% 
  select(-c(obsdate)) %>% 
  rbind(obs_M_ethnicity_unique_freq, obs_M_ethnicity_unique_genpop)

# find those who had 'Not stated' records and no other ethnicity record
not_stated <- not_stated %>%
  group_by(M_patid) %>% 
  summarise(M_patid = first(M_patid), eth16_2001 = first(eth16_2001), 
            eth18_2011 = first(eth18_2011), eth19_2021 = first(eth19_2021)) %>% 
  filter(!(M_patid %in% obs_M_ethnicity_unique$M_patid))

# join this dataframe of 'Not stated' records to the records that have been processed through the algorithm
obs_M_ethnicity_unique <- rbind(obs_M_ethnicity_unique, not_stated)

# check whether the unqiue number of patients after running the algorithm is the same as the original dataframe
n_distinct(obs_M_ethnicity_unique$M_patid) == n_distinct(obs_M_ethnicity$M_patid)

# check the 'not stated' category has been retained in this final dataframe
levels(obs_M_ethnicity_unique$eth18_2011)

# save
save(obs_M_ethnicity_clean, file="~/obs_M_ethnicity_clean.Rdata")

# Join to patid list to create an Aurum ethnicity file for the whole cohort --------------------------------------------------

# load patid list
load(file = "~/patids.Rdata")

mumpatids <- patids %>% 
  select(M_patid, C_patid) %>% 

# join to obs_M_ethnicity_clean file
M_ethnicity_aurum_clean <- mumpatids %>% 
  left_join(obs_M_ethnicity_clean, by = "M_patid")

# IF NEEDED - create 'Missing' categories for disaggregated ethnicity (distinct from 'Not stated')
levels(M_ethnicity_aurum_clean$eth18_2011) <- c(levels(M_ethnicity_aurum_clean$eth18_2011),"Missing")
M_ethnicity_aurum_clean <- M_ethnicity_aurum_clean %>% 
  mutate(eth18_2011 = replace_na(eth18_2011, "Missing")) 

# IF NEEDED - rename ethnicity fields to distinguish child and mother
M_ethnicity_aurum_clean <- M_ethnicity_aurum_clean %>% 
  rename(M_eth18_2011 = eth18_2011) 

# check that there are no more NAs
summary(is.na(M_ethnicity_aurum_clean)) 

# check that the number of patients equals the patid list
n_distinct(M_ethnicity_aurum_clean$M_patid) == n_distinct(patids$M_patid)

# save
save(M_ethnicity_aurum_clean, file="~/M_ethnicity_aurum_clean.Rdata")

# Save counts and checks -------------------------------------------------------------

# join & save counts
n_ethnicity_records <- as.data.frame(n_ethnicity_records)
n_implausible_records <- as.data.frame(n_implausible_records)
percent_implausible_records <- as.data.frame(percent_implausible_records)
n_no_obsdate_records <- as.data.frame(n_no_obsdate_records)
percent_no_obsdate_records <- as.data.frame(percent_no_obsdate_records)
n_after_fifthbday_records <- as.data.frame(n_after_fifthbday_records)
percent_after_fifthbday_records <- as.data.frame(percent_after_fifthbday_records)
n_ethnicity_records_excl <- as.data.frame(n_ethnicity_records_excl)
n_no_obsdate_records_excl <- as.data.frame(n_no_obsdate_records_excl)
percent_no_obsdate_records_excl <- as.data.frame(percent_no_obsdate_records_excl)

M_ethnicity_aurum_cleaning_counts <- cbind(n_ethnicity_records, 
                                   n_implausible_records, percent_implausible_records, 
                                   n_no_obsdate_records, percent_no_obsdate_records,
                                   n_after_fifthbday_records, percent_after_fifthbday_records,
                                   n_ethnicity_records_excl, 
                                   n_no_obsdate_records_excl,percent_no_obsdate_records_excl)

write_csv(M_ethnicity_aurum_cleaning_counts, file = "~/M_ethnicity_aurum_cleaning_counts.csv")

rm(list=ls())


#--------------------------------------------------------------------------
# DERIVE MOTHER ETHNICITY IN HES APC                        
#--------------------------------------------------------------------------
# Identify individuals with missing and 'Not stated' codes in Aurum after applying first part of the algorithm -------

# load files
load(file="~/M_ethnicity_aurum_clean.Rdata")

# split into those who have and don't have known ethnicity in aurum
M_no_aurum_eth <- M_ethnicity_aurum_clean %>% 
  filter(M_eth18_2011 == "Not stated" | M_eth18_2011 == "Missing")

M_aurum_eth <- M_ethnicity_aurum_clean %>% 
  filter(!(M_eth18_2011 == "Not stated" | M_eth18_2011 == "Missing")) # ethnicity sorted for this group, no further actions needed

# Clean ethnicity codes in HES --------------------------------------------

# load files
load(file="~/hes_M_episode_eth_clean.Rdata")
load(file="~/aurum_M_clean.Rdata")
load(file = "~/patids.Rdata")

# select and create need variables 
hes_M_ethnos <- hes_M_episode_eth_clean %>% 
  select(M_patid, ethnos, epistart) 

# create date variables for child's fifth birthday (fifth birthday of their most recent child) and the date a year before mother's DOB
# EDIT IF DERIVING ETHNICITY FOR CHILDREN OR OTHER ADULTS
dates <- patids %>% 
  mutate(fifth_bday = deldate %m+% years(5)) %>% 
  select(-c(C_patid, deldate)) %>% 
  group_by(M_patid) %>% 
  slice(which.max(fifth_bday)) %>% 
  ungroup()

dates_2 <- aurum_M_clean %>% 
  select(M_patid, M_dob_aurum) %>% 
  mutate(year_before = M_dob_aurum %m-% years(1)) %>% 
  select(-M_dob_aurum)

hes_M_ethnos <- hes_M_ethnos %>% 
  left_join(dates, by = c("M_patid" = "M_patid")) %>% 
  left_join(dates_2, by = "M_patid")

# get denominators
n_hes_eth_records <- nrow(hes_M_ethnos)  

# count how ethnicity records had implausible dates (>1 year before DOB as per Aurum cleaning)
n_implausible_hes_records <- sum(hes_M_ethnos$epistart < hes_M_ethnos$year_before, na.rm = T) 
percent_implausible_hes_records <- round(n_implausible_hes_records/n_hes_eth_records*100, 3) 

# count how many ethnicity records had no obsdate
n_no_epistart_records <- sum(is.na(hes_M_ethnos$epistart)) 
percent_no_epistart_records <- round(n_no_epistart_records/n_hes_eth_records*100, 3)

# count how many ethnicity records had dates after the child's fifth birthday
n_after_fifthbday_hes_records <- sum(hes_M_ethnos$epistart > hes_M_ethnos$fifth_bday, na.rm = T) 
percent_after_fifthbday_hes_records <- round(n_after_fifthbday_hes_records/n_hes_eth_records*100, 2)

# recode implausible dates (>1 year before dob) to NA
# drop ethnicity records that were recorded after fifth birthday 
# don't exclude no epistart It gets processed later on via the final step of the algorithm
hes_M_ethnos <- hes_M_ethnos %>% 
  mutate(epistart = replace(epistart, (epistart < year_before), NA)) %>% 
  filter(epistart <= fifth_bday | is.na(epistart)) %>%  
  select(-c(year_before, fifth_bday))

# categorise ethnicity, including preprocessed ethnicity into the most commonly occurring ini the general population (for women 15-49 years) 
# EDIT IF USING A DIFFERENT AGE GROUP/SEXX
genpop_ethnic_freq <- read_csv('~/e_ethnicity_2011_women_15to49y_freq_order.csv', col_names = TRUE)
hes_M_ethnos_clean <- hes_M_ethnos %>% 
  mutate(M_hes_eth18_2011 = ethnos) %>% 
  mutate(M_hes_eth18_2011 = replace(M_hes_eth18_2011, M_hes_eth18_2011 == "Bl_Afric", "African")) %>% 
  mutate(M_hes_eth18_2011 = replace(M_hes_eth18_2011, M_hes_eth18_2011 == "Bl_Carib", "Caribbean")) %>% 
  mutate(M_hes_eth18_2011 = replace(M_hes_eth18_2011, M_hes_eth18_2011 == "Bl_Other", "Any other Black, African or Caribbean background")) %>% 
  mutate(M_hes_eth18_2011 = replace(M_hes_eth18_2011, M_hes_eth18_2011 == "White", "English, Welsh, Scottish, Northern Irish or British")) %>% 
  mutate(M_hes_eth18_2011 = replace(M_hes_eth18_2011, M_hes_eth18_2011 == "Oth_Asian", "Any other Asian background")) %>% 
  mutate(M_hes_eth18_2011 = replace(M_hes_eth18_2011, M_hes_eth18_2011 == "Mixed", "White and Black Caribbean")) %>% 
  mutate(M_hes_eth18_2011 = replace(M_hes_eth18_2011, M_hes_eth18_2011 == "Other", "Any other ethnic group")) %>% 
  mutate(M_hes_eth18_2011 = replace(M_hes_eth18_2011, M_hes_eth18_2011 == "Unknown", "Missing")) %>% 
  mutate(M_hes_eth18_2011 = as.factor(M_hes_eth18_2011))

save(hes_M_ethnos_clean, file="~/hes_M_ethnos_clean.Rdata")

# count new n_ethnicity_records and missing obsdate after the above
n_hes_eth_records_excl <- nrow(hes_M_ethnos_clean) 
n_no_epistart_records_excl <- sum(is.na(hes_M_ethnos_clean$epistart)) 
percent_no_epistart_records_excl <- round(n_no_epistart_records_excl/n_hes_eth_records_excl*100, 3) 

# Apply HES APC part of the algorithm -------------------------------------

# Remove duplicates with same ethnicity on the same date
hes_M_ethnos_dedup <- hes_M_ethnos_clean %>% 
  select(-ethnos) %>% 
  distinct(M_patid, epistart, M_hes_eth18_2011)

# join only to those who had not stated or missing records in aurum, and deduplicate by mother's patid
hes_M_ethnicity <- M_no_aurum_eth %>% 
  left_join(hes_M_ethnos_dedup, by = c("M_patid")) %>% 
  select(-c(C_patid)) %>% 
  distinct()

# remove 'missing' records and records that didn't join but keep them in a separate dataframe to add back in at the end
not_joined <- hes_M_ethnicity %>% 
  filter(is.na(M_hes_eth18_2011))
missing <- hes_M_ethnicity %>%
  filter(M_hes_eth18_2011 == "Missing")
hes_M_ethnicity_dedup <- hes_M_ethnicity %>%
  filter(M_hes_eth18_2011 != "Missing") %>% 
  filter(!is.na(M_hes_eth18_2011))

# select the most frequently occurring 
hes_M_ethnicity_dedup_freq <- hes_M_ethnicity_dedup %>% 
  group_by(M_patid, M_hes_eth18_2011) %>%
  count() %>% 
  filter(n == max(n))

# find those who have achieved 1 unique record and separate them into a new dataframe
hes_M_ethnicity_unique_freq <- hes_M_ethnicity_dedup_freq %>% 
  group_by(M_patid) %>%
  count() %>% 
  filter(n == 1) %>% 
  select(-n) %>% 
  left_join(hes_M_ethnicity_dedup_freq, by = "M_patid") %>% 
  select(-n)

# check this new dataframe is all 1 row per patient
n_distinct(hes_M_ethnicity_unique_freq$M_patid) == nrow(hes_M_ethnicity_unique_freq)

# find patients who have not achieved 1 unique record (will expand out number of rows to attach epistart again) and also separate them into a new dataframe
hes_M_ethnicity_nonunique <- hes_M_ethnicity_dedup_freq %>% 
  anti_join(hes_M_ethnicity_unique_freq, by = "M_patid") %>% 
  select(-n) %>% 
  left_join(hes_M_ethnicity_dedup, by = c("M_patid", "M_hes_eth18_2011"))

# check the 2 new dataframes still add up to the same number of patients as the original
n_distinct(hes_M_ethnicity_nonunique$M_patid) + n_distinct(hes_M_ethnicity_unique_freq$M_patid) == n_distinct(hes_M_ethnicity_dedup_freq$M_patid)

# for patients who have not achieved 1 unique record, select their most recent ethnicity record 
hes_M_ethnicity_unique_recent <- hes_M_ethnicity_nonunique %>%
  group_by(M_patid) %>%
  slice(which.max(epistart))

# check whether all patients have been processed at this stage of the algorithm and whether there are any NA obsdates to deal with
sum(is.na(hes_M_ethnicity_unique_recent$epistart)) # 0
nrow(hes_M_ethnicity_unique_recent) == n_distinct(hes_M_ethnicity_nonunique$M_patid) 
n_distinct(hes_M_ethnicity_nonunique$M_patid) - nrow(hes_M_ethnicity_unique_recent) 

# identify missing patients from hes_M_ethnicity_unique_recent, i.e. those who are proceeding to the next stage of the algorithm
hes_M_ethnicity_nonunique2 <- hes_M_ethnicity_nonunique %>% 
  anti_join(hes_M_ethnicity_unique_recent, by = "M_patid") 

# select category most frequently occurring in the general population 
# EDIT IF USING A DIFFERENT AGE GROUP/SEXX
genpop_ethnic_freq <- read_csv('~/e_ethnicity_2011_women_15to49y_freq_order.csv', col_names = TRUE) 
genpop_ethnic_freq <- genpop_ethniM_freq %>% 
  mutate(eth18 = as.factor(eth18))
hes_M_ethnicity_unique_genpop <- hes_M_ethnicity_nonunique2 %>% 
  left_join(genpop_ethnic_freq, by=c(M_hes_eth18_2011 = "eth18")) %>% 
  group_by(M_patid) %>% 
  slice(which.min(order)) %>% 
  select(-c(epistart, order))

# check whether unqiue number of patients after algorithm applied is the same as before it was applied
n_distinct(hes_M_ethnicity_unique_recent$M_patid) + n_distinct(hes_M_ethnicity_unique_freq$M_patid) + n_distinct(hes_M_ethnicity_unique_genpop$M_patid) == n_distinct(hes_M_ethnicity_dedup$M_patid)

# join 
hes_M_ethnicity_unique <- hes_M_ethnicity_unique_recent %>% 
  select(-c(epistart)) %>% 
  rbind(hes_M_ethnicity_unique_freq)

# find those who had 'missing' records and no other ethnicity record
missing <- missing %>%
  group_by(M_patid) %>% 
  summarise(M_patid = first(M_patid), M_hes_eth18_2011 = first(M_hes_eth18_2011),
            M_eth16_2001 = first(M_eth16_2001),
            M_eth18_2011 = first(M_eth18_2011),
            M_eth19_2021 = first(M_eth19_2021)) %>% 
  filter(!(M_patid %in% hes_M_ethnicity_unique$M_patid))

# deduplicate not joined records (one per mother instead of one per child)
not_joined <- not_joined %>% 
  select(-epistart) %>% 
  distinct(M_patid, M_hes_eth18_2011, .keep_all = T) 
  
# join this dataframe of 'Not stated' records to the records that have been processed through the algorithm
hes_M_ethnicity_unique <- rbind(hes_M_ethnicity_unique, missing, not_joined) %>% ungroup()

# check whether the unqiue number of patients after running the algorithm is the same as the original dataframe
n_distinct(hes_M_ethnicity_unique$M_patid) == n_distinct(hes_M_ethnicity$M_patid)

# replace aurum eth categories with hes ones
M_ethnicity_hes <- hes_M_ethnicity_unique %>% 
  ungroup() %>% 
  mutate(across(everything(), as.character)) %>% 
  mutate(M_eth18_2011 = ifelse(is.na(M_hes_eth18_2011), M_eth18_2011, M_hes_eth18_2011)) %>% 
  mutate(M_eth18_2011 = replace(M_eth18_2011, M_eth18_2011 == "Bangladesi", "Bangladeshi")) %>% # fix CPRD HES APC's typo
  select(-c(M_hes_eth18_2011)) 

# Create aggregated ethnic groups -----------------------------------------

# convert into factor variables
M_ethnicity_hes_clean <- M_ethnicity_hes %>%
  mutate(M_eth18_2011 = as.factor(M_eth18_2011)) 

# REMOVE IF NOT NEEDED
# join to children's patids (mothers' ethnicity will be repeated across all their children)
C_patids <- patids %>% 
  filter(M_patid %in% M_ethnicity_hes_clean$M_patid) %>% 
  select(-deldate)

M_ethnicity_hes_clean <- C_patids %>% 
  left_join(M_ethnicity_hes_clean, by = "M_patid")

nrow(M_no_aurum_eth) == nrow(M_ethnicity_hes_clean)

# join to mothers with aurum ethnicity
M_ethnicity_clean <- rbind(M_ethnicity_hes_clean, M_aurum_eth)

# add factor levels that aren't present because of pre-processed hes ethnos
levels(M_ethnicity_clean$M_eth18_2011) <-  c(levels(M_ethnicity_clean$M_eth18_2011),"Any other White background",  
                                             "Irish", "Gypsy or Irish Traveller", "White and Asian", 
                                             "White and Black African", "Any other Mixed or multiple ethnic background", 
                                             "Arab")

# relevel disaggregated ethnic groups
M_ethnicity_clean <- M_ethnicity_clean %>%
  mutate(M_eth18_2011 = fct_relevel(M_eth18_2011, "English, Welsh, Scottish, Northern Irish or British", "Irish", "Gypsy or Irish Traveller", "Any other White background",  
                                    "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background",
                                    "Caribbean", "African", "Any other Black, African or Caribbean background",
                                    "White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed or multiple ethnic background",
                                    "Arab", "Any other ethnic group", "Not stated")) 

# create 5 level 2011 categories
M_ethnicity_clean <- M_ethnicity_clean %>%
  mutate(M_eth5_2011 = M_eth18_2011)
levels(M_ethnicity_clean$M_eth5_2011) <-  c(levels(M_ethnicity_clean$M_eth5_2011),"White",  
                                            "Asian or Asian British", "Black, African, Caribbean or Black British", 
                                            "Mixed or multiple ethnic groups", "Other ethnic group", "Unknown")
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "English, Welsh, Scottish, Northern Irish or British"] <- "White"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Irish"] <- "White"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Gypsy or Irish Traveller"] <- "White"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Any other White background"] <- "White"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "White and Black Caribbean"] <- "Mixed or multiple ethnic groups"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "White and Black African"] <- "Mixed or multiple ethnic groups"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "White and Asian"] <- "Mixed or multiple ethnic groups"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Any other Mixed or multiple ethnic background"] <- "Mixed or multiple ethnic groups"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Indian"] <- "Asian or Asian British"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Pakistani"] <- "Asian or Asian British"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Bangladeshi"] <- "Asian or Asian British"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Chinese"] <- "Asian or Asian British"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Any other Asian background"] <- "Asian or Asian British"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "African"] <- "Black, African, Caribbean or Black British"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Caribbean"] <- "Black, African, Caribbean or Black British"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Any other Black, African or Caribbean background"] <- "Black, African, Caribbean or Black British"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Arab"] <- "Other ethnic group"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Any other ethnic group"] <- "Other ethnic group"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Not stated"] <- "Unknown"
M_ethnicity_clean$M_eth5_2011[M_ethnicity_clean$M_eth5_2011 == "Missing"] <- "Unknown"

# drop and check levels
M_ethnicity_clean <- droplevels(M_ethnicity_clean)
levels(M_ethnicity_clean$M_eth18_2011)
levels(M_ethnicity_clean$M_eth5_2011)

# check for NAs
summary(is.na(M_ethnicity_clean)) 

# save
save(M_ethnicity_clean, file="~/M_ethnicity_clean.Rdata")

# Set missing and not stated ethnicity to unknown -------------------------

levels(M_ethnicity_clean$M_eth18_2011) <-  c(levels(M_ethnicity_clean$M_eth18_2011),"Unknown") 

M_ethnicity_clean$M_eth18_2011[M_ethnicity_clean$M_eth18_2011 == "Missing"] <- "Unknown"
M_ethnicity_clean$M_eth18_2011[M_ethnicity_clean$M_eth18_2011 == "Not stated"] <- "Unknown"
M_ethnicity_clean <- droplevels(M_ethnicity_clean) 

# save
save(M_ethnicity_clean, file="~/M_ethnicity_clean.Rdata")

# Save counts and checks -------------------------------------------------------------

# join & save counts
n_hes_eth_records <- as.data.frame(n_hes_eth_records)
n_implausible_hes_records <- as.data.frame(n_implausible_hes_records)
percent_implausible_hes_records <- as.data.frame(percent_implausible_hes_records)
n_no_epistart_records <- as.data.frame(n_no_epistart_records)
percent_no_epistart_records <- as.data.frame(percent_no_epistart_records)
n_after_fifthbday_hes_records <- as.data.frame(n_after_fifthbday_hes_records)
percent_after_fifthbday_hes_records <- as.data.frame(percent_after_fifthbday_hes_records)
n_hes_eth_records_excl <- as.data.frame(n_hes_eth_records_excl)
n_no_epistart_records_excl <- as.data.frame(n_no_epistart_records_excl)
percent_no_epistart_records_excl <- as.data.frame(percent_no_epistart_records_excl)

M_ethnicity_hes_cleaning_counts <- cbind(n_hes_eth_records, 
                                         n_implausible_hes_records, percent_implausible_hes_records, 
                                         n_no_epistart_records, percent_no_epistart_records,
                                         n_after_fifthbday_hes_records, percent_after_fifthbday_hes_records,
                                         n_hes_eth_records_excl, 
                                         n_no_epistart_records_excl,percent_no_epistart_records_excl)

write_csv(M_ethnicity_hes_cleaning_counts, file = "~/M_ethnicity_hes_cleaning_counts.csv")

rm(list=ls())
