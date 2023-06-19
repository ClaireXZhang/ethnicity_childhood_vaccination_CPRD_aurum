#--------------------------------------------------------------------------
# SET UP                          
#--------------------------------------------------------------------------
# Description -------------------------------------------------------------

# Derive routine childhood vaccinations in CPRD Aurum (except influenza)
# BEFORE RUNNING SCRIPT - turn each Excel sheet into a txt file

# Set working directory ---------------------------------------------------

setwd('X:/')

# Load packages -----------------------------------------------------------

pacman::p_load(tidyverse)

# Disable scientific notation ---------------------------------------------

options(scipen = 999)

#--------------------------------------------------------------------------
# CREATE BIRTHDAY COHORTS                     
#--------------------------------------------------------------------------

# count loss to follow-up  -----------------------

# load aurum source population file
load(file="~/sourcepop.Rdata")

# count how many children have follow-up end on or after 1, 2 and 5 years
# i.e. how many children remain in the cohort by 1, 2 and 5 years
total_child <- nrow(sourcepop)
fu_end_onafter1y <- sum(sourcepop$C_pc_fu_end >= (sourcepop$dob %m+% years(1)))
percent_fu_end_onafter1y <- round(fu_end_onafter1y/total_child*100, 2)  
fu_end_onafter2y <- sum(sourcepop$C_pc_fu_end >= (sourcepop$dob %m+% years(2)))
percent_fu_end_onafter2y <- round(fu_end_onafter2y/total_child*100, 2)  
fu_end_onafter5y <- sum(sourcepop$C_pc_fu_end >= (sourcepop$dob %m+% years(5)))
percent_fu_end_onafter5y <- round(fu_end_onafter5y/total_child*100, 2)  

# save counts
fu_end_onafter1y <- as.data.frame(fu_end_onafter1y)
percent_fu_end_onafter1y <- as.data.frame(percent_fu_end_onafter1y)
fu_end_onafter2y <- as.data.frame(fu_end_onafter2y)
percent_fu_end_onafter2y <- as.data.frame(percent_fu_end_onafter2y)
fu_end_onafter5y <- as.data.frame(fu_end_onafter5y)
percent_fu_end_onafter5y <- as.data.frame(percent_fu_end_onafter5y)

vax_bday_cohort_counts <- cbind(fu_end_onafter1y, percent_fu_end_onafter1y,
                       fu_end_onafter2y, percent_fu_end_onafter2y,
                       fu_end_onafter5y, percent_fu_end_onafter5y)

write_csv(vax_bday_cohort_counts, file = "~/vax_bday_cohort_counts.csv")

# create birthday cohorts -------------------------------------------------

# make first, second and fifth birthday variables
first_bday_patids <- sourcepop %>% 
  select(C_patid, dob, C_pc_fu_end) %>% 
  mutate(first_bday = dob %m+% years(1)) %>%
  filter(C_pc_fu_end >= first_bday)

second_bday_patids <- sourcepop %>% 
  select(C_patid, dob, C_pc_fu_end) %>% 
  mutate(second_bday = dob %m+% years(2)) %>%
  filter(C_pc_fu_end >= second_bday)

fifth_bday_patids <- sourcepop %>% 
  select(C_patid, dob, C_pc_fu_end) %>% 
  mutate(fifth_bday = dob %m+% years(5)) %>%
  filter(C_pc_fu_end >= fifth_bday)

# save
save(first_bday_patids, file = "~/first_bday_patids.Rdata")
save(second_bday_patids, file = "~/second_bday_patids.Rdata")
save(fifth_bday_patids, file = "~/fifth_bday_patids.Rdata")

#--------------------------------------------------------------------------
# CLEAN MMR (1y & 3y4m)                   
#--------------------------------------------------------------------------

# Identify all obs with medcodes and prodcodes  ---------------------------------------------

# load Aurum observations and drug issue files
load(file="~/obs_C_clean.Rdata")
load(file="~/drug_C_clean.Rdata")
load(file="~/patids.Rdata")

# medcode file 
mmr_medcodes <- read_tsv('~/mmr_medcodes.txt', 
                         col_names = TRUE) %>% 
  mutate(medcodeid = as.character(medcodeid)) %>% 
  select(medcodeid, type) 

# prodcode file
mmr_prodcodes <- read_tsv('~/mmr_prodcodes.txt', 
                          col_names = TRUE) %>% 
  mutate(prodcodeid = as.character(prodcodeid))  %>% 
  select(prodcodeid, type) 

# identify all medcode obs
obs_mmr_medcodes <- obs_C_clean %>% 
  inner_join(mmr_medcodes, by = "medcodeid") %>% 
  mutate(mmr = 1) %>% 
  mutate(code = "medcode") %>% 
  select(C_patid, obsdate, type, mmr, code, consid)

# identify all prodcode obs
obs_mmr_prodcodes <- drug_C_clean %>% 
  inner_join(mmr_prodcodes, by = "prodcodeid") %>% 
  mutate(mmr = 1) %>% 
  rename(obsdate = issuedate) %>% 
  mutate(code = "prodcode") %>% 
  mutate(consid = NA) %>% 
  select(C_patid, obsdate, type, mmr, code, consid)

# join medcodes and prodcodes
obs_mmr <- rbind(obs_mmr_medcodes, obs_mmr_prodcodes) 

# save for quality check analysis
save(obs_mmr, file="~/obs_mmr.Rdata")

# Counts & exclusions -----------------------------------------------------

# load C_patid list 
load(file = "~/sourcepop.Rdata")

babypatids <- sourcepop %>% 
  select(C_patid, dob, C_pc_fu_end) 

total_obs <- nrow(obs_mmr)

# count missing obsdates
nodate_obs <- sum(is.na(obs_mmr$obsdate)) 
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# count obsdates before DOB, after fifth birthday, and after follow up end
obs_mmr <- obs_mmr %>% 
  left_join(babypatids, by = "C_patid") %>% 
  mutate(fifth_bday = dob %m+% years(5)) 
  
obs_before_dob <- sum(obs_mmr$obsdate < obs_mmr$dob, na.rm = T)
percent_obs_before_dob <- round(obs_before_dob/total_obs*100, 3)

obs_after5y <- sum(obs_mmr$obsdate > obs_mmr$fifth_bday, na.rm = T)
percent_obs_after5y <- round(obs_after5y/total_obs*100, 3)

obs_after_fuend <- sum(obs_mmr$obsdate > obs_mmr$C_pc_fu_end, na.rm = T)
percent_obs_after_fuend<- round(obs_after_fuend/total_obs*100, 3)

# drop obs with no obsdate, obsdate before DOB, after fifth birthday and after follow up end
obs_mmr_full <- obs_mmr %>% 
  filter(!is.na(obsdate)) %>% 
  filter(obsdate >= dob) %>% 
  filter(obsdate <= fifth_bday) %>% 
  filter(obsdate <= C_pc_fu_end) 

# count number of med+prodcode obs after exclusions
obs_mmr_excl <- nrow(obs_mmr_full) 

# save
save(obs_mmr_full, file="~/obs_mmr_full.Rdata")

# Cut down to 1 record per person per obsdate  ------------------------------------

# prioritise full over partial codes, then medcode over prodcode
# if multiple of the same priority, prioritise ones that have consid
obs_mmr_full <- obs_mmr_full %>%
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, type == "full" & code == "prodcode", 2)) %>% 
  mutate(priority = replace(priority, type == "partial" & code == "medcode", 3)) %>% 
  mutate(priority = replace(priority, type == "partial" & code == "prodcode", 4)) %>% 
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  ungroup() %>% 
  select(-priority) %>% 
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, is.na(consid), 2)) %>%  
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  group_by(C_patid, obsdate) %>% 
  summarise(mmr = 1,
            consid = first(consid))
n_distinct(obs_mmr_full$C_patid) != nrow(obs_mmr_full) # should not yet be one row per patient

# Join to second bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/second_bday_patids.Rdata")

# join second bday
second_bday_mmr <- second_bday_patids %>% 
  left_join(obs_mmr_full, by = "C_patid") %>% 
  rename(obsdate_mmr = obsdate) %>% 
  mutate(mmr = replace_na(mmr, 0))

# count obs occurring after second birthday
obs_second_bday_mmr <- nrow(second_bday_mmr) # denom is number of mmr obs in this bday cohort 
mmr_obs_after2y <- sum(second_bday_mmr$obsdate_mmr > second_bday_mmr$second_bday, na.rm = T)
percent_mmr_obs_after2y <- round(mmr_obs_after2y/obs_second_bday_mmr*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
second_bday_mmr <- second_bday_mmr %>%
  mutate(mmr = replace(mmr, obsdate_mmr > second_bday, 0)) %>% 
  mutate(obsdate_mmr = replace(obsdate_mmr, obsdate_mmr > second_bday, NA)) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
second_bday_mmr <- second_bday_mmr %>% 
  distinct(C_patid, mmr, obsdate_mmr, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(second_bday_mmr$C_patid) == n_distinct(second_bday_patids$C_patid)

# save
save(second_bday_mmr, file="~/second_bday_mmr.Rdata")

# Join to fifth bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/fifth_bday_patids.Rdata")

# join fifth bday
fifth_bday_mmr <- fifth_bday_patids %>% 
  left_join(obs_mmr_full, by = "C_patid") %>% 
  rename(obsdate_mmr = obsdate) %>% 
  mutate(mmr = replace_na(mmr, 0))

obs_fifth_bday_mmr <- nrow(fifth_bday_mmr) # denom is number of mmr obs in this bday cohort 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
fifth_bday_mmr <- fifth_bday_mmr %>% 
  distinct(C_patid, mmr, obsdate_mmr, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(fifth_bday_mmr$C_patid) == n_distinct(fifth_bday_patids$C_patid)

# save
save(fifth_bday_mmr, file="~/fifth_bday_mmr.Rdata")

# Save counts and checks -------------------------------------------------------------

total_obs <- as.data.frame(total_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
obs_before_dob <- as.data.frame(obs_before_dob)
percent_obs_before_dob <- as.data.frame(percent_obs_before_dob)
obs_after5y <- as.data.frame(obs_after5y)
percent_obs_after5y <- as.data.frame(percent_obs_after5y)
obs_after_fuend <- as.data.frame(obs_after_fuend)
percent_obs_after_fuend <- as.data.frame(percent_obs_after_fuend)
obs_mmr_excl <- as.data.frame(obs_mmr_excl)

obs_second_bday_mmr <- as.data.frame(obs_second_bday_mmr)
mmr_obs_after2y <- as.data.frame(mmr_obs_after2y)
percent_mmr_obs_after2y <- as.data.frame(percent_mmr_obs_after2y)

obs_fifth_bday_mmr <- as.data.frame(obs_fifth_bday_mmr)

mmr_cleaning_counts <- cbind(total_obs, 
                             nodate_obs, percent_nodate_obs,
                             obs_before_dob, percent_obs_before_dob,
                             obs_after5y, percent_obs_after5y,
                             obs_after_fuend, percent_obs_after_fuend,
                             obs_mmr_excl,
                             obs_second_bday_mmr, mmr_obs_after2y, percent_mmr_obs_after2y,
                             obs_fifth_bday_mmr)

write_csv(mmr_cleaning_counts, file = "~/mmr_cleaning_counts.csv")

# Create completeness & timing variables for second bday cohort ---------------------------------------------------

# this is the cumulative number of doses expected for a given vaccine by each birthday, as per the national schedule
vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
second_bday_mmr_completeness <- second_bday_mmr %>% 
  group_by(C_patid) %>%
  summarise(mmr = sum(mmr))

# create expected dose counts
second_bday_mmr_completeness <- second_bday_mmr_completeness %>% 
  mutate(mmr_exp = as.numeric(vaccine_exp_dose$second_bday_exp[vaccine_exp_dose$vaccine == "mmr"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
second_bday_mmr_completeness$mmr_pc[second_bday_mmr_completeness$mmr_exp>second_bday_mmr_completeness$mmr] <- 0
second_bday_mmr_completeness$mmr_pc[second_bday_mmr_completeness$mmr_exp<=second_bday_mmr_completeness$mmr] <- 1

# check for NAs
sum(is.na(second_bday_mmr_completeness$mmr_pc))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
second_bday_mmr_analysis <- second_bday_mmr_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
second_bday_mmr_analysis <- second_bday_mmr_analysis %>% 
  mutate(second_bday = dob %m+% years(2)) %>%
  filter(second_bday >= "2006-04-01") %>% 
  filter(second_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (second_bday >= "2006-04-01" & second_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (second_bday >= "2007-04-01" & second_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (second_bday >= "2008-04-01" & second_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (second_bday >= "2009-04-01" & second_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (second_bday >= "2010-04-01" & second_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (second_bday >= "2011-04-01" & second_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (second_bday >= "2012-04-01" & second_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (second_bday >= "2013-04-01" & second_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (second_bday >= "2014-04-01" & second_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (second_bday >= "2015-04-01" & second_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (second_bday >= "2016-04-01" & second_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (second_bday >= "2017-04-01" & second_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (second_bday >= "2018-04-01" & second_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (second_bday >= "2019-04-01" & second_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (second_bday >= "2020-04-01" & second_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(second_bday_mmr_analysis, file="~/second_bday_mmr_analysis.Rdata")

# Create completeness & timing variables for fifth bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
fifth_bday_mmr_completeness <- fifth_bday_mmr %>% 
  group_by(C_patid) %>%
  summarise(mmr = sum(mmr))

# create expected dose counts
fifth_bday_mmr_completeness <- fifth_bday_mmr_completeness %>% 
  mutate(mmr_exp_pc = as.numeric(vaccine_exp_dose$second_bday_exp[vaccine_exp_dose$vaccine == "mmr"])) %>% 
  mutate(mmr_exp_booster = as.numeric(vaccine_exp_dose$fifth_bday_exp[vaccine_exp_dose$vaccine == "mmr"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
fifth_bday_mmr_completeness$mmr_pc[fifth_bday_mmr_completeness$mmr_exp_pc>fifth_bday_mmr_completeness$mmr] <- 0
fifth_bday_mmr_completeness$mmr_pc[fifth_bday_mmr_completeness$mmr_exp_pc<=fifth_bday_mmr_completeness$mmr] <- 1

# derive booster variable (complete = 1 vs incomplete = 0)
fifth_bday_mmr_completeness$mmr_booster[fifth_bday_mmr_completeness$mmr_exp_booster>fifth_bday_mmr_completeness$mmr] <- 0
fifth_bday_mmr_completeness$mmr_booster[fifth_bday_mmr_completeness$mmr_exp_booster<=fifth_bday_mmr_completeness$mmr] <- 1

# check for NAs
sum(is.na(fifth_bday_mmr_completeness$mmr_pc))
sum(is.na(fifth_bday_mmr_completeness$mmr_booster))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
fifth_bday_mmr_analysis <- fifth_bday_mmr_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
fifth_bday_mmr_analysis <- fifth_bday_mmr_analysis %>% 
  mutate(fifth_bday = dob %m+% years(5)) %>%
  filter(fifth_bday >= "2006-04-01") %>% 
  filter(fifth_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (fifth_bday >= "2006-04-01" & fifth_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2007-04-01" & fifth_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2008-04-01" & fifth_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2009-04-01" & fifth_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2010-04-01" & fifth_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2011-04-01" & fifth_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2012-04-01" & fifth_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2013-04-01" & fifth_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2014-04-01" & fifth_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2015-04-01" & fifth_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2016-04-01" & fifth_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2017-04-01" & fifth_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2018-04-01" & fifth_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2019-04-01" & fifth_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2020-04-01" & fifth_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(fifth_bday_mmr_analysis, file="~/fifth_bday_mmr_analysis.Rdata")

#--------------------------------------------------------------------------
# CLEAN ROTAVIRUS (8w & 12w since July 2013)                    
#--------------------------------------------------------------------------

# Identify all obs with medcodes and prodcodes  ---------------------------------------------

load(file="~/obs_C_clean.Rdata")
load(file="~/drug_C_clean.Rdata")
load(file="~/patids.Rdata")

# medcode file 
rotavirus_medcodes <- read_tsv('~/rotavirus_medcodes.txt', 
                         col_names = TRUE) %>% 
  mutate(medcodeid = as.character(medcodeid)) %>% 
  select(medcodeid)

# prodcode file
rotavirus_prodcodes <- read_tsv('~/rotavirus_prodcodes.txt', 
                          col_names = TRUE) %>% 
  mutate(prodcodeid = as.character(prodcodeid)) %>% 
  select(prodcodeid)

# identify all medcode obs
obs_rotavirus_medcodes <- obs_C_clean %>% 
  inner_join(rotavirus_medcodes, by = "medcodeid") %>% 
  mutate(rotavirus = 1) %>% 
  mutate(code = "medcode") %>% 
  select(C_patid, obsdate, rotavirus, code, consid)

# identify all prodcode obs
obs_rotavirus_prodcodes <- drug_C_clean %>% 
  inner_join(rotavirus_prodcodes, by = "prodcodeid") %>% 
  mutate(rotavirus = 1) %>% 
  rename(obsdate = issuedate) %>% 
  mutate(code = "prodcode") %>% 
  mutate(consid = NA) %>% 
  select(C_patid, obsdate, rotavirus, code, consid)

# join medcodes and prodcodes
obs_rotavirus <- rbind(obs_rotavirus_medcodes, obs_rotavirus_prodcodes) 

# save for quality check analysis
save(obs_rotavirus, file="~/obs_rotavirus.Rdata")

# Counts & exclusions -----------------------------------------------------

# load C_patid list 
load(file = "~/sourcepop.Rdata")

babypatids <- sourcepop %>% 
  select(C_patid, dob, C_pc_fu_end) 

total_obs <- nrow(obs_rotavirus)

# count missing obsdates
nodate_obs <- sum(is.na(obs_rotavirus$obsdate)) 
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# count obsdates before DOB, before July 2013, after fifth birthday, and after follow up end
obs_rotavirus <- obs_rotavirus %>% 
  left_join(babypatids, by = "C_patid") %>% 
  mutate(fifth_bday = dob %m+% years(5)) 

obs_before_dob <- sum(obs_rotavirus$obsdate < obs_rotavirus$dob, na.rm = T)
percent_obs_before_dob <- round(obs_before_dob/total_obs*100, 3)

obs_before_jul2013<- sum(obs_rotavirus$obsdate < "2013-07-01", na.rm = T)
percent_obs_before_jul2013 <- round(obs_before_jul2013/total_obs*100, 3)

obs_after5y <- sum(obs_rotavirus$obsdate > obs_rotavirus$fifth_bday, na.rm = T)
percent_obs_after5y <- round(obs_after5y/total_obs*100, 3)

obs_after_fuend <- sum(obs_rotavirus$obsdate > obs_rotavirus$C_pc_fu_end, na.rm = T)
percent_obs_after_fuend<- round(obs_after_fuend/total_obs*100, 3)

# drop obs with no obsdate, obsdate before DOB, before July 2013, after fifth birthday and after follow up end
obs_rotavirus_full <- obs_rotavirus %>% 
  filter(!is.na(obsdate)) %>% 
  filter(obsdate >= dob) %>% 
  filter(obsdate >= "2013-07-01") %>% 
  filter(obsdate <= fifth_bday) %>% 
  filter(obsdate <= C_pc_fu_end) 

# count number of med+prodcode obs after exclusions
obs_rotavirus_excl <- nrow(obs_rotavirus_full) 

# save 
save(obs_rotavirus_full, file="~/obs_rotavirus_full.Rdata")

# Cut down to 1 record per person per obsdate   ------------------------------------

# prioritise ones that have consid
obs_rotavirus_full <- obs_rotavirus_full %>%
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, is.na(consid), 2)) %>%  
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  group_by(C_patid, obsdate) %>% 
  summarise(rotavirus = first(rotavirus),
            consid = first(consid))
n_distinct(obs_rotavirus_full$C_patid) != nrow(obs_rotavirus_full) # should not yet be one row per patient

# Join to first bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/first_bday_patids.Rdata")

# join first bday for those born on or after July 2013
first_bday_patids_jul2013 <- first_bday_patids %>% 
  filter(dob >= "2013-07-01") 

first_bday_rotavirus <- first_bday_patids_jul2013 %>% 
  left_join(obs_rotavirus_full, by = "C_patid") %>% 
  rename(obsdate_rotavirus = obsdate) %>% 
  mutate(rotavirus = replace_na(rotavirus, 0))

# count obs occurring after first birthday
obs_first_bday_rotavirus <- nrow(first_bday_rotavirus) # denom is number of rotavirus obs in this bday cohort 
rotavirus_obs_after1y <- sum(first_bday_rotavirus$obsdate_rotavirus > first_bday_rotavirus$first_bday, na.rm = T)
percent_rotavirus_obs_after1y <- round(rotavirus_obs_after1y/obs_first_bday_rotavirus*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
first_bday_rotavirus <- first_bday_rotavirus %>%
  mutate(rotavirus = replace(rotavirus, obsdate_rotavirus > first_bday, 0)) %>% 
  mutate(obsdate_rotavirus = replace(obsdate_rotavirus, obsdate_rotavirus > first_bday, NA))  

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
first_bday_rotavirus <- first_bday_rotavirus %>% 
  distinct(C_patid, rotavirus, obsdate_rotavirus, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(first_bday_rotavirus$C_patid) == n_distinct(first_bday_patids$C_patid[first_bday_patids$dob >= "2013-07-01"])

# save
save(first_bday_rotavirus, file="~/first_bday_rotavirus.Rdata")

# Save counts and checks -------------------------------------------------------------

total_obs <- as.data.frame(total_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
obs_before_dob <- as.data.frame(obs_before_dob)
percent_obs_before_dob <- as.data.frame(percent_obs_before_dob)
obs_after5y <- as.data.frame(obs_after5y)
percent_obs_after5y <- as.data.frame(percent_obs_after5y)
obs_after_fuend <- as.data.frame(obs_after_fuend)
percent_obs_after_fuend <- as.data.frame(percent_obs_after_fuend)
obs_rotavirus_excl <- as.data.frame(obs_rotavirus_excl)

obs_first_bday_rotavirus <- as.data.frame(obs_first_bday_rotavirus)
rotavirus_obs_after1y <- as.data.frame(rotavirus_obs_after1y)
percent_rotavirus_obs_after1y <- as.data.frame(percent_rotavirus_obs_after1y)

rotavirus_cleaning_counts <- cbind(total_obs, 
                             nodate_obs, percent_nodate_obs,
                             obs_before_dob, percent_obs_before_dob,
                             obs_after5y, percent_obs_after5y,
                             obs_after_fuend, percent_obs_after_fuend,
                             obs_rotavirus_excl,
                             obs_first_bday_rotavirus, rotavirus_obs_after1y, percent_rotavirus_obs_after1y)

write_csv(rotavirus_cleaning_counts, file = "~/rotavirus_cleaning_counts.csv")

# Create completeness & timing variables for first bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
first_bday_rotavirus_completeness <- first_bday_rotavirus %>% 
  group_by(C_patid) %>%
  summarise(rotavirus = sum(rotavirus))

# create expected dose counts
first_bday_rotavirus_completeness <- first_bday_rotavirus_completeness %>% 
  mutate(rotavirus_exp = as.numeric(vaccine_exp_dose$first_bday_exp[vaccine_exp_dose$vaccine == "rotavirus"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
first_bday_rotavirus_completeness$rotavirus_pc[first_bday_rotavirus_completeness$rotavirus_exp>first_bday_rotavirus_completeness$rotavirus] <- 0
first_bday_rotavirus_completeness$rotavirus_pc[first_bday_rotavirus_completeness$rotavirus_exp<=first_bday_rotavirus_completeness$rotavirus] <- 1

# check for NAs
sum(is.na(first_bday_rotavirus_completeness$rotavirus_pc))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
first_bday_rotavirus_analysis <- first_bday_rotavirus_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
first_bday_rotavirus_analysis <- first_bday_rotavirus_analysis %>% 
  mutate(first_bday = dob %m+% years(1)) %>%
  filter(first_bday >= "2006-04-01") %>% 
  filter(first_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (first_bday >= "2006-04-01" & first_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (first_bday >= "2007-04-01" & first_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (first_bday >= "2008-04-01" & first_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (first_bday >= "2009-04-01" & first_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (first_bday >= "2010-04-01" & first_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (first_bday >= "2011-04-01" & first_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (first_bday >= "2012-04-01" & first_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (first_bday >= "2013-04-01" & first_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (first_bday >= "2014-04-01" & first_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (first_bday >= "2015-04-01" & first_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (first_bday >= "2016-04-01" & first_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (first_bday >= "2017-04-01" & first_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (first_bday >= "2018-04-01" & first_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (first_bday >= "2019-04-01" & first_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (first_bday >= "2020-04-01" & first_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(first_bday_rotavirus_analysis, file="~/first_bday_rotavirus_analysis.Rdata")

#--------------------------------------------------------------------------
# CLEAN MENB (8w, 16w & 1y since September 2015)                    
#--------------------------------------------------------------------------

# Identify all obs with medcodes and prodcodes  ---------------------------------------------

# load child observations and drug issue file
load(file="~/obs_C_clean.Rdata")
load(file="~/drug_C_clean.Rdata")
load(file="~/patids.Rdata")

# medcode file 
menb_medcodes <- read_tsv('~/menb_medcodes.txt', 
                               col_names = TRUE) %>% 
  mutate(medcodeid = as.character(medcodeid)) %>% 
  select(medcodeid, type)

# prodcode file
menb_prodcodes <- read_tsv('~/menb_prodcodes.txt', 
                                col_names = TRUE) %>% 
  mutate(prodcodeid = as.character(prodcodeid)) %>% 
  select(prodcodeid, type)

# identify all medcode obs
obs_menb_medcodes <- obs_C_clean %>% 
  inner_join(menb_medcodes, by = "medcodeid") %>% 
  mutate(menb = 1) %>% 
  mutate(code = "medcode") %>% 
  select(C_patid, obsdate, menb, code, consid, type)

# identify all prodcode obs
obs_menb_prodcodes <- drug_C_clean %>% 
  inner_join(menb_prodcodes, by = "prodcodeid") %>% 
  mutate(menb = 1) %>% 
  rename(obsdate = issuedate) %>% 
  mutate(code = "prodcode") %>% 
  mutate(consid = NA) %>% 
  select(C_patid, obsdate, menb, code, consid, type)

# join medcodes and prodcodes
obs_menb <- rbind(obs_menb_medcodes, obs_menb_prodcodes) 

# save for quality check analysis
save(obs_menb, file="~/obs_menb.Rdata")

# Counts & exclusions -----------------------------------------------------

# load C_patid list 
load(file = "~/sourcepop.Rdata")

babypatids <- sourcepop %>% 
  select(C_patid, dob, C_pc_fu_end) 

total_obs <- nrow(obs_menb)

# count missing obsdates
nodate_obs <- sum(is.na(obs_menb$obsdate))
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# count obsdates before DOB, before September 2015, after fifth birthday, and after follow up end
obs_menb <- obs_menb %>% 
  left_join(babypatids, by = "C_patid") %>% 
  mutate(fifth_bday = dob %m+% years(5)) 

obs_before_dob <- sum(obs_menb$obsdate < obs_menb$dob, na.rm = T)
percent_obs_before_dob <- round(obs_before_dob/total_obs*100, 3)

obs_before_sep2015<- sum(obs_menb$obsdate < "2015-09-01", na.rm = T)
percent_obs_before_sep2015 <- round(obs_before_sep2015/total_obs*100, 3)

obs_after5y <- sum(obs_menb$obsdate > obs_menb$fifth_bday, na.rm = T)
percent_obs_after5y <- round(obs_after5y/total_obs*100, 3)

obs_after_fuend <- sum(obs_menb$obsdate > obs_menb$C_pc_fu_end, na.rm = T)
percent_obs_after_fuend<- round(obs_after_fuend/total_obs*100, 3)

# drop obs with no obsdate, obsdate before DOB, before September 2015, after fifth birthday and after follow up end
obs_menb_full <- obs_menb %>% 
  filter(!is.na(obsdate)) %>% 
  filter(obsdate >= dob) %>% 
  filter(obsdate >= "2015-09-01") %>% 
  filter(obsdate <= fifth_bday) %>% 
  filter(obsdate <= C_pc_fu_end) 

# count number of med+prodcode obs after exclusions
obs_menb_excl <- nrow(obs_menb_full) 

# save
save(obs_menb_full, file="~/obs_menb_full.Rdata")

# Cut down to 1 record per person per obsdate   ------------------------------------

# prioritise ones that have consid 
obs_menb_full <- obs_menb_full %>%
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, is.na(consid), 2)) %>%  
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  group_by(C_patid, obsdate) %>% 
  summarise(menb = first(menb),
            consid = first(consid))
n_distinct(obs_menb_full$C_patid) != nrow(obs_menb_full) # should not yet be one row per patient

# Join to first bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/first_bday_patids.Rdata")

# join first bday for those born on or after September 2015
first_bday_patids_sep2015 <- first_bday_patids %>% 
  filter(dob >= "2015-09-01") 

first_bday_menb <- first_bday_patids_sep2015 %>% 
  left_join(obs_menb_full, by = "C_patid") %>% 
  rename(obsdate_menb = obsdate) %>% 
  mutate(menb = replace_na(menb, 0))

# count obs occurring after first birthday
obs_first_bday_menb <- nrow(first_bday_menb) # denom is number of menb obs in this bday cohort 
menb_obs_after1y <- sum(first_bday_menb$obsdate_menb > first_bday_menb$first_bday, na.rm = T)
percent_menb_obs_after1y <- round(menb_obs_after1y/obs_first_bday_menb*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
first_bday_menb <- first_bday_menb %>%
  mutate(menb = replace(menb, obsdate_menb > first_bday, 0)) %>% 
  mutate(obsdate_menb = replace(obsdate_menb, obsdate_menb > first_bday, NA)) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
first_bday_menb <- first_bday_menb %>% 
  distinct(C_patid, menb, obsdate_menb, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(first_bday_menb$C_patid) == n_distinct(first_bday_patids$C_patid[first_bday_patids$dob >= "2015-09-01"])

# save
save(first_bday_menb, file="~/first_bday_menb.Rdata")

# Join to second bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/second_bday_patids.Rdata")

# join second bday for those born on or after September 2015
second_bday_patids_sep2015 <- second_bday_patids %>% 
  filter(dob >= "2015-09-01") 

second_bday_menb <- second_bday_patids_sep2015 %>% 
  left_join(obs_menb_full, by = "C_patid") %>% 
  rename(obsdate_menb = obsdate) %>% 
  mutate(menb = replace_na(menb, 0))

# count obs occurring after second birthday
obs_second_bday_menb <- nrow(second_bday_menb) # denom is number of menb obs in this bday cohort 
menb_obs_after2y <- sum(second_bday_menb$obsdate_menb > second_bday_menb$second_bday, na.rm = T)
percent_menb_obs_after2y <- round(menb_obs_after2y/obs_second_bday_menb*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
second_bday_menb <- second_bday_menb %>%
  mutate(menb = replace(menb, obsdate_menb > second_bday, 0)) %>% 
  mutate(obsdate_menb = replace(obsdate_menb, obsdate_menb > second_bday, NA)) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
second_bday_menb <- second_bday_menb %>% 
  distinct(C_patid, menb, obsdate_menb, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(second_bday_menb$C_patid) == n_distinct(second_bday_patids$C_patid[second_bday_patids$dob >= "2015-09-01"])

# save
save(second_bday_menb, file="~/second_bday_menb.Rdata")

# Save counts and checks -------------------------------------------------------------

total_obs <- as.data.frame(total_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
obs_before_dob <- as.data.frame(obs_before_dob)
percent_obs_before_dob <- as.data.frame(percent_obs_before_dob)
obs_after5y <- as.data.frame(obs_after5y)
percent_obs_after5y <- as.data.frame(percent_obs_after5y)
obs_after_fuend <- as.data.frame(obs_after_fuend)
percent_obs_after_fuend <- as.data.frame(percent_obs_after_fuend)
obs_menb_excl <- as.data.frame(obs_menb_excl)

obs_first_bday_menb <- as.data.frame(obs_first_bday_menb)
menb_obs_after1y <- as.data.frame(menb_obs_after1y)
percent_menb_obs_after1y <- as.data.frame(percent_menb_obs_after1y)
obs_second_bday_menb <- as.data.frame(obs_second_bday_menb)
menb_obs_after2y <- as.data.frame(menb_obs_after2y)
percent_menb_obs_after2y <- as.data.frame(percent_menb_obs_after2y)

menb_cleaning_counts <- cbind(total_obs, 
                              nodate_obs, percent_nodate_obs,
                              obs_before_dob, percent_obs_before_dob,
                              obs_after5y, percent_obs_after5y,
                              obs_after_fuend, percent_obs_after_fuend,
                              obs_menb_excl,
                              obs_first_bday_menb, menb_obs_after1y, percent_menb_obs_after1y,
                              obs_second_bday_menb, menb_obs_after2y, percent_menb_obs_after2y)

write_csv(menb_cleaning_counts, file = "~/menb_cleaning_counts.csv")


# Create completeness & timing variables for first bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
first_bday_menb_completeness <- first_bday_menb %>% 
  group_by(C_patid) %>%
  summarise(menb = sum(menb))

# create expected dose counts
first_bday_menb_completeness <- first_bday_menb_completeness %>% 
  mutate(menb_exp = as.numeric(vaccine_exp_dose$first_bday_exp[vaccine_exp_dose$vaccine == "menb"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
first_bday_menb_completeness$menb_pc[first_bday_menb_completeness$menb_exp>first_bday_menb_completeness$menb] <- 0
first_bday_menb_completeness$menb_pc[first_bday_menb_completeness$menb_exp<=first_bday_menb_completeness$menb] <- 1

# check for NAs
sum(is.na(first_bday_menb_completeness$menb_pc))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
first_bday_menb_analysis <- first_bday_menb_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
first_bday_menb_analysis <- first_bday_menb_analysis %>% 
  mutate(first_bday = dob %m+% years(1)) %>%
  filter(first_bday >= "2006-04-01") %>% 
  filter(first_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (first_bday >= "2006-04-01" & first_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (first_bday >= "2007-04-01" & first_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (first_bday >= "2008-04-01" & first_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (first_bday >= "2009-04-01" & first_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (first_bday >= "2010-04-01" & first_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (first_bday >= "2011-04-01" & first_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (first_bday >= "2012-04-01" & first_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (first_bday >= "2013-04-01" & first_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (first_bday >= "2014-04-01" & first_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (first_bday >= "2015-04-01" & first_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (first_bday >= "2016-04-01" & first_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (first_bday >= "2017-04-01" & first_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (first_bday >= "2018-04-01" & first_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (first_bday >= "2019-04-01" & first_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (first_bday >= "2020-04-01" & first_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(first_bday_menb_analysis, file="~/first_bday_menb_analysis.Rdata")

# Create completeness & timing variables for second bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
second_bday_menb_completeness <- second_bday_menb %>% 
  group_by(C_patid) %>%
  summarise(menb = sum(menb))

# create expected dose counts
second_bday_menb_completeness <- second_bday_menb_completeness %>% 
  mutate(menb_exp_pc = as.numeric(vaccine_exp_dose$first_bday_exp[vaccine_exp_dose$vaccine == "menb"])) %>% 
  mutate(menb_exp_booster = as.numeric(vaccine_exp_dose$second_bday_exp[vaccine_exp_dose$vaccine == "menb"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
second_bday_menb_completeness$menb_pc[second_bday_menb_completeness$menb_exp_pc>second_bday_menb_completeness$menb] <- 0
second_bday_menb_completeness$menb_pc[second_bday_menb_completeness$menb_exp_pc<=second_bday_menb_completeness$menb] <- 1

# derive booster variable (complete = 1 vs incomplete = 0)
second_bday_menb_completeness$menb_booster[second_bday_menb_completeness$menb_exp_booster>second_bday_menb_completeness$menb] <- 0
second_bday_menb_completeness$menb_booster[second_bday_menb_completeness$menb_exp_booster<=second_bday_menb_completeness$menb] <- 1

# check for NAs
sum(is.na(second_bday_menb_completeness$menb_pc))
sum(is.na(second_bday_menb_completeness$menb_booster))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
second_bday_menb_analysis <- second_bday_menb_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
second_bday_menb_analysis <- second_bday_menb_analysis %>% 
  mutate(second_bday = dob %m+% years(2)) %>%
  filter(second_bday >= "2006-04-01") %>% 
  filter(second_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (second_bday >= "2006-04-01" & second_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (second_bday >= "2007-04-01" & second_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (second_bday >= "2008-04-01" & second_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (second_bday >= "2009-04-01" & second_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (second_bday >= "2010-04-01" & second_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (second_bday >= "2011-04-01" & second_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (second_bday >= "2012-04-01" & second_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (second_bday >= "2013-04-01" & second_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (second_bday >= "2014-04-01" & second_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (second_bday >= "2015-04-01" & second_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (second_bday >= "2016-04-01" & second_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (second_bday >= "2017-04-01" & second_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (second_bday >= "2018-04-01" & second_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (second_bday >= "2019-04-01" & second_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (second_bday >= "2020-04-01" & second_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(second_bday_menb_analysis, file="~/second_bday_menb_analysis.Rdata")

#--------------------------------------------------------------------------
# CLEAN PNEUMOCOCCAL (8w, 16w & 1y, then 12w & 1y since January 2020)                    
#--------------------------------------------------------------------------

# Identify all obs with medcodes and prodcodes  ---------------------------------------------

# load child observations and drug issue file
load(file="~/obs_C_clean.Rdata")
load(file="~/drug_C_clean.Rdata")
load(file="~/patids.Rdata")

# medcode file 
pneumococcal_medcodes <- read_tsv('~/pneumococcal_medcodes.txt', 
                          col_names = TRUE) %>% 
  mutate(medcodeid = as.character(medcodeid)) %>% 
  select(medcodeid)

# prodcode file
pneumococcal_prodcodes <- read_tsv('~/pneumococcal_prodcodes.txt', 
                           col_names = TRUE) %>% 
  mutate(prodcodeid = as.character(prodcodeid)) %>% 
  select(prodcodeid)

# identify all medcode obs
obs_pneumococcal_medcodes <- obs_C_clean %>% 
  inner_join(pneumococcal_medcodes, by = "medcodeid") %>% 
  mutate(pneumococcal = 1) %>% 
  mutate(code = "medcode") %>% 
  select(C_patid, obsdate, pneumococcal, code, consid)

# identify all prodcode obs
obs_pneumococcal_prodcodes <- drug_C_clean %>% 
  inner_join(pneumococcal_prodcodes, by = "prodcodeid") %>% 
  mutate(pneumococcal = 1) %>% 
  rename(obsdate = issuedate) %>% 
  mutate(code = "prodcode") %>% 
  mutate(consid = NA) %>% 
  select(C_patid, obsdate, pneumococcal, code, consid)

# join medcodes and prodcodes
obs_pneumococcal <- rbind(obs_pneumococcal_medcodes, obs_pneumococcal_prodcodes) 

# save for quality check analysis
save(obs_pneumococcal, file="~/obs_pneumococcal.Rdata")

# Counts & exclusions -----------------------------------------------------

# load C_patid list 
load(file = "~/sourcepop.Rdata")

babypatids <- sourcepop %>% 
  select(C_patid, dob, C_pc_fu_end) 

total_obs <- nrow(obs_pneumococcal)

# count missing obsdates
nodate_obs <- sum(is.na(obs_pneumococcal$obsdate)) 
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# count obsdates before DOB, after fifth birthday, and after follow up end
obs_pneumococcal <- obs_pneumococcal %>% 
  left_join(babypatids, by = "C_patid") %>% 
  mutate(fifth_bday = dob %m+% years(5)) 

obs_before_dob <- sum(obs_pneumococcal$obsdate < obs_pneumococcal$dob, na.rm = T)
percent_obs_before_dob <- round(obs_before_dob/total_obs*100, 3)

obs_after5y <- sum(obs_pneumococcal$obsdate > obs_pneumococcal$fifth_bday, na.rm = T)
percent_obs_after5y <- round(obs_after5y/total_obs*100, 3)

obs_after_fuend <- sum(obs_pneumococcal$obsdate > obs_pneumococcal$C_pc_fu_end, na.rm = T)
percent_obs_after_fuend<- round(obs_after_fuend/total_obs*100, 3)

# drop obs with no obsdate, obsdate before DOB, after fifth birthday and after follow up end
obs_pneumococcal_full <- obs_pneumococcal %>% 
  filter(!is.na(obsdate)) %>% 
  filter(obsdate >= dob) %>% 
  filter(obsdate <= fifth_bday) %>% 
  filter(obsdate <= C_pc_fu_end) 

# count number of med+prodcode obs after exclusions
obs_pneumococcal_excl <- nrow(obs_pneumococcal_full) 

# save 
save(obs_pneumococcal_full, file="~/obs_pneumococcal_full.Rdata")

# Cut down to 1 record per person per obsdate  ------------------------------------

# prioritise ones that have consid 
obs_pneumococcal_full <- obs_pneumococcal_full %>%
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, is.na(consid), 2)) %>%  
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  group_by(C_patid, obsdate) %>% 
  summarise(pneumococcal = first(pneumococcal),
            consid = first(consid))
n_distinct(obs_pneumococcal_full$C_patid) != nrow(obs_pneumococcal_full) # should not yet be one row per patient

# Join to first bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/first_bday_patids.Rdata")

# join first bday 
first_bday_pneumococcal <- first_bday_patids %>% 
  left_join(obs_pneumococcal_full, by = "C_patid") %>% 
  rename(obsdate_pneumococcal = obsdate) %>% 
  mutate(pneumococcal = replace_na(pneumococcal, 0))

# count obs occurring after first birthday
obs_first_bday_pneumococcal <- nrow(first_bday_pneumococcal) # denom is number of pneumococcal obs in this bday cohort 
pneumococcal_obs_after1y <- sum(first_bday_pneumococcal$obsdate_pneumococcal > first_bday_pneumococcal$first_bday, na.rm = T)
percent_pneumococcal_obs_after1y <- round(pneumococcal_obs_after1y/obs_first_bday_pneumococcal*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
first_bday_pneumococcal <- first_bday_pneumococcal %>%
  mutate(pneumococcal = replace(pneumococcal, obsdate_pneumococcal > first_bday, 0)) %>% 
  mutate(obsdate_pneumococcal = replace(obsdate_pneumococcal, obsdate_pneumococcal > first_bday, NA)) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
first_bday_pneumococcal <- first_bday_pneumococcal %>% 
  distinct(C_patid, pneumococcal, obsdate_pneumococcal, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(first_bday_pneumococcal$C_patid) == n_distinct(first_bday_patids$C_patid)

# save
save(first_bday_pneumococcal, file="~/first_bday_pneumococcal.Rdata")

# Join to second bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/second_bday_patids.Rdata")

# join second bday 
second_bday_pneumococcal <- second_bday_patids %>% 
  left_join(obs_pneumococcal_full, by = "C_patid") %>% 
  rename(obsdate_pneumococcal = obsdate) %>% 
  mutate(pneumococcal = replace_na(pneumococcal, 0))

# count obs occurring after second birthday
obs_second_bday_pneumococcal <- nrow(second_bday_pneumococcal) # denom is number of pneumococcal obs in this bday cohort 
pneumococcal_obs_after2y <- sum(second_bday_pneumococcal$obsdate_pneumococcal > second_bday_pneumococcal$second_bday, na.rm = T)
percent_pneumococcal_obs_after2y <- round(pneumococcal_obs_after2y/obs_second_bday_pneumococcal*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
second_bday_pneumococcal <- second_bday_pneumococcal %>%
  mutate(pneumococcal = replace(pneumococcal, obsdate_pneumococcal > second_bday, 0)) %>% 
  mutate(obsdate_pneumococcal = replace(obsdate_pneumococcal, obsdate_pneumococcal > second_bday, NA)) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
second_bday_pneumococcal <- second_bday_pneumococcal %>% 
  distinct(C_patid, pneumococcal, obsdate_pneumococcal, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(second_bday_pneumococcal$C_patid) == n_distinct(second_bday_patids$C_patid)

# save
save(second_bday_pneumococcal, file="~/second_bday_pneumococcal.Rdata")

# Save counts and checks -------------------------------------------------------------

total_obs <- as.data.frame(total_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
obs_before_dob <- as.data.frame(obs_before_dob)
percent_obs_before_dob <- as.data.frame(percent_obs_before_dob)
obs_after5y <- as.data.frame(obs_after5y)
percent_obs_after5y <- as.data.frame(percent_obs_after5y)
obs_after_fuend <- as.data.frame(obs_after_fuend)
percent_obs_after_fuend <- as.data.frame(percent_obs_after_fuend)
obs_pneumococcal_excl <- as.data.frame(obs_pneumococcal_excl)

obs_first_bday_pneumococcal <- as.data.frame(obs_first_bday_pneumococcal)
pneumococcal_obs_after1y <- as.data.frame(pneumococcal_obs_after1y)
percent_pneumococcal_obs_after1y <- as.data.frame(percent_pneumococcal_obs_after1y)
obs_second_bday_pneumococcal <- as.data.frame(obs_second_bday_pneumococcal)
pneumococcal_obs_after2y <- as.data.frame(pneumococcal_obs_after2y)
percent_pneumococcal_obs_after2y <- as.data.frame(percent_pneumococcal_obs_after2y)

pneumococcal_cleaning_counts <- cbind(total_obs, 
                              nodate_obs, percent_nodate_obs,
                              obs_before_dob, percent_obs_before_dob,
                              obs_after5y, percent_obs_after5y,
                              obs_after_fuend, percent_obs_after_fuend,
                              obs_pneumococcal_excl,
                              obs_first_bday_pneumococcal, pneumococcal_obs_after1y, percent_pneumococcal_obs_after1y,
                              obs_second_bday_pneumococcal, pneumococcal_obs_after2y, percent_pneumococcal_obs_after2y)

write_csv(pneumococcal_cleaning_counts, file = "~/pneumococcal_cleaning_counts.csv")


# Create completeness & timing variables for first bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# divide cohort by when the schedule changed (January 2020)
first_bday_pneumococcal <- first_bday_pneumococcal %>% 
  mutate(schedule_change = "before") %>% 
  mutate(schedule_change = replace(schedule_change, dob >= "2020-01-01", "after"))

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
first_bday_pneumococcal_completeness <- first_bday_pneumococcal %>% 
  group_by(C_patid) %>%
  summarise(pneumococcal = sum(pneumococcal),
            schedule_change = first(schedule_change))

# create expected dose counts
first_bday_pneumococcal_completeness <- first_bday_pneumococcal_completeness %>% 
  mutate(pneumococcal_exp = as.numeric(vaccine_exp_dose$first_bday_exp[vaccine_exp_dose$vaccine == "pneumococcal_before_2020"])) %>% 
  mutate(pneumococcal_exp = replace(pneumococcal_exp, schedule_change == "after", as.numeric(vaccine_exp_dose$first_bday_exp[vaccine_exp_dose$vaccine == "pneumococcal_after_2020"])))

# derive primary course variable (complete = 1 vs incomplete = 0)
first_bday_pneumococcal_completeness$pneumococcal_pc[first_bday_pneumococcal_completeness$pneumococcal_exp>first_bday_pneumococcal_completeness$pneumococcal] <- 0
first_bday_pneumococcal_completeness$pneumococcal_pc[first_bday_pneumococcal_completeness$pneumococcal_exp<=first_bday_pneumococcal_completeness$pneumococcal] <- 1

# check for NAs
sum(is.na(first_bday_pneumococcal_completeness$pneumococcal_pc))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
first_bday_pneumococcal_analysis <- second_bday_menb_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
first_bday_pneumococcal_analysis <- first_bday_pneumococcal_analysis %>% 
  mutate(first_bday = dob %m+% years(1)) %>%
  filter(first_bday >= "2006-04-01") %>% 
  filter(first_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (first_bday >= "2006-04-01" & first_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (first_bday >= "2007-04-01" & first_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (first_bday >= "2008-04-01" & first_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (first_bday >= "2009-04-01" & first_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (first_bday >= "2010-04-01" & first_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (first_bday >= "2011-04-01" & first_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (first_bday >= "2012-04-01" & first_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (first_bday >= "2013-04-01" & first_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (first_bday >= "2014-04-01" & first_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (first_bday >= "2015-04-01" & first_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (first_bday >= "2016-04-01" & first_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (first_bday >= "2017-04-01" & first_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (first_bday >= "2018-04-01" & first_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (first_bday >= "2019-04-01" & first_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (first_bday >= "2020-04-01" & first_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(first_bday_pneumococcal_analysis, file="~/first_bday_pneumococcal_analysis.Rdata")

# Create completeness & timing variables for second bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# no need to divide cohort by when the schedule changed since no one born after January 2020 would have reached their second birthday by end of the study period
# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
second_bday_pneumococcal_completeness <- second_bday_pneumococcal %>% 
  group_by(C_patid) %>%
  summarise(pneumococcal = sum(pneumococcal))

# create expected dose counts
second_bday_pneumococcal_completeness <- second_bday_pneumococcal_completeness %>% 
  mutate(pneumococcal_exp_pc = as.numeric(vaccine_exp_dose$first_bday_exp[vaccine_exp_dose$vaccine == "pneumococcal_before_2020"])) %>% 
  mutate(pneumococcal_exp_booster = as.numeric(vaccine_exp_dose$second_bday_exp[vaccine_exp_dose$vaccine == "pneumococcal_before_2020"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
second_bday_pneumococcal_completeness$pneumococcal_pc[second_bday_pneumococcal_completeness$pneumococcal_exp_pc>second_bday_pneumococcal_completeness$pneumococcal] <- 0
second_bday_pneumococcal_completeness$pneumococcal_pc[second_bday_pneumococcal_completeness$pneumococcal_exp_pc<=second_bday_pneumococcal_completeness$pneumococcal] <- 1

# derive booster variable (complete = 1 vs incomplete = 0)
second_bday_pneumococcal_completeness$pneumococcal_booster[second_bday_pneumococcal_completeness$pneumococcal_exp_booster>second_bday_pneumococcal_completeness$pneumococcal] <- 0
second_bday_pneumococcal_completeness$pneumococcal_booster[second_bday_pneumococcal_completeness$pneumococcal_exp_booster<=second_bday_pneumococcal_completeness$pneumococcal] <- 1

# check for NAs
sum(is.na(second_bday_pneumococcal_completeness$pneumococcal_pc))
sum(is.na(second_bday_pneumococcal_completeness$pneumococcal_booster))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
second_bday_pneumococcal_analysis <- second_bday_pneumococcal_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
second_bday_pneumococcal_analysis <- second_bday_pneumococcal_analysis %>% 
  mutate(second_bday = dob %m+% years(2)) %>%
  filter(second_bday >= "2006-04-01") %>% 
  filter(second_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (second_bday >= "2006-04-01" & second_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (second_bday >= "2007-04-01" & second_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (second_bday >= "2008-04-01" & second_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (second_bday >= "2009-04-01" & second_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (second_bday >= "2010-04-01" & second_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (second_bday >= "2011-04-01" & second_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (second_bday >= "2012-04-01" & second_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (second_bday >= "2013-04-01" & second_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (second_bday >= "2014-04-01" & second_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (second_bday >= "2015-04-01" & second_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (second_bday >= "2016-04-01" & second_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (second_bday >= "2017-04-01" & second_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (second_bday >= "2018-04-01" & second_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (second_bday >= "2019-04-01" & second_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (second_bday >= "2020-04-01" & second_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(second_bday_pneumococcal_analysis, file="~/second_bday_pneumococcal_analysis.Rdata")

#--------------------------------------------------------------------------
# PRELIM CLEANING FOR COMBINED INFANT VAX                    
#--------------------------------------------------------------------------

# Identify all obs with medcodes and prodcodes  ---------------------------------------------

# load child observations and drug issue file
load(file="~/obs_C_clean.Rdata")
load(file="~/drug_C_clean.Rdata")
load(file="~/patids.Rdata")

# medcode file 
combined_medcodes <- read_tsv('~/combined_vaccination_medcodes.txt', 
                                  col_names = TRUE) %>% 
  mutate(medcodeid = as.character(medcodeid)) %>% 
  select(-c(term, snomedct_conceptid, snomedct_descriptionid)) 

# prodcode file
combined_prodcodes <- read_tsv('~/combined_vaccination_prodcodes.txt', 
                                   col_names = TRUE) %>% 
  mutate(prodcodeid = as.character(prodcodeid)) %>% 
  select(-c(term, dmdid)) 

# identify all medcode obs
obs_combined_medcodes <- obs_C_clean %>% 
  inner_join(combined_medcodes, by = "medcodeid") %>% 
  mutate(code = "medcode") %>% 
  mutate(combined = 1) %>% 
  select(C_patid, obsdate, combined, 
         diptheria, tetanus, pertussis, polio, hepb, hib, menc, 
         comment, code, consid)

# identify all prodcode obs
obs_combined_prodcodes <- drug_C_clean %>% 
  inner_join(combined_prodcodes, by = "prodcodeid") %>% 
  mutate(code = "prodcode") %>% 
  rename(obsdate = issuedate) %>% 
  mutate(consid = NA) %>% 
  mutate(combined = 1) %>% 
  select(C_patid, obsdate, combined, 
         diptheria, tetanus, pertussis, polio, hepb, hib, menc, 
         comment, code, consid)

# join medcodes and prodcodes
obs_combined <- rbind(obs_combined_medcodes, obs_combined_prodcodes) %>%  
  mutate(id = row_number()) 
  
# load C_patid list 
load(file = "~/sourcepop.Rdata")

# create first and third birthday variables (needed in the algorithm)
babypatids <- sourcepop %>% 
  select(C_patid, dob, C_pc_fu_end) 

obs_combined <- obs_combined %>% 
  left_join(babypatids, by = "C_patid") %>% 
  mutate(first_bday = dob %m+% years(1)) %>% 
  mutate(third_bday = dob %m+% years(3)) 

combined_total <- nrow(obs_combined)

# save
save(obs_combined, file = "~/obs_combined.Rdata")

# Assign vaccines with clear codes using algorithm  ---------------------------

## clear indication of 6-in-1 
obs_6in1_clear <- obs_combined %>% 
  filter(diptheria == 1 & 
           tetanus == 1 &
           pertussis == 1 &
           polio == 1 &
           hepb == 1 &
           hib == 1 &
           is.na(menc))
obs_6in1_clear$vaccine <- "6-in-1"
obs_6in1_clear$type <- "full"

## clear indication of 5-in-1 
obs_5in1_clear <- obs_combined %>% 
  filter(diptheria == 1 & 
           tetanus == 1 &
           pertussis == 1 &
           polio == 1 &
           is.na(hepb) &
           hib == 1 &
           is.na(menc))
obs_5in1_clear$vaccine <- "5-in-1"
obs_5in1_clear$type <- "full"

## clear indication of 4-in-1 
obs_4in1_clear <- obs_combined %>% 
  filter(obsdate >= third_bday) %>% 
  filter(diptheria == 1 & 
           tetanus == 1 &
           pertussis == 1 &
           polio == 1 &
           is.na(hib) &
           is.na(hepb) &
           is.na(menc))
obs_4in1_clear$vaccine <- "4-in-1"
obs_4in1_clear$type <- "full"

## clear indication of hib/menc 
obs_hibmenc_clear <- obs_combined %>% 
  filter(is.na(diptheria) & 
           is.na(tetanus) &
           is.na(pertussis) &
           is.na(polio) &
           is.na(hepb) &
           hib == 1 &
           menc == 1)
obs_hibmenc_clear$vaccine <- "hibmenc"
obs_hibmenc_clear$type <- "full"

# Assign vaccines with unclear codes using algorithm ---------------------------

# remove 'clear' codes from master files before finding 'unclear' codes
obs_clear <- rbind(obs_6in1_clear, obs_5in1_clear, obs_4in1_clear, obs_hibmenc_clear)

obs_unclear <- obs_combined %>% 
  anti_join(obs_clear, by = "id")

nrow(obs_unclear) + nrow(obs_clear) == nrow(obs_combined)

## 6-in-1
obs_6in1_unclear <- obs_unclear %>% 
  filter(obsdate < third_bday) %>% 
  filter(dob >= "2017-08-01") %>% 
  filter(diptheria == 1 | 
           tetanus == 1 |
           pertussis == 1 |
           polio == 1 |
           (hepb == 1 & obsdate >= "2017-08-01"))
obs_6in1_unclear$vaccine <- "6-in-1"
obs_6in1_unclear$type <- "partial"

# remove the above codes from master file to avoid duplicate filtering at the next step
obs_unclear <- obs_unclear %>% 
  anti_join(obs_6in1_unclear, by = "id")

obs_6in1_hib_unclear <- obs_unclear %>% 
  filter(obsdate < first_bday) %>% 
  filter(dob >= "2017-08-01") %>% 
  filter(hib == 1) 
obs_6in1_hib_unclear$vaccine <- "6-in-1"
obs_6in1_hib_unclear$type <- "partial"

obs_unclear <- obs_unclear %>% 
  anti_join(obs_6in1_hib_unclear, by = "id")

## 5-in-1
obs_5in1_unclear <- obs_unclear %>% 
  filter(obsdate < third_bday) %>% 
  filter(dob < "2017-08-01") %>% 
  filter(diptheria == 1 | 
           tetanus == 1 |
           pertussis == 1 |
           polio == 1)
obs_5in1_unclear$vaccine <- "5-in-1"
obs_5in1_unclear$type <- "partial"

obs_unclear <- obs_unclear %>% 
  anti_join(obs_5in1_unclear, by = "id")

obs_5in1_hib_unclear <- obs_unclear %>% 
  filter(obsdate < first_bday) %>% 
  filter(obsdate < "2017-08-01") %>% 
  filter(hib == 1)
obs_5in1_hib_unclear$vaccine <- "5-in-1"
obs_5in1_hib_unclear$type <- "partial"

obs_unclear <- obs_unclear %>% 
  anti_join(obs_5in1_hib_unclear, by = "id")

## 4-in-1
obs_4in1_unclear <- obs_unclear %>% 
  filter(obsdate >= third_bday) %>% 
  filter(diptheria == 1 | 
           tetanus == 1 |
           pertussis == 1 |
           polio == 1)
obs_4in1_unclear$vaccine <- "4-in-1"
obs_4in1_unclear$type <- "partial"

obs_unclear <- obs_unclear %>% 
  anti_join(obs_4in1_unclear, by = "id")

## hib/menc
obs_hibmenc_unclear <- obs_unclear %>% 
  filter(obsdate >= first_bday) %>% 
  filter(hib == 1)
obs_hibmenc_unclear$vaccine <- "hibmenc"
obs_hibmenc_unclear$type <- "partial"

obs_hibmenc_unclear2 <- obs_unclear %>% 
  filter(dob >= "2016-07-01") %>% 
  filter(menc == 1)
obs_hibmenc_unclear2$vaccine <- "hibmenc"
obs_hibmenc_unclear2$type <- "partial"

obs_hibmenc_unclear3 <- obs_unclear %>% 
  filter(dob < "2016-07-01") %>% 
  filter(obsdate >= first_bday) %>% 
  filter(menc == 1)
obs_hibmenc_unclear3$vaccine <- "hibmenc"
obs_hibmenc_unclear3$type <- "partial"

## menc only
obs_menc_unclear <- obs_unclear %>% 
  filter(dob < "2016-07-01") %>% 
  filter(obsdate < first_bday) %>% 
  filter(menc == 1)
obs_menc_unclear$vaccine <- "menc"
obs_menc_unclear$type <- "full"

# Join clear and unclear codes for each vaccine ---------------------------

# join medcodes and check no duplicate obs ids
obs_6in1 <- rbind(obs_6in1_clear, obs_6in1_unclear, obs_6in1_hib_unclear)
obs_6in1 <- unique(obs_6in1)

obs_5in1 <- rbind(obs_5in1_clear, obs_5in1_unclear, obs_5in1_hib_unclear)
obs_5in1 <- unique(obs_5in1)

obs_4in1 <- rbind(obs_4in1_clear, obs_4in1_unclear)
obs_4in1 <- unique(obs_4in1)

obs_hibmenc <- rbind(obs_hibmenc_clear, obs_hibmenc_unclear, obs_hibmenc_unclear2, obs_hibmenc_unclear3)
obs_hibmenc <- unique(obs_hibmenc)

obs_menc <- obs_menc_unclear

# check codes that did not get picked up by the algorithm
obs_algorithm <- rbind(obs_6in1, obs_5in1, obs_4in1, obs_hibmenc, obs_menc)
obs_remaining <- obs_combined %>% 
  anti_join(obs_algorithm, by = "id") # all single dose hepb or missing obsdate so couldn't be processed using algorithm
combined_processed <- nrow(obs_algorithm)
combined_notprocessed <- nrow(obs_remaining)

# save for quality check analysis
save(obs_6in1, file = "~/obs_6in1.Rdata")
save(obs_5in1, file = "~/obs_5in1.Rdata")
save(obs_4in1, file = "~/obs_4in1.Rdata")
save(obs_hibmenc, file = "~/obs_hibmenc.Rdata")
save(obs_menc, file = "~/obs_menc.Rdata")

# save counts
combined_total <- as.data.frame(combined_total)
combined_processed <- as.data.frame(combined_processed)
combined_notprocessed <- as.data.frame(combined_notprocessed)

combined_cleaning_counts <- cbind(combined_total, combined_processed, combined_notprocessed)

write_csv(combined_cleaning_counts, file = "~/combined_cleaning_counts.csv")

#--------------------------------------------------------------------------
# 6-in-1 (8, 12 & 16w since August 2017)                   
#--------------------------------------------------------------------------

# Counts & exclusions -----------------------------------------------------

load(file = "~/obs_6in1.Rdata")

total_obs <- nrow(obs_6in1)

# count missing obsdates
nodate_obs <- sum(is.na(obs_6in1$obsdate)) 
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# count obsdates before DOB, before August 2017, after fifth birthday, and after follow up end
obs_6in1 <- obs_6in1 %>% 
  mutate(fifth_bday = dob %m+% years(5)) 

obs_before_dob <- sum(obs_6in1$obsdate < obs_6in1$dob, na.rm = T)
percent_obs_before_dob <- round(obs_before_dob/total_obs*100, 3)

obs_before_aug2017<- sum(obs_6in1$obsdate < "2017-08-01", na.rm = T)
percent_obs_before_aug2017 <- round(obs_before_aug2017/total_obs*100, 3)

obs_after5y <- sum(obs_6in1$obsdate > obs_6in1$fifth_bday, na.rm = T)
percent_obs_after5y <- round(obs_after5y/total_obs*100, 3)

obs_after_fuend <- sum(obs_6in1$obsdate > obs_6in1$C_pc_fu_end, na.rm = T)
percent_obs_after_fuend<- round(obs_after_fuend/total_obs*100, 3)

# drop obs with no obsdate, obsdate before DOB, before August 2017, after fifth birthday and after follow up end
obs_6in1_full <- obs_6in1 %>% 
  filter(!is.na(obsdate)) %>% 
  filter(obsdate >= dob) %>% 
  filter(obsdate >= "2017-08-01") %>% 
  filter(obsdate <= fifth_bday) %>% 
  filter(obsdate <= C_pc_fu_end) 

# count number of med+prodcode obs after exclusions
obs_6in1_excl <- nrow(obs_6in1_full) 

# save 
save(obs_6in1_full, file="~/obs_6in1_full.Rdata")

# Cut down to 1 record per person per obsdate  ------------------------------------

# prioritise full over partial, then medcode over prodcode
# if multiple of the same priority, prioritise ones that have consid 
obs_6in1_full <- obs_6in1_full %>%
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, type == "full" & code == "prodcode", 2)) %>% 
  mutate(priority = replace(priority, type == "partial" & code == "medcode", 3)) %>% 
  mutate(priority = replace(priority, type == "partial" & code == "prodcode", 4)) %>% 
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  ungroup() %>% 
  select(-priority) %>% 
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, is.na(consid), 2)) %>%  
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  group_by(C_patid, obsdate) %>% 
  summarise(`6in1` = 1,
            consid = first(consid))
n_distinct(obs_6in1_full$C_patid) != nrow(obs_6in1_full) # should not yet be one row per patient

# Join to first bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/first_bday_patids.Rdata")

# join first bday for those born on or after August 2017
first_bday_patids_aug2017 <- first_bday_patids %>%
  filter(dob >= "2017-08-01")

first_bday_6in1 <- first_bday_patids_aug2017 %>% 
  left_join(obs_6in1_full, by = "C_patid") %>% 
  rename(obsdate_6in1 = obsdate) %>% 
  mutate(`6in1` = replace_na(`6in1`, 0))

# count obs occurring after first birthday
obs_first_bday_6in1 <- nrow(first_bday_6in1) # denom is number of `6in1` obs in this bday cohort 
obs_6in1_after1y <- sum(first_bday_6in1$obsdate_6in1 > first_bday_6in1$first_bday, na.rm = T)
percent_obs_6in1_after1y <- round(obs_6in1_after1y/obs_first_bday_6in1*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
first_bday_6in1 <- first_bday_6in1 %>%
  mutate(`6in1` = replace(`6in1`, obsdate_6in1 > first_bday, 0)) %>% 
  mutate(obsdate_6in1 = replace(obsdate_6in1, obsdate_6in1 > first_bday, NA)) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
first_bday_6in1 <- first_bday_6in1 %>% 
  distinct(C_patid, `6in1`, obsdate_6in1, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(first_bday_6in1$C_patid) == n_distinct(first_bday_patids$C_patid[first_bday_patids$dob >= "2017-08-01"])

# save
save(first_bday_6in1, file="~/first_bday_6in1.Rdata")

# Join to second bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/second_bday_patids.Rdata")

# join second bday for those born on or after August 2017
second_bday_patids_aug2017 <- second_bday_patids %>%
  filter(dob >= "2017-08-01")

second_bday_6in1 <- second_bday_patids_aug2017 %>% 
  left_join(obs_6in1_full, by = "C_patid") %>% 
  rename(obsdate_6in1 = obsdate) %>% 
  mutate(`6in1` = replace_na(`6in1`, 0))

# count obs occurring after second birthday
obs_second_bday_6in1 <- nrow(second_bday_6in1) # denom is number of `6in1` obs in this bday cohort 
obs_6in1_after2y <- sum(second_bday_6in1$obsdate_6in1 > second_bday_6in1$second_bday, na.rm = T)
percent_obs_6in1_after2y <- round(obs_6in1_after2y/obs_second_bday_6in1*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
second_bday_6in1 <- second_bday_6in1 %>%
  mutate(`6in1` = replace(`6in1`, obsdate_6in1 > second_bday, 0)) %>% 
  mutate(obsdate_6in1 = replace(obsdate_6in1, obsdate_6in1 > second_bday, NA)) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
second_bday_6in1 <- second_bday_6in1 %>% 
  distinct(C_patid, `6in1`, obsdate_6in1, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(second_bday_6in1$C_patid) == n_distinct(second_bday_patids$C_patid[second_bday_patids$dob >= "2017-08-01"])

# save
save(second_bday_6in1, file="~/second_bday_6in1.Rdata")

# Save counts and checks -------------------------------------------------------------

total_obs <- as.data.frame(total_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
obs_before_dob <- as.data.frame(obs_before_dob)
percent_obs_before_dob <- as.data.frame(percent_obs_before_dob)
obs_before_aug2017 <- as.data.frame(obs_before_aug2017)
percent_obs_before_aug2017 <- as.data.frame(percent_obs_before_aug2017)
obs_after5y <- as.data.frame(obs_after5y)
percent_obs_after5y <- as.data.frame(percent_obs_after5y)
obs_after_fuend <- as.data.frame(obs_after_fuend)
percent_obs_after_fuend <- as.data.frame(percent_obs_after_fuend)
obs_6in1_excl <- as.data.frame(obs_6in1_excl)

obs_first_bday_6in1 <- as.data.frame(obs_first_bday_6in1)
obs_6in1_after1y <- as.data.frame(obs_6in1_after1y)
percent_obs_6in1_after1y <- as.data.frame(percent_obs_6in1_after1y)
obs_second_bday_6in1 <- as.data.frame(obs_second_bday_6in1)
obs_6in1_after2y <- as.data.frame(obs_6in1_after2y)
percent_obs_6in1_after2y <- as.data.frame(percent_obs_6in1_after2y)

cleaning_6in1_counts <- cbind(total_obs, 
                              nodate_obs, percent_nodate_obs,
                              obs_before_dob, percent_obs_before_dob,
                              obs_before_aug2017, percent_obs_before_aug2017,
                              obs_after5y, percent_obs_after5y,
                              obs_after_fuend, percent_obs_after_fuend,
                              obs_6in1_excl,
                              obs_first_bday_6in1, obs_6in1_after1y, percent_obs_6in1_after1y,
                              obs_second_bday_6in1, obs_6in1_after2y, percent_obs_6in1_after2y)

write_csv(cleaning_6in1_counts, file = "~/cleaning_6in1_counts.csv")


# Create completeness & timing variables for first bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
first_bday_6in1_completeness <- first_bday_6in1 %>% 
  group_by(C_patid) %>%
  summarise(`6in1` = sum(`6in1`))

# create expected dose counts
first_bday_6in1_completeness <- first_bday_6in1_completeness %>% 
  mutate(`6in1_exp_pc` = as.numeric(vaccine_exp_dose$first_bday_exp[vaccine_exp_dose$vaccine == "6_5in1"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
first_bday_6in1_completeness$`6in1_pc`[first_bday_6in1_completeness$`6in1_exp_pc`>first_bday_6in1_completeness$`6in1`] <- 0
first_bday_6in1_completeness$`6in1_pc`[first_bday_6in1_completeness$`6in1_exp_pc`<=first_bday_6in1_completeness$`6in1`] <- 1

# check for NAs
sum(is.na(first_bday_6in1_completeness$`6in1_pc`))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
first_bday_6in1_analysis <- first_bday_6in1_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
first_bday_6in1_analysis <- first_bday_6in1_analysis %>% 
  mutate(first_bday = dob %m+% years(1)) %>%
  filter(first_bday >= "2006-04-01") %>% 
  filter(first_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (first_bday >= "2006-04-01" & first_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (first_bday >= "2007-04-01" & first_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (first_bday >= "2008-04-01" & first_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (first_bday >= "2009-04-01" & first_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (first_bday >= "2010-04-01" & first_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (first_bday >= "2011-04-01" & first_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (first_bday >= "2012-04-01" & first_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (first_bday >= "2013-04-01" & first_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (first_bday >= "2014-04-01" & first_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (first_bday >= "2015-04-01" & first_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (first_bday >= "2016-04-01" & first_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (first_bday >= "2017-04-01" & first_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (first_bday >= "2018-04-01" & first_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (first_bday >= "2019-04-01" & first_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (first_bday >= "2020-04-01" & first_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(first_bday_6in1_analysis, file="~/first_bday_6in1_analysis.Rdata")

# Create completeness & timing variables for second bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
second_bday_6in1_completeness <- second_bday_6in1 %>% 
  group_by(C_patid) %>%
  summarise(`6in1` = sum(`6in1`))

# create expected dose counts
second_bday_6in1_completeness <- second_bday_6in1_completeness %>% 
  mutate(`6in1_exp` = as.numeric(vaccine_exp_dose$second_bday_exp[vaccine_exp_dose$vaccine == "6_5in1"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
second_bday_6in1_completeness$`6in1_pc`[second_bday_6in1_completeness$`6in1_exp`>second_bday_6in1_completeness$`6in1`] <- 0
second_bday_6in1_completeness$`6in1_pc`[second_bday_6in1_completeness$`6in1_exp`<=second_bday_6in1_completeness$`6in1`] <- 1

# check for NAs
sum(is.na(second_bday_6in1_completeness$`6in1_pc`))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
second_bday_6in1_analysis <- second_bday_6in1_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
second_bday_6in1_analysis <- second_bday_6in1_analysis %>% 
  mutate(second_bday = dob %m+% years(2)) %>%
  filter(second_bday >= "2006-04-01") %>% 
  filter(second_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (second_bday >= "2006-04-01" & second_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (second_bday >= "2007-04-01" & second_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (second_bday >= "2008-04-01" & second_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (second_bday >= "2009-04-01" & second_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (second_bday >= "2010-04-01" & second_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (second_bday >= "2011-04-01" & second_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (second_bday >= "2012-04-01" & second_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (second_bday >= "2013-04-01" & second_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (second_bday >= "2014-04-01" & second_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (second_bday >= "2015-04-01" & second_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (second_bday >= "2016-04-01" & second_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (second_bday >= "2017-04-01" & second_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (second_bday >= "2018-04-01" & second_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (second_bday >= "2019-04-01" & second_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (second_bday >= "2020-04-01" & second_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(second_bday_6in1_analysis, file="~/second_bday_6in1_analysis.Rdata")

#--------------------------------------------------------------------------
# 5-in-1 (8, 12 & 16w until July 2017)                   
#--------------------------------------------------------------------------

# Counts & exclusions -----------------------------------------------------

load(file = "~/obs_5in1.Rdata")

total_obs <- nrow(obs_5in1)

# count missing obsdates
nodate_obs <- sum(is.na(obs_5in1$obsdate)) # 87
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# count obsdates before DOB, after July 2017, after fifth birthday, and after follow up end
obs_5in1 <- obs_5in1 %>% 
  mutate(fifth_bday = dob %m+% years(5)) 

obs_before_dob <- sum(obs_5in1$obsdate < obs_5in1$dob, na.rm = T)
percent_obs_before_dob <- round(obs_before_dob/total_obs*100, 3)

obs_after_jul2017<- sum(obs_5in1$obsdate > "2017-07-31", na.rm = T)
percent_obs_after_jul2017 <- round(obs_after_jul2017/total_obs*100, 3)

obs_after5y <- sum(obs_5in1$obsdate > obs_5in1$fifth_bday, na.rm = T)
percent_obs_after5y <- round(obs_after5y/total_obs*100, 3)

obs_after_fuend <- sum(obs_5in1$obsdate > obs_5in1$C_pc_fu_end, na.rm = T)
percent_obs_after_fuend<- round(obs_after_fuend/total_obs*100, 3)

# drop obs with no obsdate, obsdate before DOB, after July 2017, after fifth birthday and after follow up end
obs_5in1_full <- obs_5in1 %>% 
  filter(!is.na(obsdate)) %>% 
  filter(obsdate >= dob) %>% 
  filter(obsdate <= "2017-07-31") %>% 
  filter(obsdate <= fifth_bday) %>% 
  filter(obsdate <= C_pc_fu_end) 

# count number of med+prodcode obs after exclusions
obs_5in1_excl <- nrow(obs_5in1_full) 

# save 
save(obs_5in1_full, file="~/obs_5in1_full.Rdata")

# Cut down to 1 record per person per obsdate  ------------------------------------

# prioritise full over partial, then medcode over prodcode
# if multiple of the same priority, prioritise ones that have consid 
obs_5in1_full <- obs_5in1_full %>%
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, type == "full" & code == "prodcode", 2)) %>% 
  mutate(priority = replace(priority, type == "partial" & code == "medcode", 3)) %>% 
  mutate(priority = replace(priority, type == "partial" & code == "prodcode", 4)) %>% 
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  ungroup() %>% 
  select(-priority) %>% 
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, is.na(consid), 2)) %>%  
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  group_by(C_patid, obsdate) %>% 
  summarise(`5in1` = 1,
            consid = first(consid))
n_distinct(obs_5in1_full$C_patid) != nrow(obs_5in1_full) # should not yet be one row per patient

# Join to first bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/first_bday_patids.Rdata")

# join first bday for those born before August 2017
first_bday_patids_befaug2017 <- first_bday_patids %>%
  filter(dob < "2017-08-01")

first_bday_5in1 <- first_bday_patids_befaug2017 %>% 
  left_join(obs_5in1_full, by = "C_patid") %>% 
  rename(obsdate_5in1 = obsdate) %>% 
  mutate(`5in1` = replace_na(`5in1`, 0))

# count obs occurring after first birthday
obs_first_bday_5in1 <- nrow(first_bday_5in1) # denom is number of `5in1` obs in this bday cohort 
obs_5in1_after1y <- sum(first_bday_5in1$obsdate_5in1 > first_bday_5in1$first_bday, na.rm = T)
percent_obs_5in1_after1y <- round(obs_5in1_after1y/obs_first_bday_5in1*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
first_bday_5in1 <- first_bday_5in1 %>%
  mutate(`5in1` = replace(`5in1`, obsdate_5in1 > first_bday, 0)) %>% 
  mutate(obsdate_5in1 = replace(obsdate_5in1, obsdate_5in1 > first_bday, NA)) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
first_bday_5in1 <- first_bday_5in1 %>% 
  distinct(C_patid, `5in1`, obsdate_5in1, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(first_bday_5in1$C_patid) == n_distinct(first_bday_patids$C_patid[first_bday_patids$dob < "2017-08-01"])

# save
save(first_bday_5in1, file="~/first_bday_5in1.Rdata")

# Join to second bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/second_bday_patids.Rdata")

# join second bday for those born before August 2017
second_bday_patids_befaug2017 <- second_bday_patids %>%
  filter(dob < "2017-08-01")

second_bday_5in1 <- second_bday_patids_befaug2017 %>% 
  left_join(obs_5in1_full, by = "C_patid") %>% 
  rename(obsdate_5in1 = obsdate) %>% 
  mutate(`5in1` = replace_na(`5in1`, 0))

# count obs occurring after second birthday
obs_second_bday_5in1 <- nrow(second_bday_5in1) # denom is number of `5in1` obs in this bday cohort 
obs_5in1_after2y <- sum(second_bday_5in1$obsdate_5in1 > second_bday_5in1$second_bday, na.rm = T)
percent_obs_5in1_after2y <- round(obs_5in1_after2y/obs_second_bday_5in1*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
second_bday_5in1 <- second_bday_5in1 %>%
  mutate(`5in1` = replace(`5in1`, obsdate_5in1 > second_bday, 0)) %>% 
  mutate(obsdate_5in1 = replace(obsdate_5in1, obsdate_5in1 > second_bday, NA)) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
second_bday_5in1 <- second_bday_5in1 %>% 
  distinct(C_patid, `5in1`, obsdate_5in1, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(second_bday_5in1$C_patid) == n_distinct(second_bday_patids$C_patid[second_bday_patids$dob < "2017-08-01"])

# save
save(second_bday_5in1, file="~/second_bday_5in1.Rdata")

# Join to fifth bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/fifth_bday_patids.Rdata")

# join fifth bday for those born before August 2017 
fifth_bday_patids_befaug2017 <- fifth_bday_patids %>% 
  filter(dob < "2017-08-01") 

fifth_bday_5in1 <- fifth_bday_patids_befaug2017 %>% 
  left_join(obs_5in1_full, by = "C_patid") %>% 
  rename(obsdate_5in1 = obsdate) %>% 
  mutate(`5in1` = replace_na(`5in1`, 0))

obs_fifth_bday_5in1 <- nrow(fifth_bday_5in1) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
fifth_bday_5in1 <- fifth_bday_5in1 %>% 
  distinct(C_patid, `5in1`, obsdate_5in1, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(fifth_bday_5in1$C_patid) == n_distinct(fifth_bday_patids$C_patid[fifth_bday_patids$dob < "2017-08-01"])

# save
save(fifth_bday_5in1, file="~/fifth_bday_5in1.Rdata")

# Save counts and checks -------------------------------------------------------------

total_obs <- as.data.frame(total_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
obs_before_dob <- as.data.frame(obs_before_dob)
percent_obs_before_dob <- as.data.frame(percent_obs_before_dob)
obs_after_jul2017 <- as.data.frame(obs_after_jul2017)
percent_obs_after_jul2017 <- as.data.frame(percent_obs_after_jul2017)
obs_after5y <- as.data.frame(obs_after5y)
percent_obs_after5y <- as.data.frame(percent_obs_after5y)
obs_after_fuend <- as.data.frame(obs_after_fuend)
percent_obs_after_fuend <- as.data.frame(percent_obs_after_fuend)
obs_5in1_excl <- as.data.frame(obs_5in1_excl)

obs_first_bday_5in1 <- as.data.frame(obs_first_bday_5in1)
obs_5in1_after1y <- as.data.frame(obs_5in1_after1y)
percent_obs_5in1_after1y <- as.data.frame(percent_obs_5in1_after1y)
obs_second_bday_5in1 <- as.data.frame(obs_second_bday_5in1)
obs_5in1_after2y <- as.data.frame(obs_5in1_after2y)
percent_obs_5in1_after2y <- as.data.frame(percent_obs_5in1_after2y)
obs_fifth_bday_5in1 <- as.data.frame(obs_fifth_bday_5in1)

cleaning_5in1_counts <- cbind(total_obs, 
                              nodate_obs, percent_nodate_obs,
                              obs_before_dob, percent_obs_before_dob,
                              obs_after_jul2017, percent_obs_after_jul2017,
                              obs_after5y, percent_obs_after5y,
                              obs_after_fuend, percent_obs_after_fuend,
                              obs_5in1_excl,
                              obs_first_bday_5in1, obs_5in1_after1y, percent_obs_5in1_after1y,
                              obs_second_bday_5in1, obs_5in1_after2y, percent_obs_5in1_after2y,
                              obs_fifth_bday_5in1)

write_csv(cleaning_5in1_counts, file = "~/cleaning_5in1_counts.csv")

# Create completeness & timing variables for first bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
first_bday_5in1_completeness <- first_bday_5in1 %>% 
  group_by(C_patid) %>%
  summarise(`5in1` = sum(`5in1`))

# create expected dose counts
first_bday_5in1_completeness <- first_bday_5in1_completeness %>% 
  mutate(`5in1_exp` = as.numeric(vaccine_exp_dose$first_bday_exp[vaccine_exp_dose$vaccine == "6_5in1"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
first_bday_5in1_completeness$`5in1_pc`[first_bday_5in1_completeness$`5in1_exp`>first_bday_5in1_completeness$`5in1`] <- 0
first_bday_5in1_completeness$`5in1_pc`[first_bday_5in1_completeness$`5in1_exp`<=first_bday_5in1_completeness$`5in1`] <- 1

# check for NAs
sum(is.na(first_bday_5in1_completeness$`5in1_pc`))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
first_bday_5in1_analysis <- first_bday_5in1_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
first_bday_5in1_analysis <- first_bday_5in1_analysis %>% 
  mutate(first_bday = dob %m+% years(1)) %>%
  filter(first_bday >= "2006-04-01") %>% 
  filter(first_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (first_bday >= "2006-04-01" & first_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (first_bday >= "2007-04-01" & first_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (first_bday >= "2008-04-01" & first_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (first_bday >= "2009-04-01" & first_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (first_bday >= "2010-04-01" & first_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (first_bday >= "2011-04-01" & first_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (first_bday >= "2012-04-01" & first_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (first_bday >= "2013-04-01" & first_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (first_bday >= "2014-04-01" & first_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (first_bday >= "2015-04-01" & first_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (first_bday >= "2016-04-01" & first_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (first_bday >= "2017-04-01" & first_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (first_bday >= "2018-04-01" & first_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (first_bday >= "2019-04-01" & first_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (first_bday >= "2020-04-01" & first_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(first_bday_5in1_analysis, file="~/first_bday_5in1_analysis.Rdata")

# Create completeness & timing variables for second bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
second_bday_5in1_completeness <- second_bday_5in1 %>% 
  group_by(C_patid) %>%
  summarise(`5in1` = sum(`5in1`))

# create expected dose counts
second_bday_5in1_completeness <- second_bday_5in1_completeness %>% 
  mutate(`5in1_exp` = as.numeric(vaccine_exp_dose$second_bday_exp[vaccine_exp_dose$vaccine == "6_5in1"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
second_bday_5in1_completeness$`5in1_pc`[second_bday_5in1_completeness$`5in1_exp`>second_bday_5in1_completeness$`5in1`] <- 0
second_bday_5in1_completeness$`5in1_pc`[second_bday_5in1_completeness$`5in1_exp`<=second_bday_5in1_completeness$`5in1`] <- 1

# check for NAs
sum(is.na(second_bday_5in1_completeness$`5in1_pc`))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
second_bday_5in1_analysis <- second_bday_5in1_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
second_bday_5in1_analysis <- second_bday_5in1_analysis %>% 
  mutate(second_bday = dob %m+% years(2)) %>%
  filter(second_bday >= "2006-04-01") %>% 
  filter(second_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (second_bday >= "2006-04-01" & second_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (second_bday >= "2007-04-01" & second_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (second_bday >= "2008-04-01" & second_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (second_bday >= "2009-04-01" & second_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (second_bday >= "2010-04-01" & second_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (second_bday >= "2011-04-01" & second_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (second_bday >= "2012-04-01" & second_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (second_bday >= "2013-04-01" & second_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (second_bday >= "2014-04-01" & second_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (second_bday >= "2015-04-01" & second_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (second_bday >= "2016-04-01" & second_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (second_bday >= "2017-04-01" & second_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (second_bday >= "2018-04-01" & second_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (second_bday >= "2019-04-01" & second_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (second_bday >= "2020-04-01" & second_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(second_bday_5in1_analysis, file="~/second_bday_5in1_analysis.Rdata")

# Create completeness & timing variables for fifth bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
fifth_bday_5in1_completeness <- fifth_bday_5in1 %>% 
  group_by(C_patid) %>%
  summarise(`5in1` = sum(`5in1`))

# create expected dose counts
fifth_bday_5in1_completeness <- fifth_bday_5in1_completeness %>% 
  mutate(`5in1_exp_pc` = as.numeric(vaccine_exp_dose$fifth_bday_exp[vaccine_exp_dose$vaccine == "6_5in1"])) 

# derive primary course variable (complete = 1 vs incomplete = 0)
fifth_bday_5in1_completeness$`5in1_pc`[fifth_bday_5in1_completeness$`5in1_exp_pc`>fifth_bday_5in1_completeness$`5in1`] <- 0
fifth_bday_5in1_completeness$`5in1_pc`[fifth_bday_5in1_completeness$`5in1_exp_pc`<=fifth_bday_5in1_completeness$`5in1`] <- 1

# check for NAs
sum(is.na(fifth_bday_5in1_completeness$`5in1_pc`))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
fifth_bday_5in1_analysis <- fifth_bday_5in1_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
fifth_bday_5in1_analysis <- fifth_bday_5in1_analysis %>% 
  mutate(fifth_bday = dob %m+% years(5)) %>%
  filter(fifth_bday >= "2006-04-01") %>% 
  filter(fifth_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (fifth_bday >= "2006-04-01" & fifth_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2007-04-01" & fifth_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2008-04-01" & fifth_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2009-04-01" & fifth_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2010-04-01" & fifth_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2011-04-01" & fifth_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2012-04-01" & fifth_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2013-04-01" & fifth_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2014-04-01" & fifth_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2015-04-01" & fifth_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2016-04-01" & fifth_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2017-04-01" & fifth_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2018-04-01" & fifth_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2019-04-01" & fifth_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2020-04-01" & fifth_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(fifth_bday_5in1_analysis, file="~/fifth_bday_5in1_analysis.Rdata")

#--------------------------------------------------------------------------
# 4-in-1 PRESCHOOL BOOSTER (3y4m)                   
#--------------------------------------------------------------------------

# Counts & exclusions -----------------------------------------------------

load(file = "~/obs_4in1.Rdata")

total_obs <- nrow(obs_4in1)

# count missing obsdates
nodate_obs <- sum(is.na(obs_4in1$obsdate)) 
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# count obsdates before DOB, before 3 years (since booster isn't due until 3y4m), after fifth birthday, and after follow up end
obs_4in1 <- obs_4in1 %>% 
  mutate(fifth_bday = dob %m+% years(5)) %>% 
  mutate(third_bday = dob %m+% years(3)) 

obs_before_dob <- sum(obs_4in1$obsdate < obs_4in1$dob, na.rm = T)
percent_obs_before_dob <- round(obs_before_dob/total_obs*100, 3)

obs_before_3y <- sum(obs_4in1$obsdate < obs_4in1$third_bday, na.rm = T)
percent_obs_before_3y <- round(obs_before_3y/total_obs*100, 3)

obs_after5y <- sum(obs_4in1$obsdate > obs_4in1$fifth_bday, na.rm = T)
percent_obs_after5y <- round(obs_after5y/total_obs*100, 3)

obs_after_fuend <- sum(obs_4in1$obsdate > obs_4in1$C_pc_fu_end, na.rm = T)
percent_obs_after_fuend<- round(obs_after_fuend/total_obs*100, 3)

# drop obs with no obsdate, obsdate before DOB, after fifth birthday and after follow up end
obs_4in1_full <- obs_4in1 %>% 
  filter(!is.na(obsdate)) %>% 
  filter(obsdate >= dob) %>% 
  filter(obsdate >= third_bday) %>% 
  filter(obsdate <= fifth_bday) %>% 
  filter(obsdate <= C_pc_fu_end) 

# count number of med+prodcode obs after exclusions
obs_4in1_excl <- nrow(obs_4in1_full) 

# save 
save(obs_4in1_full, file="~/obs_4in1_full.Rdata")

# Cut down to 1 record per person per obsdate  ------------------------------------

# prioritise full over partial, then medcode over prodcode
# if multiple of the same priority, prioritise ones that have consid 
obs_4in1_full <- obs_4in1_full %>%
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, type == "full" & code == "prodcode", 2)) %>% 
  mutate(priority = replace(priority, type == "partial" & code == "medcode", 3)) %>% 
  mutate(priority = replace(priority, type == "partial" & code == "prodcode", 4)) %>% 
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  ungroup() %>% 
  select(-priority) %>% 
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, is.na(consid), 2)) %>%  
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  group_by(C_patid, obsdate) %>% 
  summarise(`4in1` = 1,
            consid = first(consid))
n_distinct(obs_4in1_full$C_patid) != nrow(obs_4in1_full) # should not yet be one row per patient

# Join to fifth bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/fifth_bday_patids.Rdata")

# join fifth bday 
fifth_bday_4in1 <- fifth_bday_patids %>% 
  left_join(obs_4in1_full, by = "C_patid") %>% 
  rename(obsdate_4in1 = obsdate) %>% 
  mutate(`4in1` = replace_na(`4in1`, 0))

obs_fifth_bday_4in1 <- nrow(fifth_bday_4in1) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
fifth_bday_4in1 <- fifth_bday_4in1 %>% 
  distinct(C_patid, `4in1`, obsdate_4in1, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(fifth_bday_4in1$C_patid) == n_distinct(fifth_bday_patids$C_patid)

# save
save(fifth_bday_4in1, file="~/fifth_bday_4in1.Rdata")

# Save counts and checks -------------------------------------------------------------

total_obs <- as.data.frame(total_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
obs_before_dob <- as.data.frame(obs_before_dob)
percent_obs_before_dob <- as.data.frame(percent_obs_before_dob)
obs_before_3y <- as.data.frame(obs_before_3y)
percent_obs_before_3y <- as.data.frame(percent_obs_before_3y)
obs_after5y <- as.data.frame(obs_after5y)
percent_obs_after5y <- as.data.frame(percent_obs_after5y)
obs_after_fuend <- as.data.frame(obs_after_fuend)
percent_obs_after_fuend <- as.data.frame(percent_obs_after_fuend)
obs_4in1_excl <- as.data.frame(obs_4in1_excl)

obs_fifth_bday_4in1 <- as.data.frame(obs_fifth_bday_4in1)

cleaning_4in1_counts <- cbind(total_obs, 
                              nodate_obs, percent_nodate_obs,
                              obs_before_dob, percent_obs_before_dob,
                              obs_before_3y, percent_obs_before_3y,
                              obs_after5y, percent_obs_after5y,
                              obs_after_fuend, percent_obs_after_fuend,
                              obs_4in1_excl,
                              obs_fifth_bday_4in1)

write_csv(cleaning_4in1_counts, file = "~/cleaning_4in1_counts.csv")

# Create completeness & timing variables for fifth bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
fifth_bday_4in1_completeness <- fifth_bday_4in1 %>% 
  group_by(C_patid) %>%
  summarise(`4in1` = sum(`4in1`))

# create expected dose counts
fifth_bday_4in1_completeness <- fifth_bday_4in1_completeness %>% 
  mutate(`4in1_exp_booster` = as.numeric(vaccine_exp_dose$fifth_bday_exp[vaccine_exp_dose$vaccine == "4in1"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
fifth_bday_4in1_completeness$`4in1_booster`[fifth_bday_4in1_completeness$`4in1_exp_booster`>fifth_bday_4in1_completeness$`4in1`] <- 0
fifth_bday_4in1_completeness$`4in1_booster`[fifth_bday_4in1_completeness$`4in1_exp_booster`<=fifth_bday_4in1_completeness$`4in1`] <- 1

# check for NAs
sum(is.na(fifth_bday_4in1_completeness$`4in1_booster`))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
fifth_bday_4in1_analysis <- fifth_bday_4in1_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
fifth_bday_4in1_analysis <- fifth_bday_4in1_analysis %>% 
  mutate(fifth_bday = dob %m+% years(5)) %>%
  filter(fifth_bday >= "2006-04-01") %>% 
  filter(fifth_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (fifth_bday >= "2006-04-01" & fifth_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2007-04-01" & fifth_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2008-04-01" & fifth_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2009-04-01" & fifth_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2010-04-01" & fifth_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2011-04-01" & fifth_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2012-04-01" & fifth_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2013-04-01" & fifth_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2014-04-01" & fifth_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2015-04-01" & fifth_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2016-04-01" & fifth_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2017-04-01" & fifth_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2018-04-01" & fifth_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2019-04-01" & fifth_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2020-04-01" & fifth_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(fifth_bday_4in1_analysis, file="~/fifth_bday_4in1_analysis.Rdata")

#--------------------------------------------------------------------------
# COMBINE 6/5/4-in-1                  
#--------------------------------------------------------------------------

# Load --------------------------------------------------------------------

load(file="~/first_bday_6in1_analysis.Rdata")
load(file="~/second_bday_6in1_analysis.Rdata")

load(file="~/first_bday_5in1_analysis.Rdata")
load(file="~/second_bday_5in1_analysis.Rdata")
load(file="~/fifth_bday_5in1_analysis.Rdata")

load(file="~/fifth_bday_4in1_analysis.Rdata")

# Combine pc variables ----------------------------------------

# combine first birthday pc
first_bday_6in1_analysis <- first_bday_6in1_analysis %>% 
  rename(`654in1_pc` = `6in1_pc`)

first_bday_5in1_analysis <- first_bday_5in1_analysis %>% 
  rename(`654in1_pc` = `5in1_pc`)

first_bday_654in1_analysis <- rbind(first_bday_5in1_analysis, first_bday_6in1_analysis)

# for those who were born in 2017, deduplicate and select 1 over 0 (to deal with overlap in the year that 5-in-1 switched to 6-in-1)
first_bday_654in1_analysis <- first_bday_654in1_analysis %>% 
  group_by(C_patid) %>% 
  slice(which.max(`654in1_pc`)) %>% 
  ungroup()
   
n_distinct(first_bday_654in1_analysis$C_patid) == nrow(first_bday_654in1_analysis)

# combine second birthday pc
second_bday_6in1_analysis <- second_bday_6in1_analysis %>% 
  rename(`654in1_pc` = `6in1_pc`)

second_bday_5in1_analysis <- second_bday_5in1_analysis %>% 
  rename(`654in1_pc` = `5in1_pc`)

second_bday_654in1_analysis <- rbind(second_bday_5in1_analysis, second_bday_6in1_analysis)

# for those who were born in 2017, deduplicate and select 1 over 0
second_bday_654in1_analysis <- second_bday_654in1_analysis %>% 
  group_by(C_patid) %>% 
  slice(which.max(`654in1_pc`)) %>% 
  ungroup()

n_distinct(second_bday_654in1_analysis$C_patid) == nrow(second_bday_654in1_analysis)

# combine fifth birthday pc and booster
fifth_bday_4in1_analysis <- fifth_bday_4in1_analysis %>% 
  select(C_patid, `4in1_booster`, delayed)

fifth_bday_654in1_analysis <- fifth_bday_5in1_analysis %>% 
  select(-delayed) %>% 
  left_join(fifth_bday_4in1_analysis, by = "C_patid")

# create overall completeness variable for this group of vax, i.e. whether child had completed pc AND booster
fifth_bday_654in1_analysis <- fifth_bday_654in1_analysis %>% 
  mutate(`54in1_booster` = ifelse(`4in1_booster` == 1 & `5in1_pc` == 1, 1, 0))

# save
save(first_bday_654in1_analysis, file="~/first_bday_654in1_analysis.Rdata")
save(second_bday_654in1_analysis, file="~/second_bday_654in1_analysis.Rdata")
save(fifth_bday_654in1_analysis, file="~/fifth_bday_654in1_analysis.Rdata")

#--------------------------------------------------------------------------
# MenC (12 & 16w, only 12w from June 2013, then none from July 2016)                   
#--------------------------------------------------------------------------

# Counts & exclusions -----------------------------------------------------

load(file = "~/obs_menc.Rdata")

total_obs <- nrow(obs_menc)

# count missing obsdates
nodate_obs <- sum(is.na(obs_menc$obsdate)) 
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# count obsdates before DOB, after June 2016, after fifth birthday, and after follow up end
obs_menc <- obs_menc %>% 
  mutate(fifth_bday = dob %m+% years(5)) 

obs_before_dob <- sum(obs_menc$obsdate < obs_menc$dob, na.rm = T)
percent_obs_before_dob <- round(obs_before_dob/total_obs*100, 3)

obs_after_jun2016 <- sum(obs_menc$obsdate > '2016-06-30', na.rm = T)
percent_obs_after_jun2016 <- round(obs_after_jun2016/total_obs*100, 3)

obs_after5y <- sum(obs_menc$obsdate > obs_menc$fifth_bday, na.rm = T)
percent_obs_after5y <- round(obs_after5y/total_obs*100, 3)

obs_after_fuend <- sum(obs_menc$obsdate > obs_menc$C_pc_fu_end, na.rm = T)
percent_obs_after_fuend<- round(obs_after_fuend/total_obs*100, 3)

# drop obs with no obsdate, obsdate before DOB, after June 2016, after fifth birthday and after follow up end
obs_menc_full <- obs_menc %>% 
  filter(!is.na(obsdate)) %>% 
  filter(obsdate >= dob) %>% 
  filter(obsdate < '2016-07-01') %>% 
  filter(obsdate <= fifth_bday) %>% 
  filter(obsdate <= C_pc_fu_end) 

# count number of med+prodcode obs after exclusions
obs_menc_excl <- nrow(obs_menc_full) 

# save 
save(obs_menc_full, file="~/obs_menc_full.Rdata")

# Cut down to 1 record per person per obsdate  ------------------------------------

# prioritise ones that have consid 
obs_menc_full <- obs_menc_full %>%
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, is.na(consid), 2)) %>%  
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  group_by(C_patid, obsdate) %>% 
  summarise(menc = first(menc),
            consid = first(consid))
n_distinct(obs_menc_full$C_patid) != nrow(obs_menc_full) # should not yet be one row per patient

# Join to first bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/first_bday_patids.Rdata")

# join first bday for those born before July 2016
first_bday_patids_befjul2016 <- first_bday_patids %>% 
  filter(dob < "2016-07-01") 

first_bday_menc <- first_bday_patids_befjul2016 %>% 
  left_join(obs_menc_full, by = "C_patid") %>% 
  rename(obsdate_menc = obsdate) %>% 
  mutate(menc = as.numeric(menc)) %>% 
  mutate(menc = replace_na(menc, 0))

# count obs occurring after first birthday
obs_first_bday_menc <- nrow(first_bday_menc) # denom is number of menc obs in this bday cohort 
menc_obs_after1y <- sum(first_bday_menc$obsdate_menc > first_bday_menc$first_bday, na.rm = T)
percent_menc_obs_after1y <- round(menc_obs_after1y/obs_first_bday_menc*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
first_bday_menc <- first_bday_menc %>%
  mutate(menc = replace(menc, obsdate_menc > first_bday, 0)) %>% 
  mutate(obsdate_menc = replace(obsdate_menc, obsdate_menc > first_bday, NA)) 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
first_bday_menc <- first_bday_menc %>% 
  distinct(C_patid, menc, obsdate_menc, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(first_bday_menc$C_patid) == n_distinct(first_bday_patids$C_patid[first_bday_patids$dob < "2016-07-01"])

# save
save(first_bday_menc, file="~/first_bday_menc.Rdata")

# Save counts and checks -------------------------------------------------------------

total_obs <- as.data.frame(total_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
obs_before_dob <- as.data.frame(obs_before_dob)
percent_obs_before_dob <- as.data.frame(percent_obs_before_dob)
obs_after_jun2016 <- as.data.frame(obs_after_jun2016)
percent_obs_after_jun2016 <- as.data.frame(percent_obs_after_jun2016)
obs_after5y <- as.data.frame(obs_after5y)
percent_obs_after5y <- as.data.frame(percent_obs_after5y)
obs_after_fuend <- as.data.frame(obs_after_fuend)
percent_obs_after_fuend <- as.data.frame(percent_obs_after_fuend)
obs_menc_excl <- as.data.frame(obs_menc_excl)

obs_first_bday_menc <- as.data.frame(obs_first_bday_menc)
menc_obs_after1y <- as.data.frame(menc_obs_after1y)
percent_menc_obs_after1y <- as.data.frame(percent_menc_obs_after1y)

menc_cleaning_counts <- cbind(total_obs, 
                              nodate_obs, percent_nodate_obs,
                              obs_before_dob, percent_obs_before_dob,
                              obs_after_jun2016, percent_obs_after_jun2016,
                              obs_after5y, percent_obs_after5y,
                              obs_after_fuend, percent_obs_after_fuend,
                              obs_menc_excl,
                              obs_first_bday_menc, menc_obs_after1y, percent_menc_obs_after1y)

write_csv(menc_cleaning_counts, file = "~/menc_cleaning_counts.csv")


# Create completeness & timing variables for first bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# divide cohort by when the schedule changed (June 2013)
first_bday_menc <- first_bday_menc %>% 
  mutate(schedule_change = "before") %>% 
  mutate(schedule_change = replace(schedule_change, dob >= "2013-06-01", "after"))

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
first_bday_menc_completeness <- first_bday_menc %>% 
  group_by(C_patid) %>%
  summarise(menc = sum(menc),
            schedule_change = first(schedule_change))

# create expected dose counts
first_bday_menc_completeness <- first_bday_menc_completeness %>% 
  mutate(menc_exp = as.numeric(vaccine_exp_dose$first_bday_exp[vaccine_exp_dose$vaccine == "menc_before_jun2013"])) %>% 
  mutate(menc_exp = replace(menc_exp, schedule_change == "after", as.numeric(vaccine_exp_dose$first_bday_exp[vaccine_exp_dose$vaccine == "menc_bw_jun2013_jul2016"])))

# derive primary course variable (complete = 1 vs incomplete = 0)
first_bday_menc_completeness$menc_pc[first_bday_menc_completeness$menc_exp>first_bday_menc_completeness$menc] <- 0
first_bday_menc_completeness$menc_pc[first_bday_menc_completeness$menc_exp<=first_bday_menc_completeness$menc] <- 1

# check for NAs
sum(is.na(first_bday_menc_completeness$menc_pc))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
first_bday_menc_analysis <- fifth_bday_4in1_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
first_bday_menc_analysis <- first_bday_menc_analysis %>% 
  mutate(first_bday = dob %m+% years(1)) %>%
  filter(first_bday >= "2006-04-01") %>% 
  filter(first_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (first_bday >= "2006-04-01" & first_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (first_bday >= "2007-04-01" & first_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (first_bday >= "2008-04-01" & first_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (first_bday >= "2009-04-01" & first_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (first_bday >= "2010-04-01" & first_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (first_bday >= "2011-04-01" & first_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (first_bday >= "2012-04-01" & first_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (first_bday >= "2013-04-01" & first_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (first_bday >= "2014-04-01" & first_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (first_bday >= "2015-04-01" & first_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (first_bday >= "2016-04-01" & first_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (first_bday >= "2017-04-01" & first_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (first_bday >= "2018-04-01" & first_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (first_bday >= "2019-04-01" & first_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (first_bday >= "2020-04-01" & first_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

# exclude those who had their first birthday in 2017-2018 (too few would have had the MenC vax after it was removed from the schedule)
first_bday_menc_analysis <- first_bday_menc_analysis %>% 
  filter(year != "2017-18")

save(first_bday_menc_analysis, file="~/first_bday_menc_analysis.Rdata")

#--------------------------------------------------------------------------
# Hib/MenC (1y)                   
#--------------------------------------------------------------------------

# Counts & exclusions -----------------------------------------------------

load(file = "~/obs_hibmenc.Rdata")

total_obs <- nrow(obs_hibmenc)

# count missing obsdates
nodate_obs <- sum(is.na(obs_hibmenc$obsdate)) 
percent_nodate_obs <- round(nodate_obs/total_obs*100, 3) 

# count obsdates before DOB, after fifth birthday, and after follow up end
obs_hibmenc <- obs_hibmenc %>% 
  mutate(fifth_bday = dob %m+% years(5)) 

obs_before_dob <- sum(obs_hibmenc$obsdate < obs_hibmenc$dob, na.rm = T)
percent_obs_before_dob <- round(obs_before_dob/total_obs*100, 3)

obs_after5y <- sum(obs_hibmenc$obsdate > obs_hibmenc$fifth_bday, na.rm = T)
percent_obs_after5y <- round(obs_after5y/total_obs*100, 3)

obs_after_fuend <- sum(obs_hibmenc$obsdate > obs_hibmenc$C_pc_fu_end, na.rm = T)
percent_obs_after_fuend<- round(obs_after_fuend/total_obs*100, 3)

# drop obs with no obsdate, obsdate before DOB, after fifth birthday and after follow up end
obs_hibmenc_full <- obs_hibmenc %>% 
  filter(!is.na(obsdate)) %>% 
  filter(obsdate >= dob) %>% 
  filter(obsdate <= fifth_bday) %>% 
  filter(obsdate <= C_pc_fu_end) 

# count number of med+prodcode obs after exclusions
obs_hibmenc_excl <- nrow(obs_hibmenc_full) 

# save 
save(obs_hibmenc_full, file="~/obs_hibmenc_full.Rdata")

# Cut down to 1 record per person per obsdate  ------------------------------------

# prioritise full over partial, then medcode over prodcode
# if multiple of the same priority, prioritise ones that have consid 
obs_hibmenc_full <- obs_hibmenc_full %>%
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, type == "full" & code == "prodcode", 2)) %>% 
  mutate(priority = replace(priority, type == "partial" & code == "medcode", 3)) %>% 
  mutate(priority = replace(priority, type == "partial" & code == "prodcode", 4)) %>% 
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  ungroup() %>% 
  select(-priority) %>% 
  mutate(priority = 1) %>% 
  mutate(priority = replace(priority, is.na(consid), 2)) %>%  
  group_by(C_patid, obsdate) %>% 
  slice(which.min(priority)) %>% 
  group_by(C_patid, obsdate) %>% 
  summarise(hibmenc = 1,
            consid = first(consid))
n_distinct(obs_hibmenc_full$C_patid) != nrow(obs_hibmenc_full) # should not yet be one row per patient

# Join to second bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/second_bday_patids.Rdata")

# join second bday
second_bday_hibmenc <- second_bday_patids %>% 
  left_join(obs_hibmenc_full, by = "C_patid") %>% 
  rename(obsdate_hibmenc = obsdate) %>% 
  mutate(hibmenc = replace_na(hibmenc, 0))

# count obs occurring after second birthday
obs_second_bday_hibmenc <- nrow(second_bday_hibmenc) # denom is number of hibmenc obs in this bday cohort 
hibmenc_obs_after2y <- sum(second_bday_hibmenc$obsdate_hibmenc > second_bday_hibmenc$second_bday, na.rm = T)
percent_hibmenc_obs_after2y <- round(hibmenc_obs_after2y/obs_second_bday_hibmenc*100, 3)

# drop these obs (not by removing the row but by setting the dose to 0, obsdate to NA)
second_bday_hibmenc <- second_bday_hibmenc %>%
  mutate(hibmenc = replace(hibmenc, obsdate_hibmenc > second_bday, 0)) %>% 
  mutate(obsdate_hibmenc = replace(obsdate_hibmenc, obsdate_hibmenc > second_bday, NA))

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
second_bday_hibmenc <- second_bday_hibmenc %>% 
  distinct(C_patid, hibmenc, obsdate_hibmenc, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(second_bday_hibmenc$C_patid) == n_distinct(second_bday_patids$C_patid)

# save
save(second_bday_hibmenc, file="~/second_bday_hibmenc.Rdata")

# Join to fifth bday patid list -----------------------------------------------

# load bday cohort patid file
load(file="~/fifth_bday_patids.Rdata")

# join fifth bday
fifth_bday_hibmenc <- fifth_bday_patids %>% 
  left_join(obs_hibmenc_full, by = "C_patid") %>% 
  rename(obsdate_hibmenc = obsdate) %>% 
  mutate(hibmenc = replace_na(hibmenc, 0))

obs_fifth_bday_hibmenc <- nrow(fifth_bday_hibmenc) # denom is number of hibmenc obs in this bday cohort 

# remove duplicates where a patient has multiple rows of vaccine dose = 0 and obsdate = NA
fifth_bday_hibmenc <- fifth_bday_hibmenc %>% 
  distinct(C_patid, hibmenc, obsdate_hibmenc, .keep_all = TRUE)

# check distinct patients still the same as the original cohort file before obs were added
n_distinct(fifth_bday_hibmenc$C_patid) == n_distinct(fifth_bday_patids$C_patid)

# save
save(fifth_bday_hibmenc, file="~/fifth_bday_hibmenc.Rdata")

# Save counts and checks -------------------------------------------------------------

total_obs <- as.data.frame(total_obs)
nodate_obs <- as.data.frame(nodate_obs)
percent_nodate_obs <- as.data.frame(percent_nodate_obs)
obs_before_dob <- as.data.frame(obs_before_dob)
percent_obs_before_dob <- as.data.frame(percent_obs_before_dob)
obs_after5y <- as.data.frame(obs_after5y)
percent_obs_after5y <- as.data.frame(percent_obs_after5y)
obs_after_fuend <- as.data.frame(obs_after_fuend)
percent_obs_after_fuend <- as.data.frame(percent_obs_after_fuend)
obs_hibmenc_excl <- as.data.frame(obs_hibmenc_excl)

obs_second_bday_hibmenc <- as.data.frame(obs_second_bday_hibmenc)
hibmenc_obs_after2y <- as.data.frame(hibmenc_obs_after2y)
percent_hibmenc_obs_after2y <- as.data.frame(percent_hibmenc_obs_after2y)

obs_fifth_bday_hibmenc <- as.data.frame(obs_fifth_bday_hibmenc)

hibmenc_cleaning_counts <- cbind(total_obs, 
                             nodate_obs, percent_nodate_obs,
                             obs_before_dob, percent_obs_before_dob,
                             obs_after5y, percent_obs_after5y,
                             obs_after_fuend, percent_obs_after_fuend,
                             obs_hibmenc_excl,
                             obs_second_bday_hibmenc, hibmenc_obs_after2y, percent_hibmenc_obs_after2y,
                             obs_fifth_bday_hibmenc)

write_csv(hibmenc_cleaning_counts, file = "~/hibmenc_cleaning_counts.csv")

# Create completeness & timing variables for second bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
second_bday_hibmenc_completeness <- second_bday_hibmenc %>% 
  group_by(C_patid) %>%
  summarise(hibmenc = sum(hibmenc))

# create expected dose counts
second_bday_hibmenc_completeness <- second_bday_hibmenc_completeness %>% 
  mutate(hibmenc_exp = as.numeric(vaccine_exp_dose$second_bday_exp[vaccine_exp_dose$vaccine == "hibmenc"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
second_bday_hibmenc_completeness$hibmenc_booster[second_bday_hibmenc_completeness$hibmenc_exp>second_bday_hibmenc_completeness$hibmenc] <- 0
second_bday_hibmenc_completeness$hibmenc_booster[second_bday_hibmenc_completeness$hibmenc_exp<=second_bday_hibmenc_completeness$hibmenc] <- 1

# check for NAs
sum(is.na(second_bday_hibmenc_completeness$hibmenc_booster))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
second_bday_hibmenc_analysis <- second_bday_hibmenc_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
second_bday_hibmenc_analysis <- second_bday_hibmenc_analysis %>% 
  mutate(second_bday = dob %m+% years(2)) %>%
  filter(second_bday >= "2006-04-01") %>% 
  filter(second_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (second_bday >= "2006-04-01" & second_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (second_bday >= "2007-04-01" & second_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (second_bday >= "2008-04-01" & second_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (second_bday >= "2009-04-01" & second_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (second_bday >= "2010-04-01" & second_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (second_bday >= "2011-04-01" & second_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (second_bday >= "2012-04-01" & second_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (second_bday >= "2013-04-01" & second_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (second_bday >= "2014-04-01" & second_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (second_bday >= "2015-04-01" & second_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (second_bday >= "2016-04-01" & second_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (second_bday >= "2017-04-01" & second_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (second_bday >= "2018-04-01" & second_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (second_bday >= "2019-04-01" & second_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (second_bday >= "2020-04-01" & second_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(second_bday_hibmenc_analysis, file="~/second_bday_hibmenc_analysis.Rdata")

# Create completeness & timing variables for fifth bday cohort ---------------------------------------------------

vaccine_exp_dose <- read_csv('~/vaccine_exp_dose.csv', 
                             col_names = TRUE)

# to derive completeness, group by patid to get one row per patient and a count of all observed doses up to this timepoint 
fifth_bday_hibmenc_completeness <- fifth_bday_hibmenc %>% 
  group_by(C_patid) %>%
  summarise(hibmenc = sum(hibmenc))

# create expected dose counts
fifth_bday_hibmenc_completeness <- fifth_bday_hibmenc_completeness %>% 
  mutate(hibmenc_exp_pc = as.numeric(vaccine_exp_dose$second_bday_exp[vaccine_exp_dose$vaccine == "hibmenc"])) %>% 
  mutate(hibmenc_exp_booster = as.numeric(vaccine_exp_dose$fifth_bday_exp[vaccine_exp_dose$vaccine == "hibmenc"]))

# derive primary course variable (complete = 1 vs incomplete = 0)
fifth_bday_hibmenc_completeness$hibmenc_booster[fifth_bday_hibmenc_completeness$hibmenc_exp_pc>fifth_bday_hibmenc_completeness$hibmenc] <- 0
fifth_bday_hibmenc_completeness$hibmenc_booster[fifth_bday_hibmenc_completeness$hibmenc_exp_pc<=fifth_bday_hibmenc_completeness$hibmenc] <- 1

# check for NAs
sum(is.na(fifth_bday_hibmenc_completeness$hibmenc_booster))

# add rest of variables needed for the analysis
load(file="~/sourcepop.Rdata")
fifth_bday_hibmenc_analysis <- fifth_bday_hibmenc_completeness %>% 
  left_join(sourcepop, by = "C_patid")

# create financial year in which birthday took place and differentiate it from financial year of birth
# divide bday cohort into years (i.e. children whose bday happened within that year)
fifth_bday_hibmenc_analysis <- fifth_bday_hibmenc_analysis %>% 
  mutate(fifth_bday = dob %m+% years(5)) %>%
  filter(fifth_bday >= "2006-04-01") %>% 
  filter(fifth_bday <= "2021-03-31") %>% 
  mutate(year = " ") %>% 
  mutate(year = replace(year, (fifth_bday >= "2006-04-01" & fifth_bday < "2007-04-01"), "2006-07")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2007-04-01" & fifth_bday < "2008-04-01"), "2007-08")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2008-04-01" & fifth_bday < "2009-04-01"), "2008-09")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2009-04-01" & fifth_bday < "2010-04-01"), "2009-10")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2010-04-01" & fifth_bday < "2011-04-01"), "2010-11")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2011-04-01" & fifth_bday < "2012-04-01"), "2011-12")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2012-04-01" & fifth_bday < "2013-04-01"), "2012-13")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2013-04-01" & fifth_bday < "2014-04-01"), "2013-14")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2014-04-01" & fifth_bday < "2015-04-01"), "2014-15")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2015-04-01" & fifth_bday < "2016-04-01"), "2015-16")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2016-04-01" & fifth_bday < "2017-04-01"), "2016-17")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2017-04-01" & fifth_bday < "2018-04-01"), "2017-18")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2018-04-01" & fifth_bday < "2019-04-01"), "2018-19")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2019-04-01" & fifth_bday < "2020-04-01"), "2019-20")) %>% 
  mutate(year = replace(year, (fifth_bday >= "2020-04-01" & fifth_bday < "2021-04-01"), "2020-21")) %>% 
  mutate(year = as.factor(year)) %>% 
  rename(financialyob = financialyear)

save(fifth_bday_hibmenc_analysis, file="~/fifth_bday_hibmenc_analysis.Rdata")

