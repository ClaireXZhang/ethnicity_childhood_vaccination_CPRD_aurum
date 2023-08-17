
#--------------------------------------------------------------------------
# SET UP                          
#--------------------------------------------------------------------------
# Description -------------------------------------------------------------

# Validate and analyse routine childhood vaccinations in CPRD Aurum (except influenza)
# Example: MMR (repeat on all other vaccines using the appropriate birthday cohorts)

# Set working directory ---------------------------------------------------

setwd('X:/')

# Load packages -----------------------------------------------------------

pacman::p_load(tidyverse, lubridate, 
               forestplot, patchwork, ggplotify,
               epitools, RColorBrewer, pals,
               gghighlight, # highlighting specific groups in graphs
               lmtest, sandwich) # modified Poisson, sandwich and clustering

# Disable scientific notation ---------------------------------------------

options(scipen = 999)

# Write functions ----------------------------------------------------------

# create results table from modified poisson regression - adapted from Zhang, Boukari et al 2021
extract_mpoisson_results <- function (x, x_sandwich, output_name, dataset, variable) {
  ethnicity <- names(coef(x))
  ethnicity <- sub("M_eth18_2011", "", ethnicity) # extracts all text after ')' from ethnicity to avoid having to recode later (excluding the intercept, the below lines of code are deal with this)
  ethnicity[1] <- levels(eval(substitute(variable), dataset))[1] # takes the reference level of a specified factor variable and replaces it with 'Intercept' in ethnicity
  estimate<- x_sandwich$RR
  lower <- x_sandwich$LCI
  upper <- x_sandwich$UCI
  p <- x_sandwich$P
  output <- cbind(ethnicity,estimate,lower, upper,p)
  output <- tibble::as_tibble(output)
  output$estimate <- as.numeric(output$estimate)
  output$upper <- as.numeric(output$upper)
  output$lower <- as.numeric(output$lower)
  output$p <- as.numeric(output$p)
  output <- output %>% dplyr::mutate(across(where(is.numeric), ~ round(.,2)))
  output$ci <- paste(output$lower, output$upper, sep ="-")
  output$estimate[1] <- 1.00
  output$lower[1] <- 1.00
  output$upper[1] <- 1.00
  output$ci[1] <- "1.00-1.00"
  output$rr_ci <- paste0(output$estimate," (", output$ci,")")
  output_table <- dplyr::select(output,ethnicity, rr_ci )
  assign(x = output_name, value = output, envir = globalenv()) # changes the name of the df output to the name you want and saves it in the global environment
  assign(x = paste0(output_name,'_table'), value = output_table, envir = globalenv()) # changes the name of the df output to the name you want and saves it in the global environment
}

# calculate stratum specific CIs using sandwich estimators
get_sandwich_stratum_ci <- function (x_sandwich, ethnicity_est, ethnicity_year_int) {
  round(exp(confint(x_sandwich, level = 0.95)[ethnicity_est,] + confint(x_sandwich, level = 0.95)[ethnicity_year_int,]), 2)
}

#--------------------------------------------------------------------------
# MMR - DESCRIPTIVE ANALYSIS                   
#--------------------------------------------------------------------------

load(file="~/second_bday_mmr_analysis.Rdata")
load(file="~/fifth_bday_mmr_analysis.Rdata")

# Validation: Comparison to UKHSA/NHSE routine statistical reports -------------------------------------------------

# load nhs reference dataset
nhs_ref <- read_csv("~/ukhsa_nhse_data_clean.csv")

# filter by mmr and available years and bday cohorts
nhs_mmr <- nhs_ref %>% 
  filter(vaccine == "mmr") %>% 
  filter(!is.na(pc_nhs) | !is.na(booster_nhs)) %>% 
  select(-vaccine)

# create denominators for each year
mmr_second_bday_nhs_denom <- second_bday_mmr_analysis %>% 
  group_by(year) %>% 
  count() %>% 
  rename(N = n)

# count number who have completed primary course (i.e. all expected doses) for each year
mmr_second_pc_nhs <- second_bday_mmr_analysis %>% 
  group_by(year, mmr_pc) %>% 
  count() %>% 
  ungroup() %>% 
  filter(mmr_pc == 1) %>% 
  dplyr::select(-mmr_pc) 

# join and create percentages
mmr_second_pc_nhs <- mmr_second_pc_nhs %>% 
  left_join(mmr_second_bday_nhs_denom, by = c("year")) %>% 
  mutate(pc_cprd = round(n/N*100, 1)) %>% 
  dplyr::select(-c(n, N)) %>% 
  mutate(bday = "second") 
mmr_second_pc_nhs$booster_cprd <- NA

# create denominators for each year
mmr_fifth_bday_nhs_denom <- fifth_bday_mmr_analysis %>% 
  group_by(year) %>% 
  count() %>% 
  rename(N = n)

# count number who have completed primary course for each year
mmr_fifth_pc_nhs <- fifth_bday_mmr_analysis %>% 
  group_by(year, mmr_pc) %>% 
  count() %>% 
  ungroup() %>% 
  filter(mmr_pc == 1) %>% 
  select(-mmr_pc) 

# count number who have completed booster for each year
mmr_fifth_booster_nhs <- fifth_bday_mmr_analysis %>% 
  group_by(year, mmr_booster) %>% 
  count() %>% 
  ungroup() %>% 
  filter(mmr_booster == 1) %>% 
  select(-mmr_booster) 

# join and create percentages
mmr_fifth_bday_nhs <- mmr_fifth_pc_nhs %>% 
  left_join(mmr_fifth_bday_nhs_denom, by = c("year")) %>% 
  mutate(pc_cprd = round(n/N*100, 1)) %>% 
  dplyr::select(-c(n)) %>% 
  left_join(mmr_fifth_booster_nhs, by = c("year")) %>% 
  mutate(booster_cprd = round(n/N*100, 1)) %>% 
  dplyr::select(-c(n, N)) %>% 
  mutate(bday = "fifth")

# join to nhs reference data
mmr_second_fifth_nhs <- rbind(mmr_second_pc_nhs, mmr_fifth_bday_nhs)
mmr_second_fifth_nhs <- mmr_second_fifth_nhs %>% 
  left_join(nhs_mmr, by = c("year", "bday")) 

# create columns for % difference
mmr_second_fifth_nhs <- mmr_second_fifth_nhs %>% 
  mutate(pc_diff = pc_cprd - pc_nhs) %>% 
  mutate(booster_diff = booster_cprd - booster_nhs) %>% 
  select(year, bday, pc_cprd, pc_nhs, pc_diff, booster_cprd, booster_nhs, booster_diff)

write_csv(mmr_second_fifth_nhs, file="~/mmr_second_fifth_nhs.csv")

# Coverage - all years, all courses -----------------------------------------------------------

# Get percentages and 95%CIs (using this method: https://sphweb.bumc.bu.edu/otlt/MPH-Modules/PH717-QuantCore/PH717-Module6-RandomError/PH717-Module6-RandomError12.html)

# second bday - primary course
second_mmr_pc_total <- nrow(second_bday_mmr_analysis)
second_mmr_pc_all <- second_bday_mmr_analysis %>% 
  group_by(mmr_pc) %>% 
  count() %>% 
  mutate(percent = round(n/second_mmr_pc_total*100, 1)) %>% 
  mutate(prop = n/second_mmr_pc_total) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/second_mmr_pc_total))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/second_mmr_pc_total))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(mmr_pc == 1) %>% 
  mutate(M_eth18_2011 = "All")
second_mmr_pc_all_table <- second_mmr_pc_all %>% 
  mutate(`Second birthday: n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  ungroup() %>% 
  select(`Second birthday: n (%, 95%CI)`, M_eth18_2011)  

second_mmr_pc_denom <- second_bday_mmr_analysis %>% 
  group_by(M_eth18_2011) %>% 
  count() %>% 
  rename(N = n)
second_mmr_pc <- second_bday_mmr_analysis %>% 
  group_by(mmr_pc, M_eth18_2011) %>% 
  count() %>% 
  left_join(second_mmr_pc_denom, by = "M_eth18_2011") %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(mmr_pc == 1) 
second_mmr_pc_table <- second_mmr_pc %>% 
  mutate(`Second birthday: n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  ungroup() %>% 
  select(M_eth18_2011, `Second birthday: n (%, 95%CI)`)

second_mmr_pc <- rbind(second_mmr_pc_all, second_mmr_pc)
second_mmr_pc_table <- rbind(second_mmr_pc_all_table, second_mmr_pc_table) %>% 
  select(M_eth18_2011, `Second birthday: n (%, 95%CI)`)

write_csv(second_mmr_pc, file="~/second_mmr_pc.csv") 
write_csv(second_mmr_pc_table, file="~/second_mmr_pc_table.csv") 

# fifth bday - primary course
fifth_mmr_pc_total <- nrow(fifth_bday_mmr_analysis)
fifth_mmr_pc_all <- fifth_bday_mmr_analysis %>% 
  group_by(mmr_pc) %>% 
  count() %>% 
  mutate(percent = round(n/fifth_mmr_pc_total*100, 1)) %>% 
  mutate(prop = n/fifth_mmr_pc_total) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/fifth_mmr_pc_total))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/fifth_mmr_pc_total))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(mmr_pc == 1) %>% 
  mutate(M_eth18_2011 = "All")
fifth_mmr_pc_all_table <- fifth_mmr_pc_all %>% 
  mutate(`Fifth birthday: n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  ungroup() %>% 
  select(`Fifth birthday: n (%, 95%CI)`, M_eth18_2011)

fifth_mmr_pc_denom <- fifth_bday_mmr_analysis %>% 
  group_by(M_eth18_2011) %>% 
  count() %>% 
  rename(N = n)
fifth_mmr_pc <- fifth_bday_mmr_analysis %>% 
  group_by(mmr_pc, M_eth18_2011) %>% 
  count() %>% 
  left_join(fifth_mmr_pc_denom, by = "M_eth18_2011") %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(mmr_pc == 1) 
fifth_mmr_pc_table <- fifth_mmr_pc %>% 
  mutate(`Fifth birthday: n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  ungroup() %>% 
  select(M_eth18_2011, `Fifth birthday: n (%, 95%CI)`)

fifth_mmr_pc <- rbind(fifth_mmr_pc_all, fifth_mmr_pc) 
fifth_mmr_pc_table <- rbind(fifth_mmr_pc_all_table, fifth_mmr_pc_table) %>% 
  select(M_eth18_2011, `Fifth birthday: n (%, 95%CI)`)

write_csv(fifth_mmr_pc, file="~/fifth_mmr_pc.csv") 
write_csv(fifth_mmr_pc_table, file="~/fifth_mmr_pc_table.csv") 

# fifth bday - booster
fifth_mmr_booster_total <- nrow(fifth_bday_mmr_analysis)
fifth_mmr_booster_all <- fifth_bday_mmr_analysis %>% 
  group_by(mmr_booster) %>% 
  count() %>% 
  mutate(percent = round(n/fifth_mmr_booster_total*100, 1)) %>% 
  mutate(prop = n/fifth_mmr_booster_total) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/fifth_mmr_booster_total))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/fifth_mmr_booster_total))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(mmr_booster == 1) %>% 
  mutate(M_eth18_2011 = "All")
fifth_mmr_booster_all_table <- fifth_mmr_booster_all %>% 
  mutate(`Fifth birthday: n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  ungroup() %>% 
  select(`Fifth birthday: n (%, 95%CI)`, M_eth18_2011)

fifth_mmr_booster_denom <- fifth_bday_mmr_analysis %>% 
  group_by(M_eth18_2011) %>% 
  count() %>% 
  rename(N = n)
fifth_mmr_booster <- fifth_bday_mmr_analysis %>% 
  group_by(mmr_booster, M_eth18_2011) %>% 
  count() %>% 
  left_join(fifth_mmr_booster_denom, by = "M_eth18_2011") %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  filter(mmr_booster == 1) 
fifth_mmr_booster_table <- fifth_mmr_booster %>% 
  mutate(`Fifth birthday: n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  ungroup() %>% 
  select(M_eth18_2011, `Fifth birthday: n (%, 95%CI)`)

fifth_mmr_booster <- rbind(fifth_mmr_booster_all, fifth_mmr_booster) 
fifth_mmr_booster_table <- rbind(fifth_mmr_booster_all_table, fifth_mmr_booster_table) %>% 
  select(M_eth18_2011, `Fifth birthday: n (%, 95%CI)`)

write_csv(fifth_mmr_booster, file="~/fifth_mmr_booster.csv") 
write_csv(fifth_mmr_booster_table, file="~/fifth_mmr_booster_table.csv") 

# combine all results
second_mmr_pc <- second_mmr_pc %>% 
  mutate(`Birthday & course` = "Second birthday primary course")
fifth_mmr_pc <- fifth_mmr_pc %>% 
  mutate(`Birthday & course` = "Fifth birthday primary course")
fifth_mmr_booster <- fifth_mmr_booster %>% 
  mutate(`Birthday & course` = "Fifth birthday booster")
allbday_mmr <- rbind(second_mmr_pc, fifth_mmr_pc, fifth_mmr_booster) %>% 
  mutate(`Birthday & course` = as.factor(`Birthday & course`)) %>% 
  mutate(`Birthday & course` = fct_relevel(`Birthday & course`, 
                                           "Second birthday primary course",
                                           "Fifth birthday primary course",
                                           "Fifth birthday booster")) %>% 
  mutate(M_eth18_2011 = as.factor(M_eth18_2011)) %>% 
  mutate(M_eth18_2011 = fct_relevel(M_eth18_2011, "All",
                                    "English, Welsh, Scottish, Northern Irish or British", "Irish", "Any other White background",  
                                    "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background",
                                    "Caribbean", "African", "Any other Black, African or Caribbean background",
                                    "White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed or multiple ethnic background",
                                    "Any other ethnic group", "Unknown"))

# relabel ethnicity
allbday_mmr <- allbday_mmr %>% 
  rename(`Mother's ethnicity` = M_eth18_2011)

# dot plots 
pd <- position_dodge(0.75) 
allbday_mmr_dotplot <- ggplot(allbday_mmr, aes(x = `Birthday & course`, y = percent, colour = `Mother's ethnicity`)) + 
  geom_linerange(aes(ymin=lowerci, ymax=upperci), position=pd) +
  theme_bw() +
  geom_point(aes(shape = `Mother's ethnicity`, colour = `Mother's ethnicity`), position=pd, size = 4) +
  scale_shape_manual(values = c(16, 8, 17, 0, 16, 8, 17, 0, 18, 16, 8, 17,16, 8, 17, 0, 16, 8)) +
  scale_color_manual(values = c("All" = "#564D50",
                                "English, Welsh, Scottish, Northern Irish or British" = "#970055",
                                "Irish" = "#970055",
                                "Any other White background" = "#970055",
                                "Indian" = "#F5B521",
                                "Pakistani" = "#F5B521",
                                "Bangladeshi" = "#F5B521",
                                "Chinese" = "#F5B521",
                                "Any other Asian background" = "#F5B521",
                                "Caribbean" = "#158754",
                                "African" = "#158754",
                                "Any other Black, African or Caribbean background" = "#158754",
                                "White and Black Caribbean" = "#2F84FF",
                                "White and Black African" = "#2F84FF",
                                "White and Asian" = "#2F84FF" ,
                                "Any other Mixed or multiple ethnic background" = "#2F84FF",
                                "Any other ethnic group" = "#A438CD",
                                "Unknown" = "#564D50")) +
  theme(text = element_text(size =20)) +
  labs(y = "MMR vaccination (%)", x = "") +
  ylim(65, 100)
allbday_mmr_dotplot

ggsave(width = 18, height = 8, dpi = 450, "~/allbday_mmr_dotplot.jpg")

# Coverage by year of vaccination - primary course -----------------------------------------------------------

## second year - primary course
second_mmr_pc_total <- nrow(second_bday_mmr_analysis)

# create denominator for each year
mmr_second_pc_year_denom_all <- second_bday_mmr_analysis %>% 
  group_by(year) %>% 
  count() %>% 
  rename(N = n)

# count number who have completed primary course for each year
mmr_second_pc_year_all <- second_bday_mmr_analysis %>% 
  group_by(year, mmr_pc) %>% 
  count() %>% 
  ungroup() %>% 
  filter(mmr_pc == 1) %>% 
  select(-mmr_pc) 

# join and create percentages
mmr_second_pc_year_all <- mmr_second_pc_year_all %>% 
  left_join(mmr_second_pc_year_denom_all, by = "year") %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  mutate(M_eth18_2011 = "All")
mmr_second_pc_year_all_table <- mmr_second_pc_year_all %>% 
  mutate(`n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  select(M_eth18_2011, year, `n (%, 95%CI)`) 

# create denominators for each year - by ethnicity
mmr_second_pc_year_denom <- second_bday_mmr_analysis %>% 
  group_by(M_eth18_2011, year) %>% 
  count() %>% 
  rename(N = n)

# count number who have completed primary course for each year
mmr_second_pc_year <- second_bday_mmr_analysis %>% 
  group_by(M_eth18_2011, year, mmr_pc) %>% 
  count() %>% 
  ungroup() %>% 
  filter(mmr_pc == 1) %>% 
  select(-mmr_pc) 

# join and create percentages
mmr_second_pc_year <- mmr_second_pc_year %>% 
  left_join(mmr_second_pc_year_denom, by = c("M_eth18_2011", "year")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) 
mmr_second_pc_year_table <- mmr_second_pc_year %>% 
  mutate(`n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  select(M_eth18_2011, year, `n (%, 95%CI)`) 
  
mmr_second_pc_year <- rbind(mmr_second_pc_year_all, mmr_second_pc_year) 
mmr_second_pc_year_table <- rbind(mmr_second_pc_year_all_table, mmr_second_pc_year_table) %>% 
  select(M_eth18_2011,  year,`n (%, 95%CI)`)

write_csv(mmr_second_pc_year, file="~/mmr_second_pc_year.csv") 
write_csv(mmr_second_pc_year_table, file="~/mmr_second_pc_year_table.csv") 

# create NAs for years that do not have any data
mmr_second_pc_year <- mmr_second_pc_year %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "All") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "English, Welsh, Scottish, Northern Irish or British") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Irish") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other White background") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Indian") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Pakistani") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Bangladeshi") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Chinese") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Asian background") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Caribbean") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "African") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Black, African or Caribbean background") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Black Caribbean") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Black African") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Asian") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Mixed or multiple ethnic background") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other ethnic group") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Unknown")  

# relabel ethnicity and set factor level
mmr_second_pc_year <- mmr_second_pc_year %>% 
  mutate(M_eth18_2011 = as.factor(M_eth18_2011)) %>% 
  mutate(M_eth18_2011 = fct_relevel(M_eth18_2011, "All",
                                    "English, Welsh, Scottish, Northern Irish or British", "Irish", "Any other White background",  
                                    "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background",
                                    "Caribbean", "African", "Any other Black, African or Caribbean background",
                                    "White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed or multiple ethnic background",
                                    "Any other ethnic group", "Unknown")) %>% 
  rename(`Mother's ethnicity` = M_eth18_2011) 

# dot plots 
mmr_pc_year_dotplot <- ggplot(mmr_second_pc_year, aes(x = year, y = percent, colour = `Mother's ethnicity`)) + 
  geom_line(aes(group = `Mother's ethnicity`), size = 0.5) +
  theme_bw() +
  geom_point(aes(shape = `Mother's ethnicity`, colour = `Mother's ethnicity`), size = 4) +
  scale_shape_manual(values = c(16, 8, 17, 0, 16, 8, 17, 0, 18, 16, 8, 17,16, 8, 17, 0, 16, 8)) +
  scale_color_manual(values = c("All" = "#564D50",
                                "English, Welsh, Scottish, Northern Irish or British" = "#970055",
                                "Irish" = "#970055",
                                "Any other White background" = "#970055",
                                "Indian" = "#F5B521",
                                "Pakistani" = "#F5B521",
                                "Bangladeshi" = "#F5B521",
                                "Chinese" = "#F5B521",
                                "Any other Asian background" = "#F5B521",
                                "Caribbean" = "#158754",
                                "African" = "#158754",
                                "Any other Black, African or Caribbean background" = "#158754",
                                "White and Black Caribbean" = "#2F84FF",
                                "White and Black African" = "#2F84FF",
                                "White and Asian" = "#2F84FF" ,
                                "Any other Mixed or multiple ethnic background" = "#2F84FF",
                                "Any other ethnic group" = "#A438CD",
                                "Unknown" = "#564D50")) +
  theme(text = element_text(size =20)) +
  labs(y = "Received primary course of MMR (%)", x = "Child's second birthday (financial year)") +
  ylim(60, 100)
mmr_pc_year_dotplot

ggsave(width = 22, height = 8, dpi = 450, "~/mmr_pc_year_dotplot.jpg")

# Coverage by year of vaccination - booster -----------------------------------------------------------

## fifth year - booster
fifth_mmr_booster_total <- nrow(fifth_bday_mmr_analysis)

# create denominator for each year
mmr_fifth_booster_year_denom_all <- fifth_bday_mmr_analysis %>% 
  group_by(year) %>% 
  count() %>% 
  rename(N = n)

# count number who have completed primary course for each year
mmr_fifth_booster_year_all <- fifth_bday_mmr_analysis %>% 
  group_by(year, mmr_booster) %>% 
  count() %>% 
  ungroup() %>% 
  filter(mmr_booster == 1) %>% 
  select(-mmr_booster) 

# join and create percentages
mmr_fifth_booster_year_all <- mmr_fifth_booster_year_all %>% 
  left_join(mmr_fifth_booster_year_denom_all, by = "year") %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) %>% 
  mutate(M_eth18_2011 = "All")
mmr_fifth_booster_year_all_table <- mmr_fifth_booster_year_all %>% 
  mutate(`n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  select(M_eth18_2011, year, `n (%, 95%CI)`) 

# create denominators for each year - by ethnicity
mmr_fifth_booster_year_denom <- fifth_bday_mmr_analysis %>% 
  group_by(M_eth18_2011, year) %>% 
  count() %>% 
  rename(N = n)

# count number who have completed primary course for each year
mmr_fifth_booster_year <- fifth_bday_mmr_analysis %>% 
  group_by(M_eth18_2011, year, mmr_booster) %>% 
  count() %>% 
  ungroup() %>% 
  filter(mmr_booster == 1) %>% 
  select(-mmr_booster) 

# join and create percentages
mmr_fifth_booster_year <- mmr_fifth_booster_year %>% 
  left_join(mmr_fifth_booster_year_denom, by = c("M_eth18_2011", "year")) %>% 
  mutate(prop = n/N) %>% 
  mutate(percent = round(n/N*100, 1)) %>% 
  mutate(lowerci = prop-1.96*sqrt(prop*((1-prop)/N))) %>%
  mutate(upperci = prop+1.96*sqrt(prop*((1-prop)/N))) %>% 
  mutate(lowerci = round(lowerci*100, 1)) %>% 
  mutate(upperci = round(upperci*100, 1)) 
mmr_fifth_booster_year_table <- mmr_fifth_booster_year %>% 
  mutate(`n (%, 95%CI)` = paste0(n," (",percent,", ",lowerci,"-",upperci,")")) %>% 
  select(M_eth18_2011, year, `n (%, 95%CI)`) 

mmr_fifth_booster_year <- rbind(mmr_fifth_booster_year_all, mmr_fifth_booster_year) 
mmr_fifth_booster_year_table <- rbind(mmr_fifth_booster_year_all_table, mmr_fifth_booster_year_table) %>% 
  select(M_eth18_2011,  year,`n (%, 95%CI)`)

write_csv(mmr_fifth_booster_year, file="~/mmr_fifth_booster_year.csv") 
write_csv(mmr_fifth_booster_year_table, file="~/mmr_fifth_booster_year_table.csv") 

# create NAs for years that do not have any data
mmr_fifth_booster_year <- mmr_fifth_booster_year %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "All") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "English, Welsh, Scottish, Northern Irish or British") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Irish") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other White background") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Indian") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Pakistani") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Bangladeshi") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Chinese") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Asian background") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Caribbean") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "African") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Black, African or Caribbean background") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Black Caribbean") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Black African") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Asian") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Mixed or multiple ethnic background") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other ethnic group") %>% 
  add_row(year = "2007-08", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Unknown") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "All") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "English, Welsh, Scottish, Northern Irish or British") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Irish") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other White background") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Indian") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Pakistani") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Bangladeshi") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Chinese") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Asian background") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Caribbean") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "African") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Black, African or Caribbean background") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Black Caribbean") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Black African") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Asian") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Mixed or multiple ethnic background") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other ethnic group") %>% 
  add_row(year = "2008-09", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Unknown")  %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "All") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "English, Welsh, Scottish, Northern Irish or British") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Irish") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other White background") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Indian") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Pakistani") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Bangladeshi") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Chinese") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Asian background") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Caribbean") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "African") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Black, African or Caribbean background") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Black Caribbean") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Black African") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Asian") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Mixed or multiple ethnic background") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other ethnic group") %>% 
  add_row(year = "2009-10", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Unknown")  %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "All") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "English, Welsh, Scottish, Northern Irish or British") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Irish") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other White background") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Indian") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Pakistani") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Bangladeshi") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Chinese") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Asian background") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Caribbean") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "African") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Black, African or Caribbean background") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Black Caribbean") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Black African") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "White and Asian") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other Mixed or multiple ethnic background") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Any other ethnic group") %>% 
  add_row(year = "2010-11", n = NA, N = NA, prop = NA, percent = NA, lowerci = NA, upperci = NA, M_eth18_2011 = "Unknown")  

# relabel ethnicity and set factor level
mmr_fifth_booster_year <- mmr_fifth_booster_year %>% 
  mutate(M_eth18_2011 = as.factor(M_eth18_2011)) %>% 
  mutate(M_eth18_2011 = fct_relevel(M_eth18_2011, "All",
                                    "English, Welsh, Scottish, Northern Irish or British", "Irish", "Any other White background",  
                                    "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background",
                                    "Caribbean", "African", "Any other Black, African or Caribbean background",
                                    "White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed or multiple ethnic background",
                                    "Any other ethnic group", "Unknown")) %>% 
  rename(`Mother's ethnicity` = M_eth18_2011) 

# dot plots 
mmr_booster_year_dotplot <- ggplot(mmr_fifth_booster_year, aes(x = year, y = percent, colour = `Mother's ethnicity`)) + 
  geom_line(aes(group = `Mother's ethnicity`), size = 0.5) +
  theme_bw() +
  geom_point(aes(shape = `Mother's ethnicity`, colour = `Mother's ethnicity`), size = 4) +
  scale_shape_manual(values = c(16, 8, 17, 0, 16, 8, 17, 0, 18, 16, 8, 17,16, 8, 17, 0, 16, 8)) +
  scale_color_manual(values = c("All" = "#564D50",
                                "English, Welsh, Scottish, Northern Irish or British" = "#970055",
                                "Irish" = "#970055",
                                "Any other White background" = "#970055",
                                "Indian" = "#F5B521",
                                "Pakistani" = "#F5B521",
                                "Bangladeshi" = "#F5B521",
                                "Chinese" = "#F5B521",
                                "Any other Asian background" = "#F5B521",
                                "Caribbean" = "#158754",
                                "African" = "#158754",
                                "Any other Black, African or Caribbean background" = "#158754",
                                "White and Black Caribbean" = "#2F84FF",
                                "White and Black African" = "#2F84FF",
                                "White and Asian" = "#2F84FF" ,
                                "Any other Mixed or multiple ethnic background" = "#2F84FF",
                                "Any other ethnic group" = "#A438CD",
                                "Unknown" = "#564D50")) +
  theme(text = element_text(size =20)) +
  labs(y = "Received MMR booster (%)", x = "Child's fifth birthday (financial year)") +
  ylim(60, 100)
mmr_booster_year_dotplot

ggsave(width = 22, height = 8, dpi = 450, "~/mmr_booster_year_dotplot.jpg")

# Risk differences by year of vaccination - primary course  ---------------------------------------------------------------------

# load prop by year
second_mmr_pc <- read_csv(file="~/second_mmr_pc.csv") 
mmr_second_pc_year <- read_csv(file="~/mmr_second_pc_year.csv") 

# restructure to remove 'all' and join ref to rest
second_mmr_pc_groups <- second_mmr_pc %>% 
  filter(M_eth18_2011 != "All") %>% 
  filter(M_eth18_2011 != "English, Welsh, Scottish, Northern Irish or British") %>% 
  mutate(year = "All") %>% 
  select(year, prop, N, M_eth18_2011) %>% 
  rename(n2 = N) %>% 
  rename(p2 = prop)
second_mmr_pc_ref <- second_mmr_pc %>% 
  filter(M_eth18_2011 == "English, Welsh, Scottish, Northern Irish or British") %>% 
  mutate(year = "All") %>% 
  select(year, prop, N) %>% 
  rename(n1 = N) %>% 
  rename(p1 = prop)
second_mmr_pc_rd <- second_mmr_pc_groups %>% 
  left_join(second_mmr_pc_ref, by = c("year"))
second_mmr_pc_year_groups <- mmr_second_pc_year %>% 
  filter(M_eth18_2011 != "All") %>% 
  filter(M_eth18_2011 != "English, Welsh, Scottish, Northern Irish or British") %>% 
  select(year, prop, N, M_eth18_2011) %>% 
  rename(n2 = N) %>% 
  rename(p2 = prop)
second_mmr_pc_year_ref <- mmr_second_pc_year %>% 
  filter(M_eth18_2011 == "English, Welsh, Scottish, Northern Irish or British") %>% 
  select(year, prop, N) %>% 
  rename(n1 = N) %>% 
  rename(p1 = prop)
second_mmr_pc_year_rd <- second_mmr_pc_year_groups %>% 
  left_join(second_mmr_pc_year_ref, by = c("year"))

# join all years to years
second_mmr_pc_year_rd <- rbind(second_mmr_pc_rd, second_mmr_pc_year_rd)

# calculate RDs and 95%CIs
second_mmr_pc_year_rd <- second_mmr_pc_year_rd %>% 
  mutate(rd = p1 - p2) %>% 
  mutate(lowerci = rd - 1.96*sqrt((p1*(1-p1))/n1) + ((p2*(1-p2))/n2)) %>% 
  mutate(upperci = rd + 1.96*sqrt((p1*(1-p1))/n1) + ((p2*(1-p2))/n2)) %>% 
  mutate(rdpercent = round(rd*100, 1)) %>% 
  mutate(lowercipercent = round(lowerci*100, 1)) %>% 
  mutate(uppercipercent = round(upperci*100, 1)) %>% 
  mutate(rd = round(rd, 2)) %>% 
  mutate(lowerci = round(lowerci, 2)) %>% 
  mutate(upperci = round(upperci, 2)) %>% 
  mutate(`RD (95%CI)` = paste0(rd," (",lowerci,",",upperci,")")) %>% 
  mutate(`%RD (95%CI)` = paste0(rdpercent," (",lowercipercent,",",uppercipercent,")"))  

# combine RD with % coverage for paper
second_mmr_pc_year_rd <- second_mmr_pc_year_rd %>% 
  mutate(coverage = round(p2*100, 1)) %>% 
  mutate(`% coverage (%RD)` = paste0(coverage," (",rdpercent,")"))

write_csv(second_mmr_pc_year_rd, file="~/second_mmr_pc_year_rd.csv") 

# Risk differences by year of vaccination - booster  ---------------------------------------------------------------------

# load prop by year
fifth_mmr_booster <- read_csv(file="~/fifth_mmr_booster.csv") 
mmr_fifth_booster_year <- read_csv(file="~/mmr_fifth_booster_year.csv") 

# restructure to remove 'all' and join ref to rest
fifth_mmr_booster_groups <- fifth_mmr_booster %>% 
  filter(M_eth18_2011 != "All") %>% 
  filter(M_eth18_2011 != "English, Welsh, Scottish, Northern Irish or British") %>% 
  mutate(year = "All") %>% 
  select(year, prop, N, M_eth18_2011) %>% 
  rename(n2 = N) %>% 
  rename(p2 = prop)
fifth_mmr_booster_ref <- fifth_mmr_booster %>% 
  filter(M_eth18_2011 == "English, Welsh, Scottish, Northern Irish or British") %>% 
  mutate(year = "All") %>% 
  select(year, prop, N) %>% 
  rename(n1 = N) %>% 
  rename(p1 = prop)
fifth_mmr_booster_rd <- fifth_mmr_booster_groups %>% 
  left_join(fifth_mmr_booster_ref, by = c("year"))
fifth_mmr_booster_year_groups <- mmr_fifth_booster_year %>% 
  filter(M_eth18_2011 != "All") %>% 
  filter(M_eth18_2011 != "English, Welsh, Scottish, Northern Irish or British") %>% 
  select(year, prop, N, M_eth18_2011) %>% 
  rename(n2 = N) %>% 
  rename(p2 = prop)
fifth_mmr_booster_year_ref <- mmr_fifth_booster_year %>% 
  filter(M_eth18_2011 == "English, Welsh, Scottish, Northern Irish or British") %>% 
  select(year, prop, N) %>% 
  rename(n1 = N) %>% 
  rename(p1 = prop)
fifth_mmr_booster_year_rd <- fifth_mmr_booster_year_groups %>% 
  left_join(fifth_mmr_booster_year_ref, by = c("year"))

# join all years to years
fifth_mmr_booster_year_rd <- rbind(fifth_mmr_booster_rd, fifth_mmr_booster_year_rd)

# calculate RDs and 95%CIs
fifth_mmr_booster_year_rd <- fifth_mmr_booster_year_rd %>% 
  mutate(rd = p1 - p2) %>% 
  mutate(lowerci = rd - 1.96*sqrt((p1*(1-p1))/n1) + ((p2*(1-p2))/n2)) %>% 
  mutate(upperci = rd + 1.96*sqrt((p1*(1-p1))/n1) + ((p2*(1-p2))/n2)) %>% 
  mutate(rdpercent = round(rd*100, 1)) %>% 
  mutate(lowercipercent = round(lowerci*100, 1)) %>% 
  mutate(uppercipercent = round(upperci*100, 1)) %>% 
  mutate(rd = round(rd, 2)) %>% 
  mutate(lowerci = round(lowerci, 2)) %>% 
  mutate(upperci = round(upperci, 2)) %>% 
  mutate(`RD (95%CI)` = paste0(rd," (",lowerci,",",upperci,")")) %>% 
  mutate(`%RD (95%CI)` = paste0(rdpercent," (",lowercipercent,",",uppercipercent,")"))  

# combine RD with % coverage for paper
fifth_mmr_booster_year_rd <- fifth_mmr_booster_year_rd %>% 
  mutate(coverage = round(p2*100, 1)) %>% 
  mutate(`% coverage (%RD)` = paste0(coverage," (",rdpercent,")"))

write_csv(fifth_mmr_booster_year_rd, file="~/fifth_mmr_booster_year_rd.csv") 

#--------------------------------------------------------------------------
# MMR - STATISTICAL ANALYSIS                   
#--------------------------------------------------------------------------

load(file="~/second_bday_mmr_analysis.Rdata")
load(file="~/fifth_bday_mmr_analysis.Rdata")

# Set reference groups ----------------------------------------------------

second_bday_mmr_analysis$matage <- relevel(second_bday_mmr_analysis$matage, ref = "30-34") # as per Li et al maternal checks paper
second_bday_mmr_analysis$year <- relevel(second_bday_mmr_analysis$year, ref = "2008-09")
fifth_bday_mmr_analysis$matage <- relevel(fifth_bday_mmr_analysis$matage, ref = "30-34")
fifth_bday_mmr_analysis$year <- relevel(fifth_bday_mmr_analysis$year, ref = "2011-12")

############################## Primary course ############################# --------------------------
# Unadjusted model + effect modification by time period ----------------------

# group years into 3 year intervals, except COVID year (2020-21)
second_bday_mmr_analysis <- second_bday_mmr_analysis %>% 
  mutate(yearcat = NA) %>% 
  mutate(yearcat = replace(yearcat, year == "2008-09" | year == "2009-10" | year == "2010-11", "2008-09 to 2010-11")) %>% 
  mutate(yearcat = replace(yearcat, year == "2011-12" | year == "2012-13" | year == "2013-14", "2011-12 to 2013-14")) %>% 
  mutate(yearcat = replace(yearcat, year == "2014-15" | year == "2015-16" | year == "2016-17", "2014-15 to 2016-17")) %>% 
  mutate(yearcat = replace(yearcat, year == "2017-18" | year == "2018-19" | year == "2019-20", "2017-18 to 2019-20")) %>% 
  mutate(yearcat = replace(yearcat, year == "2020-21", "2020-21")) %>% 
  mutate(yearcat = as.factor(yearcat)) %>% 
  mutate(yearcat = fct_relevel(yearcat, "2008-09 to 2010-11", 
                               "2011-12 to 2013-14",
                               "2014-15 to 2016-17", 
                               "2017-18 to 2019-20",
                               "2020-21"))

# run modified poisson models for likelihood ratio test
model_univar_mpoisson_year_mmr_pc <- glm(mmr_pc ~ M_eth18_2011 +
                                             yearcat,
                                           data = second_bday_mmr_analysis, family = poisson(link = "log"))
model_univar_mpoisson_yearEM_mmr_pc <- glm(mmr_pc ~ M_eth18_2011*yearcat, 
                                             data = second_bday_mmr_analysis, family = poisson(link = "log"))
lrtest(model_univar_mpoisson_year_mmr_pc, model_univar_mpoisson_yearEM_mmr_pc) 

# run modified poisson with effect modification by time period and robust sandwich estimator
univar_mpoisson_yearEM_mmr_pcsandwich <- coeftest(model_univar_mpoisson_yearEM_mmr_pc, vcov = sandwich) 

# extract results and save
univar_mpoisson_yearEM_mmr_pcsandwich_clean <- round(cbind(exp(cbind(RR = univar_mpoisson_yearEM_mmr_pcsandwich[,1],
                                                                       LCI = univar_mpoisson_yearEM_mmr_pcsandwich[,1] + qnorm(0.05/2)*univar_mpoisson_yearEM_mmr_pcsandwich[,2],
                                                                       UCI = univar_mpoisson_yearEM_mmr_pcsandwich[,1] - qnorm(0.05/2)*univar_mpoisson_yearEM_mmr_pcsandwich[,2])),
                                                             P = univar_mpoisson_yearEM_mmr_pcsandwich[,4]),2)
univar_mpoisson_yearEM_mmr_pcsandwich_clean <- as.data.frame(univar_mpoisson_yearEM_mmr_pcsandwich_clean)
extract_mpoisson_results(model_univar_mpoisson_yearEM_mmr_pc, univar_mpoisson_yearEM_mmr_pcsandwich_clean, 'univar_mpoisson_yearEM_mmr_pc', second_bday_mmr_analysis, M_eth18_2011)
write_csv(univar_mpoisson_yearEM_mmr_pc_table, file="~/univar_mpoisson_yearEM_mmr_pc_table.csv") 

# for plotting
write_csv(univar_mpoisson_yearEM_mmr_pc, file="~/univar_mpoisson_yearEM_mmr_pc.csv") 


# Calculate stratum specific effects from unadjusted model --------------------------------------

ethnicity <- c("Irish", "Any other White background",
               "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background",
               "Caribbean", "African", "Any other Black, African or Caribbean background",
               "White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed or multiple ethnic background",
               "Any other ethnic group", "Unknown")

## 2008-09 to 2010-2011
univar_mpoisson_yearEM_mmr_pc_stratum_2008 <- univar_mpoisson_yearEM_mmr_pcsandwich_clean[2:17,] %>% 
  cbind(ethnicity) %>% 
  select(-P) %>% 
  rename(estimate = RR) %>% 
  rename(upper = UCI) %>% 
  rename(lower = LCI) %>% 
  mutate(period = "2008-09 to 2010-11")

univar_mpoisson_yearEM_mmr_pc_stratum_2008_table <- univar_mpoisson_yearEM_mmr_pc_stratum_2008 %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, period)

## 2011-12 to 2013-14

# calculate stratum specific effects - original model
irish_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011African' = 1, 'M_eth18_2011African:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_2011 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:yearcat2011-12 to 2013-14")
aow_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:yearcat2011-12 to 2013-14")
indian_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:yearcat2011-12 to 2013-14")
pakistani_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:yearcat2011-12 to 2013-14")
bangladeshi_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:yearcat2011-12 to 2013-14")
chinese_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:yearcat2011-12 to 2013-14")
aoa_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:yearcat2011-12 to 2013-14")
caribbean_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:yearcat2011-12 to 2013-14")
african_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011African", "M_eth18_2011African:yearcat2011-12 to 2013-14")
aob_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:yearcat2011-12 to 2013-14")
wbc_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:yearcat2011-12 to 2013-14")
wba_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:yearcat2011-12 to 2013-14")
wa_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:yearcat2011-12 to 2013-14")
aom_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2011-12 to 2013-14")
aoe_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:yearcat2011-12 to 2013-14")
unknown_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:yearcat2011-12 to 2013-14")

all_stratum_estimates <- rbind(irish_vs_ref_2011, aow_vs_ref_2011,
                               indian_vs_ref_2011, pakistani_vs_ref_2011,
                               bangladeshi_vs_ref_2011, chinese_vs_ref_2011,
                               aoa_vs_ref_2011, caribbean_vs_ref_2011,
                               african_vs_ref_2011, aob_vs_ref_2011,
                               wbc_vs_ref_2011, wba_vs_ref_2011,
                               wa_vs_ref_2011, aom_vs_ref_2011,
                               aoe_vs_ref_2011, unknown_vs_ref_2011)
univar_mpoisson_yearEM_mmr_pc_stratum_2011 <- as.data.frame(t(cbind(data.frame(irish_vs_ref_2011_sandwich), data.frame(aow_vs_ref_2011_sandwich),
                                                                      data.frame(indian_vs_ref_2011_sandwich), data.frame(pakistani_vs_ref_2011_sandwich),
                                                                      data.frame(bangladeshi_vs_ref_2011_sandwich), data.frame(chinese_vs_ref_2011_sandwich),
                                                                      data.frame(aoa_vs_ref_2011_sandwich), data.frame(caribbean_vs_ref_2011_sandwich),
                                                                      data.frame(african_vs_ref_2011_sandwich), data.frame(aob_vs_ref_2011_sandwich),
                                                                      data.frame(wbc_vs_ref_2011_sandwich), data.frame(wba_vs_ref_2011_sandwich),
                                                                      data.frame(wa_vs_ref_2011_sandwich), data.frame(aom_vs_ref_2011_sandwich),
                                                                      data.frame(aoe_vs_ref_2011_sandwich), data.frame(unknown_vs_ref_2011_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(period = "2011-12 to 2013-14")

univar_mpoisson_yearEM_mmr_pc_stratum_2011_table <- univar_mpoisson_yearEM_mmr_pc_stratum_2011 %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, period)

## 2014-15 to 2016-17

# calculate stratum specific effects - original model
irish_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011African' = 1, 'M_eth18_2011African:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_2014 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:yearcat2014-15 to 2016-17")
aow_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:yearcat2014-15 to 2016-17")
indian_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:yearcat2014-15 to 2016-17")
pakistani_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:yearcat2014-15 to 2016-17")
bangladeshi_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:yearcat2014-15 to 2016-17")
chinese_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:yearcat2014-15 to 2016-17")
aoa_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:yearcat2014-15 to 2016-17")
caribbean_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:yearcat2014-15 to 2016-17")
african_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011African", "M_eth18_2011African:yearcat2014-15 to 2016-17")
aob_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:yearcat2014-15 to 2016-17")
wbc_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:yearcat2014-15 to 2016-17")
wba_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:yearcat2014-15 to 2016-17")
wa_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:yearcat2014-15 to 2016-17")
aom_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2014-15 to 2016-17")
aoe_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:yearcat2014-15 to 2016-17")
unknown_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:yearcat2014-15 to 2016-17")

all_stratum_estimates <- rbind(irish_vs_ref_2014, aow_vs_ref_2014,
                               indian_vs_ref_2014, pakistani_vs_ref_2014,
                               bangladeshi_vs_ref_2014, chinese_vs_ref_2014,
                               aoa_vs_ref_2014, caribbean_vs_ref_2014,
                               african_vs_ref_2014, aob_vs_ref_2014,
                               wbc_vs_ref_2014, wba_vs_ref_2014,
                               wa_vs_ref_2014, aom_vs_ref_2014,
                               aoe_vs_ref_2014, unknown_vs_ref_2014)
univar_mpoisson_yearEM_mmr_pc_stratum_2014 <- as.data.frame(t(cbind(data.frame(irish_vs_ref_2014_sandwich), data.frame(aow_vs_ref_2014_sandwich),
                                                                      data.frame(indian_vs_ref_2014_sandwich), data.frame(pakistani_vs_ref_2014_sandwich),
                                                                      data.frame(bangladeshi_vs_ref_2014_sandwich), data.frame(chinese_vs_ref_2014_sandwich),
                                                                      data.frame(aoa_vs_ref_2014_sandwich), data.frame(caribbean_vs_ref_2014_sandwich),
                                                                      data.frame(african_vs_ref_2014_sandwich), data.frame(aob_vs_ref_2014_sandwich),
                                                                      data.frame(wbc_vs_ref_2014_sandwich), data.frame(wba_vs_ref_2014_sandwich),
                                                                      data.frame(wa_vs_ref_2014_sandwich), data.frame(aom_vs_ref_2014_sandwich),
                                                                      data.frame(aoe_vs_ref_2014_sandwich), data.frame(unknown_vs_ref_2014_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(period = "2014-15 to 2016-17")

univar_mpoisson_yearEM_mmr_pc_stratum_2014_table <- univar_mpoisson_yearEM_mmr_pc_stratum_2014 %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, period)

## 2017-18 to 2019-20

# calculate stratum specific effects - original model
irish_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011African' = 1, 'M_eth18_2011African:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_2017 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:yearcat2017-18 to 2019-20")
aow_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:yearcat2017-18 to 2019-20")
indian_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:yearcat2017-18 to 2019-20")
pakistani_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:yearcat2017-18 to 2019-20")
bangladeshi_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:yearcat2017-18 to 2019-20")
chinese_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:yearcat2017-18 to 2019-20")
aoa_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:yearcat2017-18 to 2019-20")
caribbean_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:yearcat2017-18 to 2019-20")
african_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011African", "M_eth18_2011African:yearcat2017-18 to 2019-20")
aob_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:yearcat2017-18 to 2019-20")
wbc_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:yearcat2017-18 to 2019-20")
wba_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:yearcat2017-18 to 2019-20")
wa_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:yearcat2017-18 to 2019-20")
aom_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2017-18 to 2019-20")
aoe_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:yearcat2017-18 to 2019-20")
unknown_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:yearcat2017-18 to 2019-20")

all_stratum_estimates <- rbind(irish_vs_ref_2017, aow_vs_ref_2017,
                               indian_vs_ref_2017, pakistani_vs_ref_2017,
                               bangladeshi_vs_ref_2017, chinese_vs_ref_2017,
                               aoa_vs_ref_2017, caribbean_vs_ref_2017,
                               african_vs_ref_2017, aob_vs_ref_2017,
                               wbc_vs_ref_2017, wba_vs_ref_2017,
                               wa_vs_ref_2017, aom_vs_ref_2017,
                               aoe_vs_ref_2017, unknown_vs_ref_2017)
univar_mpoisson_yearEM_mmr_pc_stratum_2017 <- as.data.frame(t(cbind(data.frame(irish_vs_ref_2017_sandwich), data.frame(aow_vs_ref_2017_sandwich),
                                                                      data.frame(indian_vs_ref_2017_sandwich), data.frame(pakistani_vs_ref_2017_sandwich),
                                                                      data.frame(bangladeshi_vs_ref_2017_sandwich), data.frame(chinese_vs_ref_2017_sandwich),
                                                                      data.frame(aoa_vs_ref_2017_sandwich), data.frame(caribbean_vs_ref_2017_sandwich),
                                                                      data.frame(african_vs_ref_2017_sandwich), data.frame(aob_vs_ref_2017_sandwich),
                                                                      data.frame(wbc_vs_ref_2017_sandwich), data.frame(wba_vs_ref_2017_sandwich),
                                                                      data.frame(wa_vs_ref_2017_sandwich), data.frame(aom_vs_ref_2017_sandwich),
                                                                      data.frame(aoe_vs_ref_2017_sandwich), data.frame(unknown_vs_ref_2017_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(period = "2017-18 to 2019-20")

univar_mpoisson_yearEM_mmr_pc_stratum_2017_table <- univar_mpoisson_yearEM_mmr_pc_stratum_2017 %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, period)

## 2020-21

# calculate stratum specific effects - original model
irish_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011African' = 1, 'M_eth18_2011African:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_2020 <- as.data.frame(exp(estimable(model_univar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:yearcat2020-21")
aow_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:yearcat2020-21")
indian_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:yearcat2020-21")
pakistani_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:yearcat2020-21")
bangladeshi_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:yearcat2020-21")
chinese_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:yearcat2020-21")
aoa_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:yearcat2020-21")
caribbean_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:yearcat2020-21")
african_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011African", "M_eth18_2011African:yearcat2020-21")
aob_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:yearcat2020-21")
wbc_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:yearcat2020-21")
wba_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:yearcat2020-21")
wa_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:yearcat2020-21")
aom_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2020-21")
aoe_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:yearcat2020-21")
unknown_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(univar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:yearcat2020-21")

all_stratum_estimates <- rbind(irish_vs_ref_2020, aow_vs_ref_2020,
                               indian_vs_ref_2020, pakistani_vs_ref_2020,
                               bangladeshi_vs_ref_2020, chinese_vs_ref_2020,
                               aoa_vs_ref_2020, caribbean_vs_ref_2020,
                               african_vs_ref_2020, aob_vs_ref_2020,
                               wbc_vs_ref_2020, wba_vs_ref_2020,
                               wa_vs_ref_2020, aom_vs_ref_2020,
                               aoe_vs_ref_2020, unknown_vs_ref_2020)
univar_mpoisson_yearEM_mmr_pc_stratum_2020 <- as.data.frame(t(cbind(data.frame(irish_vs_ref_2020_sandwich), data.frame(aow_vs_ref_2020_sandwich),
                                                                      data.frame(indian_vs_ref_2020_sandwich), data.frame(pakistani_vs_ref_2020_sandwich),
                                                                      data.frame(bangladeshi_vs_ref_2020_sandwich), data.frame(chinese_vs_ref_2020_sandwich),
                                                                      data.frame(aoa_vs_ref_2020_sandwich), data.frame(caribbean_vs_ref_2020_sandwich),
                                                                      data.frame(african_vs_ref_2020_sandwich), data.frame(aob_vs_ref_2020_sandwich),
                                                                      data.frame(wbc_vs_ref_2020_sandwich), data.frame(wba_vs_ref_2020_sandwich),
                                                                      data.frame(wa_vs_ref_2020_sandwich), data.frame(aom_vs_ref_2020_sandwich),
                                                                      data.frame(aoe_vs_ref_2020_sandwich), data.frame(unknown_vs_ref_2020_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(period = "2020-21")

univar_mpoisson_yearEM_mmr_pc_stratum_2020_table <- univar_mpoisson_yearEM_mmr_pc_stratum_2020 %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, period)

# combine years
univar_mpoisson_yearEM_mmr_pc_strata <- rbind(univar_mpoisson_yearEM_mmr_pc_stratum_2008,
                                                univar_mpoisson_yearEM_mmr_pc_stratum_2011,
                                                univar_mpoisson_yearEM_mmr_pc_stratum_2014,
                                                univar_mpoisson_yearEM_mmr_pc_stratum_2017,
                                                univar_mpoisson_yearEM_mmr_pc_stratum_2020)

univar_mpoisson_yearEM_mmr_pc_strata_table <-rbind(univar_mpoisson_yearEM_mmr_pc_stratum_2008_table,
                                                     univar_mpoisson_yearEM_mmr_pc_stratum_2011_table,
                                                     univar_mpoisson_yearEM_mmr_pc_stratum_2014_table,
                                                     univar_mpoisson_yearEM_mmr_pc_stratum_2017_table,
                                                     univar_mpoisson_yearEM_mmr_pc_stratum_2020_table)

# save
write_csv(univar_mpoisson_yearEM_mmr_pc_strata, file="~/univar_mpoisson_yearEM_mmr_pc_strata.csv") 
write_csv(univar_mpoisson_yearEM_mmr_pc_strata_table, file="~/univar_mpoisson_yearEM_mmr_pc_strata_table.csv") 

# remove model
rm(model_univar_mpoisson_year_mmr_pc, model_univar_mpoisson_yearEM_mmr_pc)


# Adjusted model (sociodemographic + maternal/birth factors) + effect modification by time period ----------------------

# group years into 3 year intervals, except COVID year (2020-21)
second_bday_mmr_analysis <- second_bday_mmr_analysis %>% 
  mutate(yearcat = NA) %>% 
  mutate(yearcat = replace(yearcat, year == "2008-09" | year == "2009-10" | year == "2010-11", "2008-09 to 2010-11")) %>% 
  mutate(yearcat = replace(yearcat, year == "2011-12" | year == "2012-13" | year == "2013-14", "2011-12 to 2013-14")) %>% 
  mutate(yearcat = replace(yearcat, year == "2014-15" | year == "2015-16" | year == "2016-17", "2014-15 to 2016-17")) %>% 
  mutate(yearcat = replace(yearcat, year == "2017-18" | year == "2018-19" | year == "2019-20", "2017-18 to 2019-20")) %>% 
  mutate(yearcat = replace(yearcat, year == "2020-21", "2020-21")) %>% 
  mutate(yearcat = as.factor(yearcat)) %>% 
  mutate(yearcat = fct_relevel(yearcat, "2008-09 to 2010-11", 
                               "2011-12 to 2013-14",
                               "2014-15 to 2016-17", 
                               "2017-18 to 2019-20",
                               "2020-21"))

# run modified poisson models for likelihood ratio test
model_multivar_mpoisson_year_mmr_pc <- glm(mmr_pc ~ M_eth18_2011 +
                                             yearcat + 
                                             region + imd + ruc +  
                                             preterm + modebirth + matage + first_time_mum, 
                                           data = second_bday_mmr_analysis, family = poisson(link = "log"))
model_multivar_mpoisson_yearEM_mmr_pc <- glm(mmr_pc ~ M_eth18_2011*yearcat + 
                                               region + imd + ruc +  
                                               preterm + modebirth + matage + first_time_mum, 
                                             data = second_bday_mmr_analysis, family = poisson(link = "log"))
lrtest(model_multivar_mpoisson_year_mmr_pc, model_multivar_mpoisson_yearEM_mmr_pc) 

# run modified poisson with effect modification by time period and robust sandwich estimator
multivar_mpoisson_yearEM_mmr_pcsandwich <- coeftest(model_multivar_mpoisson_yearEM_mmr_pc, vcov = sandwich) 

# extract results and save
multivar_mpoisson_yearEM_mmr_pcsandwich_clean <- round(cbind(exp(cbind(RR = multivar_mpoisson_yearEM_mmr_pcsandwich[,1],
                                                                       LCI = multivar_mpoisson_yearEM_mmr_pcsandwich[,1] + qnorm(0.05/2)*multivar_mpoisson_yearEM_mmr_pcsandwich[,2],
                                                                       UCI = multivar_mpoisson_yearEM_mmr_pcsandwich[,1] - qnorm(0.05/2)*multivar_mpoisson_yearEM_mmr_pcsandwich[,2])),
                                                             P = multivar_mpoisson_yearEM_mmr_pcsandwich[,4]),2)
multivar_mpoisson_yearEM_mmr_pcsandwich_clean <- as.data.frame(multivar_mpoisson_yearEM_mmr_pcsandwich_clean)
extract_mpoisson_results(model_multivar_mpoisson_yearEM_mmr_pc, multivar_mpoisson_yearEM_mmr_pcsandwich_clean, 'multivar_mpoisson_yearEM_mmr_pc', second_bday_mmr_analysis, M_eth18_2011)
write_csv(multivar_mpoisson_yearEM_mmr_pc_table, file="~/multivar_mpoisson_yearEM_mmr_pc_table.csv") 

# for plotting
write_csv(multivar_mpoisson_yearEM_mmr_pc, file="~/multivar_mpoisson_yearEM_mmr_pc.csv") 

# Calculate stratum specific effects from adjusted model --------------------------------------

ethnicity <- c("Irish", "Any other White background",
                  "Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background",
                  "Caribbean", "African", "Any other Black, African or Caribbean background",
                  "White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed or multiple ethnic background",
                  "Any other ethnic group", "Unknown")

## 2008-09 to 2010-2011
multivar_mpoisson_yearEM_mmr_pc_stratum_2008 <- multivar_mpoisson_yearEM_mmr_pcsandwich_clean[2:17,] %>% 
  cbind(ethnicity) %>% 
  select(-P) %>% 
  rename(estimate = RR) %>% 
  rename(upper = UCI) %>% 
  rename(lower = LCI) %>% 
  mutate(period = "2008-09 to 2010-11")

multivar_mpoisson_yearEM_mmr_pc_stratum_2008_table <- multivar_mpoisson_yearEM_mmr_pc_stratum_2008 %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, period)

## 2011-12 to 2013-14

# calculate stratum specific effects - original model
irish_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011African' = 1, 'M_eth18_2011African:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_2011 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:yearcat2011-12 to 2013-14' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:yearcat2011-12 to 2013-14")
aow_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:yearcat2011-12 to 2013-14")
indian_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:yearcat2011-12 to 2013-14")
pakistani_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:yearcat2011-12 to 2013-14")
bangladeshi_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:yearcat2011-12 to 2013-14")
chinese_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:yearcat2011-12 to 2013-14")
aoa_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:yearcat2011-12 to 2013-14")
caribbean_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:yearcat2011-12 to 2013-14")
african_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011African", "M_eth18_2011African:yearcat2011-12 to 2013-14")
aob_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:yearcat2011-12 to 2013-14")
wbc_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:yearcat2011-12 to 2013-14")
wba_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:yearcat2011-12 to 2013-14")
wa_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:yearcat2011-12 to 2013-14")
aom_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2011-12 to 2013-14")
aoe_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:yearcat2011-12 to 2013-14")
unknown_vs_ref_2011_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:yearcat2011-12 to 2013-14")

all_stratum_estimates <- rbind(irish_vs_ref_2011, aow_vs_ref_2011,
                                                indian_vs_ref_2011, pakistani_vs_ref_2011,
                                                bangladeshi_vs_ref_2011, chinese_vs_ref_2011,
                                                aoa_vs_ref_2011, caribbean_vs_ref_2011,
                                                african_vs_ref_2011, aob_vs_ref_2011,
                                                wbc_vs_ref_2011, wba_vs_ref_2011,
                                                wa_vs_ref_2011, aom_vs_ref_2011,
                                                aoe_vs_ref_2011, unknown_vs_ref_2011)
multivar_mpoisson_yearEM_mmr_pc_stratum_2011 <- as.data.frame(t(cbind(data.frame(irish_vs_ref_2011_sandwich), data.frame(aow_vs_ref_2011_sandwich),
                                        data.frame(indian_vs_ref_2011_sandwich), data.frame(pakistani_vs_ref_2011_sandwich),
                                        data.frame(bangladeshi_vs_ref_2011_sandwich), data.frame(chinese_vs_ref_2011_sandwich),
                                        data.frame(aoa_vs_ref_2011_sandwich), data.frame(caribbean_vs_ref_2011_sandwich),
                                        data.frame(african_vs_ref_2011_sandwich), data.frame(aob_vs_ref_2011_sandwich),
                                        data.frame(wbc_vs_ref_2011_sandwich), data.frame(wba_vs_ref_2011_sandwich),
                                        data.frame(wa_vs_ref_2011_sandwich), data.frame(aom_vs_ref_2011_sandwich),
                                        data.frame(aoe_vs_ref_2011_sandwich), data.frame(unknown_vs_ref_2011_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(period = "2011-12 to 2013-14")

multivar_mpoisson_yearEM_mmr_pc_stratum_2011_table <- multivar_mpoisson_yearEM_mmr_pc_stratum_2011 %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, period)

## 2014-15 to 2016-17

# calculate stratum specific effects - original model
irish_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011African' = 1, 'M_eth18_2011African:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_2014 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:yearcat2014-15 to 2016-17' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:yearcat2014-15 to 2016-17")
aow_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:yearcat2014-15 to 2016-17")
indian_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:yearcat2014-15 to 2016-17")
pakistani_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:yearcat2014-15 to 2016-17")
bangladeshi_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:yearcat2014-15 to 2016-17")
chinese_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:yearcat2014-15 to 2016-17")
aoa_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:yearcat2014-15 to 2016-17")
caribbean_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:yearcat2014-15 to 2016-17")
african_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011African", "M_eth18_2011African:yearcat2014-15 to 2016-17")
aob_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:yearcat2014-15 to 2016-17")
wbc_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:yearcat2014-15 to 2016-17")
wba_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:yearcat2014-15 to 2016-17")
wa_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:yearcat2014-15 to 2016-17")
aom_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2014-15 to 2016-17")
aoe_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:yearcat2014-15 to 2016-17")
unknown_vs_ref_2014_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:yearcat2014-15 to 2016-17")

all_stratum_estimates <- rbind(irish_vs_ref_2014, aow_vs_ref_2014,
                               indian_vs_ref_2014, pakistani_vs_ref_2014,
                               bangladeshi_vs_ref_2014, chinese_vs_ref_2014,
                               aoa_vs_ref_2014, caribbean_vs_ref_2014,
                               african_vs_ref_2014, aob_vs_ref_2014,
                               wbc_vs_ref_2014, wba_vs_ref_2014,
                               wa_vs_ref_2014, aom_vs_ref_2014,
                               aoe_vs_ref_2014, unknown_vs_ref_2014)
multivar_mpoisson_yearEM_mmr_pc_stratum_2014 <- as.data.frame(t(cbind(data.frame(irish_vs_ref_2014_sandwich), data.frame(aow_vs_ref_2014_sandwich),
                                                                      data.frame(indian_vs_ref_2014_sandwich), data.frame(pakistani_vs_ref_2014_sandwich),
                                                                      data.frame(bangladeshi_vs_ref_2014_sandwich), data.frame(chinese_vs_ref_2014_sandwich),
                                                                      data.frame(aoa_vs_ref_2014_sandwich), data.frame(caribbean_vs_ref_2014_sandwich),
                                                                      data.frame(african_vs_ref_2014_sandwich), data.frame(aob_vs_ref_2014_sandwich),
                                                                      data.frame(wbc_vs_ref_2014_sandwich), data.frame(wba_vs_ref_2014_sandwich),
                                                                      data.frame(wa_vs_ref_2014_sandwich), data.frame(aom_vs_ref_2014_sandwich),
                                                                      data.frame(aoe_vs_ref_2014_sandwich), data.frame(unknown_vs_ref_2014_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(period = "2014-15 to 2016-17")

multivar_mpoisson_yearEM_mmr_pc_stratum_2014_table <- multivar_mpoisson_yearEM_mmr_pc_stratum_2014 %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, period)

## 2017-18 to 2019-20

# calculate stratum specific effects - original model
irish_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011African' = 1, 'M_eth18_2011African:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_2017 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:yearcat2017-18 to 2019-20' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:yearcat2017-18 to 2019-20")
aow_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:yearcat2017-18 to 2019-20")
indian_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:yearcat2017-18 to 2019-20")
pakistani_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:yearcat2017-18 to 2019-20")
bangladeshi_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:yearcat2017-18 to 2019-20")
chinese_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:yearcat2017-18 to 2019-20")
aoa_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:yearcat2017-18 to 2019-20")
caribbean_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:yearcat2017-18 to 2019-20")
african_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011African", "M_eth18_2011African:yearcat2017-18 to 2019-20")
aob_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:yearcat2017-18 to 2019-20")
wbc_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:yearcat2017-18 to 2019-20")
wba_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:yearcat2017-18 to 2019-20")
wa_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:yearcat2017-18 to 2019-20")
aom_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2017-18 to 2019-20")
aoe_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:yearcat2017-18 to 2019-20")
unknown_vs_ref_2017_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:yearcat2017-18 to 2019-20")

all_stratum_estimates <- rbind(irish_vs_ref_2017, aow_vs_ref_2017,
                               indian_vs_ref_2017, pakistani_vs_ref_2017,
                               bangladeshi_vs_ref_2017, chinese_vs_ref_2017,
                               aoa_vs_ref_2017, caribbean_vs_ref_2017,
                               african_vs_ref_2017, aob_vs_ref_2017,
                               wbc_vs_ref_2017, wba_vs_ref_2017,
                               wa_vs_ref_2017, aom_vs_ref_2017,
                               aoe_vs_ref_2017, unknown_vs_ref_2017)
multivar_mpoisson_yearEM_mmr_pc_stratum_2017 <- as.data.frame(t(cbind(data.frame(irish_vs_ref_2017_sandwich), data.frame(aow_vs_ref_2017_sandwich),
                                                                      data.frame(indian_vs_ref_2017_sandwich), data.frame(pakistani_vs_ref_2017_sandwich),
                                                                      data.frame(bangladeshi_vs_ref_2017_sandwich), data.frame(chinese_vs_ref_2017_sandwich),
                                                                      data.frame(aoa_vs_ref_2017_sandwich), data.frame(caribbean_vs_ref_2017_sandwich),
                                                                      data.frame(african_vs_ref_2017_sandwich), data.frame(aob_vs_ref_2017_sandwich),
                                                                      data.frame(wbc_vs_ref_2017_sandwich), data.frame(wba_vs_ref_2017_sandwich),
                                                                      data.frame(wa_vs_ref_2017_sandwich), data.frame(aom_vs_ref_2017_sandwich),
                                                                      data.frame(aoe_vs_ref_2017_sandwich), data.frame(unknown_vs_ref_2017_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(period = "2017-18 to 2019-20")

multivar_mpoisson_yearEM_mmr_pc_stratum_2017_table <- multivar_mpoisson_yearEM_mmr_pc_stratum_2017 %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, period)

## 2020-21

# calculate stratum specific effects - original model
irish_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Irish' = 1, 'M_eth18_2011Irish:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
aow_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other White background' = 1, 'M_eth18_2011Any other White background:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
indian_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Indian' = 1, 'M_eth18_2011Indian:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
pakistani_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Pakistani' = 1, 'M_eth18_2011Pakistani:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
bangladeshi_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Bangladeshi' = 1, 'M_eth18_2011Bangladeshi:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
chinese_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Chinese' = 1, 'M_eth18_2011Chinese:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
aoa_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Asian background' = 1, 'M_eth18_2011Any other Asian background:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
caribbean_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Caribbean' = 1, 'M_eth18_2011Caribbean:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
african_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011African' = 1, 'M_eth18_2011African:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
aob_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Black, African or Caribbean background' = 1, 'M_eth18_2011Any other Black, African or Caribbean background:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
wbc_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black Caribbean' = 1, 'M_eth18_2011White and Black Caribbean:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
wba_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Black African' = 1, 'M_eth18_2011White and Black African:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
wa_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011White and Asian' = 1, 'M_eth18_2011White and Asian:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
aom_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other Mixed or multiple ethnic background' = 1, 'M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
aoe_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Any other ethnic group' = 1, 'M_eth18_2011Any other ethnic group:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)
unknown_vs_ref_2020 <- as.data.frame(exp(estimable(model_multivar_mpoisson_yearEM_mmr_pc, c('M_eth18_2011Unknown' = 1, 'M_eth18_2011Unknown:yearcat2020-21' = 1), conf=.95))) %>% select(Estimate)

# calculate stratum specific effects - sandwich
irish_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Irish", "M_eth18_2011Irish:yearcat2020-21")
aow_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other White background", "M_eth18_2011Any other White background:yearcat2020-21")
indian_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Indian", "M_eth18_2011Indian:yearcat2020-21")
pakistani_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Pakistani", "M_eth18_2011Pakistani:yearcat2020-21")
bangladeshi_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Bangladeshi", "M_eth18_2011Bangladeshi:yearcat2020-21")
chinese_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Chinese", "M_eth18_2011Chinese:yearcat2020-21")
aoa_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Asian background", "M_eth18_2011Any other Asian background:yearcat2020-21")
caribbean_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Caribbean", "M_eth18_2011Caribbean:yearcat2020-21")
african_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011African", "M_eth18_2011African:yearcat2020-21")
aob_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Black, African or Caribbean background", "M_eth18_2011Any other Black, African or Caribbean background:yearcat2020-21")
wbc_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black Caribbean", "M_eth18_2011White and Black Caribbean:yearcat2020-21")
wba_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Black African", "M_eth18_2011White and Black African:yearcat2020-21")
wa_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011White and Asian", "M_eth18_2011White and Asian:yearcat2020-21")
aom_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other Mixed or multiple ethnic background", "M_eth18_2011Any other Mixed or multiple ethnic background:yearcat2020-21")
aoe_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Any other ethnic group", "M_eth18_2011Any other ethnic group:yearcat2020-21")
unknown_vs_ref_2020_sandwich <- get_sandwich_stratum_ci(multivar_mpoisson_yearEM_mmr_pcsandwich, "M_eth18_2011Unknown", "M_eth18_2011Unknown:yearcat2020-21")

all_stratum_estimates <- rbind(irish_vs_ref_2020, aow_vs_ref_2020,
                               indian_vs_ref_2020, pakistani_vs_ref_2020,
                               bangladeshi_vs_ref_2020, chinese_vs_ref_2020,
                               aoa_vs_ref_2020, caribbean_vs_ref_2020,
                               african_vs_ref_2020, aob_vs_ref_2020,
                               wbc_vs_ref_2020, wba_vs_ref_2020,
                               wa_vs_ref_2020, aom_vs_ref_2020,
                               aoe_vs_ref_2020, unknown_vs_ref_2020)
multivar_mpoisson_yearEM_mmr_pc_stratum_2020 <- as.data.frame(t(cbind(data.frame(irish_vs_ref_2020_sandwich), data.frame(aow_vs_ref_2020_sandwich),
                                                                      data.frame(indian_vs_ref_2020_sandwich), data.frame(pakistani_vs_ref_2020_sandwich),
                                                                      data.frame(bangladeshi_vs_ref_2020_sandwich), data.frame(chinese_vs_ref_2020_sandwich),
                                                                      data.frame(aoa_vs_ref_2020_sandwich), data.frame(caribbean_vs_ref_2020_sandwich),
                                                                      data.frame(african_vs_ref_2020_sandwich), data.frame(aob_vs_ref_2020_sandwich),
                                                                      data.frame(wbc_vs_ref_2020_sandwich), data.frame(wba_vs_ref_2020_sandwich),
                                                                      data.frame(wa_vs_ref_2020_sandwich), data.frame(aom_vs_ref_2020_sandwich),
                                                                      data.frame(aoe_vs_ref_2020_sandwich), data.frame(unknown_vs_ref_2020_sandwich)))) %>%
  cbind(ethnicity, all_stratum_estimates) %>%
  rename(lower = `2.5 %`) %>%
  rename(upper = `97.5 %`) %>%
  rename(estimate = Estimate) %>%
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(period = "2020-21")

multivar_mpoisson_yearEM_mmr_pc_stratum_2020_table <- multivar_mpoisson_yearEM_mmr_pc_stratum_2020 %>%
  mutate(rr_ci = paste0(estimate," (",lower,"-",upper, ")")) %>%
  select(ethnicity, rr_ci, period)

# combine years
multivar_mpoisson_yearEM_mmr_pc_strata <- rbind(multivar_mpoisson_yearEM_mmr_pc_stratum_2008,
                                                multivar_mpoisson_yearEM_mmr_pc_stratum_2011,
                                                multivar_mpoisson_yearEM_mmr_pc_stratum_2014,
                                                multivar_mpoisson_yearEM_mmr_pc_stratum_2017,
                                                multivar_mpoisson_yearEM_mmr_pc_stratum_2020)

multivar_mpoisson_yearEM_mmr_pc_strata_table <-rbind(multivar_mpoisson_yearEM_mmr_pc_stratum_2008_table,
                                                     multivar_mpoisson_yearEM_mmr_pc_stratum_2011_table,
                                                     multivar_mpoisson_yearEM_mmr_pc_stratum_2014_table,
                                                     multivar_mpoisson_yearEM_mmr_pc_stratum_2017_table,
                                                     multivar_mpoisson_yearEM_mmr_pc_stratum_2020_table)

# save
write_csv(multivar_mpoisson_yearEM_mmr_pc_strata, file="~/multivar_mpoisson_yearEM_mmr_pc_strata.csv") 
write_csv(multivar_mpoisson_yearEM_mmr_pc_strata_table, file="~/multivar_mpoisson_yearEM_mmr_pc_strata_table.csv") 

# remove model
rm(model_multivar_mpoisson_year_mmr_pc, model_multivar_mpoisson_yearEM_mmr_pc)

# Booster - repeat modelling with fifth birthday cohort  --------------------------

# Forestplot - unadjusted models -----------------------------

# load 
mmr_strata_pc <- read_csv(file="~/univar_mpoisson_yearEM_mmr_pc_strata.csv") 
mmr_strata_booster <- read_csv(file="~/univar_mpoisson_yearEM_mmr_booster_strata.csv") 

# plot pc
dev.new()
png(file = "~/univar_fp_mmr_pc_EM.png", width = 1500, height = 2000)
fp_mmr_pc_EM <- mmr_strata_pc %>% 
  group_by(period) %>%
  forestplot(label = ethnicity,
             mean = estimate,
             lower = lower,
             upper = upper,
             fn.ci_norm = c(fpDrawNormalCI,  fpDrawCircleCI, fpDrawPointCI,fpDrawDiamondCI, fpDrawNormalCI),
             boxsize = .125, # We set the box size to better visualize the type
             line.margin = .125, # We need to add this to avoid crowding
             xlog = TRUE,
             clip = c(0.1, 1.2),
             xticks = log10(c(0.31, 0.44, 0.60, 0.79, 1, 1.24)), 
             zero = 1,
             xlab = "RR (95%CI)",
             mar = unit(rep(10, times = 4), "mm"),
             col = fpColors(box = c("grey", "black", "black", "black", "black"),
                            lines = c("grey", "black", "black", "black", "black"),
                            zero = "black"),
             txt_gp = fpTxtGp(ticks=gpar(cex=2.5), xlab=gpar(cex=2.5), cex = 2.5, summary = gpar(fontface = 'bold'),
                              label = list(gpar(fontface = 'plain', cex = 2.5),
                                           gpar(fontface = 'bold', cex = 2.5),
                                           gpar(fontface = 'bold', cex = 2.5)))) %>% 
  fp_set_zebra_style("#F8F8F8") %>% 
  fp_add_lines(h_3 = gpar(lwd = 1.5, col = "#97928C"), 
               h_8 = gpar(lwd = 1.5, col = "#97928C"),
               h_11 = gpar(lwd = 1.5, col = "#97928C"),
               h_15 = gpar(lwd = 1.5, col = "#97928C"),
               h_16 = gpar(lwd = 1.5, col = "#97928C")) %>% 
  fp_decorate_graph(grid = structure(c(0.7, 0.801, 0.903), 
                                     gp = gpar(lty = 2, col = "#D3D3D3"))) 
print(fp_mmr_pc_EM)
dev.off()

# plot booster
dev.new()
png(file = "~/univar_fp_mmr_booster_EM.png", width = 1500, height = 2000)
fp_mmr_booster_EM <- mmr_strata_booster %>% 
  group_by(period) %>%
  forestplot(label = ethnicity,
             mean = estimate,
             lower = lower,
             upper = upper,
             fn.ci_norm = c(fpDrawCircleCI, fpDrawPointCI,fpDrawDiamondCI, fpDrawNormalCI),
             boxsize = .125, # We set the box size to better visualize the type
             line.margin = .175, # We need to add this to avoid crowding
             xlog = TRUE,
             clip = c(0.1, 1.2), 
             xticks = log10(c(0.31, 0.44, 0.6, 0.79, 1.0, 1.24)), 
             zero = 1,
             xlab = "RR (95%CI)",
             mar = unit(rep(10, times = 4), "mm"),
             col = fpColors(box = c("black", "black", "black", "black"),
                            lines = c("black", "black", "black", "black"),
                            zero = "black"),
             txt_gp = fpTxtGp(ticks=gpar(cex=2.5), xlab=gpar(cex=2.5), cex = 2.5, summary = gpar(fontface = 'bold'),
                              label = list(gpar(fontface = 'plain', cex = 2.5),
                                           gpar(fontface = 'bold', cex = 2.5),
                                           gpar(fontface = 'bold', cex = 2.5)))) %>% 
  fp_set_zebra_style("#F8F8F8") %>% 
  fp_add_lines(h_3 = gpar(lwd = 1.5, col = "#97928C"), 
               h_8 = gpar(lwd = 1.5, col = "#97928C"),
               h_11 = gpar(lwd = 1.5, col = "#97928C"),
               h_15 = gpar(lwd = 1.5, col = "#97928C"),
               h_16 = gpar(lwd = 1.5, col = "#97928C")) %>% 
  fp_decorate_graph(grid = structure(c(0.7, 0.801, 0.903), 
                                     gp = gpar(lty = 2, col = "#97928C")))
print(fp_mmr_booster_EM)
dev.off()


# Repeat for adjusted models -----------------------------
