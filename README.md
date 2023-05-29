# CPRD_Aurum_ethnicity_childhood_vaccination
Ethnicity and routine childhood vaccinations in England (CPRD Aurum)

This repository contains the code lists and code used to derive ethnicity and routine childhood vaccinations in England using the Clinical Practice Research Datalink (CPRD) Aurum database. 

Ethnicity files include:
- Ethnicity 2001 Census England & Wales (16 categories)
- Ethnicity 2011 Census England & Wales (18 categories)
- Ethnicity 2021 Census England & Wales (19 categories)
- R code for assigning each individual a single most plausible ethnic category using CPRD Aurum linked to HES Admitted Patient Care data

Routine childhood vaccination files include:
- MMR
- 6/5/4-in-1
- MenC & HibMenC
- MenB
- Rotavirus
- Pneumococcal
- R code for identifying vaccinations in the CPRD Aurum Observations and Drug Issue data tables, and deriving variables indicating whether the primary course and full course (primary + booster dose) had been completed by the expected birthday (first, second and fifth birthdays, as per the UKHSA/NHS Digital annual statistical reports)

The algorithms corresponding to these code lists and R code have been published in the supplemental appendix of XXX. Validation of completeness and representativeness are also contained within this paper. 
