#HPV004 FREQUENT VS. INFREQUENT SCREENING STUDY 


#INDEX
## 0 Packages and functions
## 1 Data loading and merging
## 2 Cytology and HPV DNA findings per visit
## 3 Prevalence and prevalence ratios
## 4 Hazards
## 5 Questionaire data
## 6 Leep exclusions

# 0 PACKAGES AND FUNCTIONS ----
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readr)
library(data.table)
require(tible)
library(tidyr)
library(stringr)
library(purrr)
library(plyr)
library(rlist)
library(gridExtra)
library(caret)
library(systemfonts)
library(gtable)
library(grid)
library(Epi)
library(lubridate)
library(survival)
library(survminer)
library(ggpubr)
require(renv)

##FUNCTION CIdel: Calculates the Confidence interval (Upper and lower bounds) using Wald formula
CIdel <- function(level, n, preval, dec.digits = 5){    
  prev <- preval
  prevalence <- as.numeric(prev)
  if(level == 0.95){
    z <- 1.96
  } else if(level == 0.99) {
    z <- 2.58
  } else if(level == 0.90) {
    z <- 1.645
  } else {
    print("error in the z, only 0.90, 0.95, 0.99 possible")
  }
  upper <- prevalence + z*sqrt((prevalence*(1-prevalence)/n))
  lower <- prevalence - z*sqrt((prevalence*(1-prevalence)/n))
  upper <- as.character(round(upper, digits = dec.digits))
  lower <- as.character(round(lower, digits = dec.digits))
  paste(lower, "-", upper)
}
#FUNCTION safe_min: to find min value checking NA
safe_min <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(min(x, na.rm = TRUE))
  }
}
#FUNCTION safe_max to find max value checking NA
safe_max <- function(x) {
  if (all(is.na(x))) {
    return(NA)
  } else {
    return(max(x, na.rm = TRUE))
  }
}

#1 DATA LOADING AND MERGING ----

#22 Years Old
HPV_cohortdata_22yo_20240131 <- readxl::read_excel("HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/HPV_cohort_data_Monica_9.1.2024_acctualized_20240816.xlsx", 
                                                   col_types = c("text", "text", "numeric", 
                                                                 "numeric", "date", "numeric", "date", 
                                                                 "numeric", "date", "text", "text",
                                                                 "date", "text", "text", "text", "text", 
                                                                 "text", "text", "text", "text", "text", 
                                                                 "text", "text", "text", "text", "text", 
                                                                 "text", "text", "text", "text", "text", 
                                                                 "text", "text", "text", "text", "text", 
                                                                 "text", "text", "text", "text", "text", 
                                                                 "date", "text", "text", "date", "numeric", 
                                                                 "numeric", "numeric", "text", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "date", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "numeric", "numeric", "numeric"), sheet = "22y")

cytology.22yo.update <- readxl::read_excel("HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/UpdateCytology_edited_To Monica 23.1.2024.xlsx", 
                                           sheet = "22yo", col_types = c("text", 
                                                                         "date", "text", "text", "text", "text", 
                                                                         "text", "text", "text", "date", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "text", "text", "text"))

#25 Years Old
HPV_cohortdata_25yo_20240123 <- readxl::read_excel("HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/HPV_cohort_data_Monica_9.1.2024_acctualized_20240816.xlsx", 
                                                   sheet = "25y", col_types = c("text", 
                                                                                "text", "numeric", "numeric", "numeric", 
                                                                                "date", "date", "text", "text", "text", 
                                                                                "text", "text", "text", "text", "text", 
                                                                                "text", "text", "text", "text", "text", 
                                                                                "text", "text", "text", "text", "text", 
                                                                                "text", "text", "text", "text", "text", 
                                                                                "text", "text", "text", "text", "text", 
                                                                                "text", "text", "date", "date", "numeric", 
                                                                                "text", "text", "text", "text", "numeric", 
                                                                                "numeric", "numeric", "text", "text", 
                                                                                "text", "numeric", "numeric", "numeric", 
                                                                                "text", "text", "text", "text", "text", 
                                                                                "numeric", "numeric", "numeric", 
                                                                                "text", "numeric", "text", "text", 
                                                                                "text", "text"))
cytology.25yo.update <- readxl::read_excel("HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/UpdateCytology_edited_To Monica 23.1.2024.xlsx", 
                                           sheet = "25yo", col_types = c("text", 
                                                                         "text", "date", "text", "text", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "text", "text", "text", "text", "text", 
                                                                         "text"))

HPVDNA.25yo.update <- readxl::read_excel("HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/HPV004_HPVDNA_25yo28yo_20240126 - kopia.xlsx", 
                                         sheet = "25v HPV-DNA", col_types = c("text", "text", "text", "numeric", "text", 
                                                                              "text", "text", "text", "numeric", 
                                                                              "numeric", "numeric", "text", "text", 
                                                                              "text", "numeric", "numeric", "numeric", 
                                                                              "text", "text", "text", "text", "text", 
                                                                              "numeric", "numeric", "numeric"))

#28 Years Old
HPV_cohortdata_28yo_20240123 <- readxl::read_excel("HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/HPV_cohort_data_Monica_9.1.2024_acctualized_20240816.xlsx", 
                                                   sheet = "28y2.0", col_types = c("text","text", "numeric", "numeric", "numeric", 
                                                                                   "date", "date", "date", "text", "text", 
                                                                                   "text", "text", "text", "text", "numeric", 
                                                                                   "numeric", "numeric", "numeric", 
                                                                                   "numeric", "numeric", "numeric", 
                                                                                   "text", "text", "text", "text", "text", "text", 
                                                                                   "text", "text", "text", "numeric", 
                                                                                   "numeric", "numeric", "numeric", 
                                                                                   "numeric", "numeric", "numeric", 
                                                                                   "text", "text", "text", "date", "text", 
                                                                                   "text", "text", "text", "text", "text", 
                                                                                   "numeric", "text", "text", "text", 
                                                                                   "text", "text", "text", "text", "text", 
                                                                                   "text", "text", "text", "text", "numeric", 
                                                                                   "numeric", "numeric", "text", "text", 
                                                                                   "text", "numeric", "numeric", "numeric", 
                                                                                   "numeric", "numeric", "text", "text", 
                                                                                   "text", "text", "text", "numeric"))

##     22 years old merging----

#Check for duplicates
length(unique(cytology.22yo.update$FID))
noccur <- data.frame(table(HPV_cohortdata_22yo_20240131$FID))
noccur[noccur$Freq > 1,] #OK
noccur <- data.frame(table(cytology.22yo.update$FID))
noccur[noccur$Freq > 1,] #OK

#CYTOLOGY Update
HPV_cohortdata_22yo_20240131.1 <- HPV_cohortdata_22yo_20240131[-c(12:41)]
HPV_cohortdata_22yo_newcyt <- merge(HPV_cohortdata_22yo_20240131.1, cytology.22yo.update, by = "FID", all.x = TRUE) 


##     25 years old merging----

#Check for duplicates
length(unique(HPV_cohortdata_25yo_20240123$FID))
length(unique(cytology.25yo.update$FID))
noccur <- data.frame(table(HPV_cohortdata_25yo_20240123$FID))
noccur[noccur$Freq > 1,] 
noccur <- data.frame(table(cytology.25yo.update$FID))
noccur[noccur$Freq > 1,] 

#CYTOLOGY Update (25)
HPV_cohortdata_25yo_20240123.1 <- HPV_cohortdata_25yo_20240123[-c(7:36)]
HPV_cohortdata_25yo_20240126_newcyt <- merge(HPV_cohortdata_25yo_20240123.1, cytology.25yo.update, by = "FID", all.x = TRUE) 

#HPV DNA Update (25)
HPV_cohortdata_25yo_20240126_newcyt.1 <- HPV_cohortdata_25yo_20240126_newcyt[-c(9:31)]
HPV_cohortdata_25yo_20240126_newcytnewDNA <- merge(HPV_cohortdata_25yo_20240126_newcyt.1, HPVDNA.25yo.update, by = "FID", all.x = TRUE) 

#add arm information from the 22 year olds data
arminfo <- HPV_cohortdata_22yo_newcyt %>% dplyr::select(FID, ARM_HPV004)
HPV_cohortdata_25yo_20240126_newcytnewDNA <- merge(HPV_cohortdata_25yo_20240126_newcytnewDNA, arminfo, by="FID", all.x = TRUE)

##     28 years old merging----

#Add arm information from 22 list
HPV_cohortdata_28yo_20240126_newcytnewDNA <- merge(HPV_cohortdata_28yo_20240123, arminfo, by="FID", all.x = TRUE)


# 2 CYTOLOGY AND HPV DNA FINDINGS PER VISIT ----

## Visit 22 ----

#Summarized number of individuals per arm and number of cytology and DNA results (22yo)
HPV_cohortdata_22yo_newcyt <- HPV_cohortdata_22yo_newcyt %>% dplyr::mutate(DNA_results = if_else(is.na(HPV_16)|HPV_16==9, "NO", "YES"))
HPV_cohortdata_22yo_newcyt %>% dplyr::filter(FID=="A50952")
N_arm <- ddply(HPV_cohortdata_22yo_newcyt, .(ARM_HPV004), summarize,
               n = length(na.omit(Visit)),
               nCIT = sum(Cytology1 != "NA", na.rm =TRUE),
               nDNA = sum(DNA_results =="YES", na.rm = TRUE))

#HPV DNA baseline findings (number and prevalence)
types.list <- c("16", "18", "31", "45", "51", "52", "P1|33|58","P2|56|59|66" , "P3|35|39|68")
typessumarylist <- list()
for (i in types.list){
  HPVtype.22 <- HPV_cohortdata_22yo_newcyt %>% filter(grepl(i, `COR Positive Genotypes`)==TRUE)
  typessumarylist[[i]] <- ddply(HPVtype.22, .(ARM_HPV004), summarize,
                                nHPVtype = length(FID))
}
typessumarylist_prct <- list()
for (i in names(typessumarylist)){
  typessumarylist_prct[[i]] <- cbind(typessumarylist[[i]],N_arm["nDNA"]) %>%  cbind(N_arm["nCIT"])
  typessumarylist_prct[[i]] <- typessumarylist_prct[[i]] %>% dplyr::mutate(percentDNA = nHPVtype/nDNA*100)
}
#total hrHPV
pos22 <- HPV_cohortdata_22yo_newcyt %>% dplyr::filter(!is.na(`COR Positive Genotypes`)) 
totalHPV_22yo <- as.data.frame(table(pos22$ARM_HPV004)) %>% cbind(N_arm["nDNA"]) %>% dplyr::mutate(precent_totalHPV = Freq/nDNA*100)
#missing from HPV analysis
N_arm %>% dplyr::mutate(missingDNAan = n-nDNA) %>% dplyr::mutate(missingCYTan = n-nCIT)

#Bar plot with table of prevalence of baseline findings 22yo
HPVDNAprev_forbarplot <- readxl::read_excel("HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/HPVprevalencefindings_forbarplot_20250218.xlsx", 
                                            col_types = c("text", "text", "numeric", 
                                                          "numeric", "numeric", "numeric", 
                                                          "numeric", "numeric"))

labelstypesX <- c("HPV16","HPV18", "HPV31","HPV45", "HPV51", "HPV52" , "HPV33/58", "HPV35/39/68", "HPV56/59/66", "Total hrHPV")
HPVDNAprev_forbarplot %>% ggplot(aes(x = HPVtype, y = prevalence22, fill = Arm)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  labs(title = "hrHPV prevalence at 1st screening visit", y = "Prevalence (%)", x = "HPV Type") +
  scale_fill_manual(name = NULL , labels = c("Frequent Screening arm", "Infrequent screening arm", "Safety arm"),
                    #values = c( '#0072b2', '#e69f00', '#009e73')) +
                    values = c( '#800080', '#a0e8eb', '#f0b7ef')) +
  scale_x_discrete(limits = labelstypesX, labels = labelstypesX) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45), text = element_text(family = "serif", face=NULL, size=15))

#Cytology baseline findings
HPV_cohortdata_22yo_newcyt2 <- HPV_cohortdata_22yo_newcyt %>% dplyr::mutate(xCytology1 = case_when(grepl("Normal|normal", Cytology1)==TRUE ~ "Normal",
                                                                                                   Cytology1=="HSIL" ~ "HSIL",
                                                                                                   grepl("AGCNOS", Cytology1)==TRUE ~ "AGC-NOS",
                                                                                                   Cytology1=="ASCH" ~ "ASC-H",
                                                                                                   Cytology1=="ASCUS" ~ "ASCUS",
                                                                                                   Cytology1=="LSIL" ~ "LSIL",
                                                                                                   Cytology1 =="ASCUS (or LSIL)" ~ "ASCUS/LSIL",
                                                                                                   Cytology1 == "ASCUS_LSIL" ~ "ASCUS/LSIL",
                                                                                                   Cytology1=="Inadequate" ~ "Inadequate"))
as.data.frame(table(HPV_cohortdata_22yo_newcyt$Cytology1))
as.data.frame(table(HPV_cohortdata_22yo_newcyt2$xCytology1))

CITtypessummarylist <- list()
for (i in names(table(HPV_cohortdata_22yo_newcyt2$xCytology1))){
  CITtype.22 <- HPV_cohortdata_22yo_newcyt2 %>% dplyr::filter(xCytology1 == i)
  CITtypessummarylist[[i]] <- ddply(CITtype.22, .(ARM_HPV004), summarize,
                                    nCITtype = length(FID))
  CITtypessummarylist[[i]] <- merge(CITtypessummarylist[[i]],N_arm, by="ARM_HPV004" )
  CITtypessummarylist[[i]] <- CITtypessummarylist[[i]] %>% dplyr::select("ARM_HPV004",  "nCITtype", "nCIT")
}
CITtypessummarylist.pct <- list()
for (i in names(table(HPV_cohortdata_22yo_newcyt2$xCytology1))){
  CITtypessummarylist.pct[[i]] <- CITtypessummarylist[[i]] %>% dplyr::mutate(percentCIT = nCITtype/nCIT*100)
}

#put ASCUS/LSIL into ASCUS and LSIL separately (sum cases)
ASCUS <- cbind(CITtypessummarylist[["ASCUS"]], CITtypessummarylist[["ASCUS/LSIL"]]["nCITtype"]) 
colnames(ASCUS) <- c("ARM_HPV004", "nCITtype1","nCIT", "nCITtype2")
ASCUS <- ASCUS %>% dplyr::mutate(nCITtype = nCITtype1+nCITtype2) %>% dplyr::mutate(percentCIT = nCITtype/nCIT*100)
LSIL <- cbind(CITtypessummarylist[["LSIL"]], CITtypessummarylist[["ASCUS/LSIL"]]["nCITtype"])
colnames(LSIL) <- c("ARM_HPV004", "nCITtype1","nCIT", "nCITtype2")
LSIL <- LSIL %>% dplyr::mutate(nCITtype = nCITtype1+nCITtype2) %>% dplyr::mutate(percentCIT = nCITtype/nCIT*100)

#Total abnormal by arm
totalabnormal <- as.data.frame(table(HPV_cohortdata_22yo_newcyt2$xCytology1))[-8,] #Take out the "Normal" Row
sum(totalabnormal$Freq)

totalabnormalbyarm <- list()
for (i in c("A1", "A2", "A3")){
  dataarm <- HPV_cohortdata_22yo_newcyt2 %>% dplyr::filter(ARM_HPV004 == i)
  dataarmfreq <- as.data.frame(table(dataarm$xCytology1))
  dataarmfreq <- dataarmfreq %>% dplyr::filter(!Var1 =="Normal")
  totalabnormalbyarm[[i]] <- sum(dataarmfreq$Freq)
}
totalabnormalbyarm.df <- t(as.data.frame(totalabnormalbyarm))
cbind(N_arm, totalabnormalbyarm.df) %>% dplyr::mutate(percentTAB = totalabnormalbyarm.df/nCIT*100 )

## Visit 25 ----

#DNA results "YES" "NO"
sum(is.na(HPV_cohortdata_25yo_20240126_newcytnewDNA$`HPV HR Result`))
HPV_cohortdata_25yo_20240126_newcytnewDNA.1 <- HPV_cohortdata_25yo_20240126_newcytnewDNA %>% dplyr::mutate(DNA_results = if_else(is.na(`HPV HR Result`)== FALSE, "YES", "NO"))
table(HPV_cohortdata_25yo_20240126_newcytnewDNA.1$DNA_results)

#summarized general numbers of cytology and HPV DNA results
N_arm25 <- ddply(HPV_cohortdata_25yo_20240126_newcytnewDNA.1, .(ARM_HPV004), summarize,
                 n = length(Visit),
                 #nCIT = sum(Cytology1 != "NA", na.rm = TRUE),  ***
                 nDNA = sum(DNA_results =="YES", na.rm = TRUE))
N_arm25 <- N_arm25[-4,]
N_arm25%>%dplyr::mutate(misDNA=n-nDNA)

#***this doesn't account for the stated as "missing" or "inadequate". Real numbers later (see #478) 
as.data.frame(table(HPV_cohortdata_25yo_20240126_newcytnewDNA.1$Cytology1))

#HPV DNA findings (number and prevalence
types.list25 <- c("16", "18", "31", "45", "51", "52", "P1|33|58", "P2|56|59|66", "P3|35|39|68")
typessumarylist25 <- list()
for (i in types.list25){
  HPVtype.25 <- HPV_cohortdata_25yo_20240126_newcytnewDNA.1 %>% filter(grepl(i, `COR Positive Genotypes`)==TRUE)
  typessumarylist25[[i]] <- ddply(HPVtype.25, .(ARM_HPV004), summarize,
                                  nHPVtype = length(FID))
}
typessumarylist_prct25 <- list()
for (i in names(typessumarylist25)){
  typessumarylist_prct25[[i]] <- cbind(typessumarylist25[[i]],N_arm25["nDNA"]) 
  typessumarylist_prct25[[i]] <- typessumarylist_prct25[[i]] %>% dplyr::mutate(percentDNA = nHPVtype/nDNA*100)
}
#total hrHPV
pos25 <- HPV_cohortdata_25yo_20240126_newcytnewDNA.1 %>% dplyr::filter(!is.na(`COR Positive Genotypes`)) 
totalHPV_25yo <- as.data.frame(table(pos25$ARM_HPV004)) %>% cbind(N_arm25["nDNA"]) %>% dplyr::mutate(precent_totalHPV = Freq/nDNA*100)

#Bar plot with table of prevalence of HPV findings 25yo
HPVDNAprev_forbarplot
labelstypesX <- c("HPV16","HPV18", "HPV31","HPV45", "HPV51", "HPV52" , "HPV33/58", "HPV35/39/68", "HPV56/59/66", "Total hrHPV")
HPVDNAprev_forbarplot %>% ggplot(aes(x = HPVtype, y = prevalence25, fill = Arm)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  labs(title = "hrHPV prevalence at 2nd screening visit", y = "Prevalence (%)", x = "HPV Type") +
  scale_fill_manual(name = NULL , labels = c("Frequent Screening arm", "Infrequent screening arm", "Safety arm"),
                    #values = c( '#0072b2', '#e69f00', '#009e73')) +
                    values = c( '#800080', '#a0e8eb', '#f0b7ef')) +
  scale_x_discrete(limits = labelstypesX, labels = labelstypesX) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45), text = element_text(family = "serif", face=NULL, size=15))

## Visit 28 ----

#DNA results "YES" "NO"
HPV_cohortdata_28yo_20240126_newcytnewDNA.1 <- HPV_cohortdata_28yo_20240126_newcytnewDNA %>% 
  dplyr::mutate(DNA_results = case_when(is.na(`HPV HR Result`)== TRUE ~ "NO",`HPV HR Result` == "POS"~ "YES",
                                        `HPV HR Result` == "NEG"~ "YES", `HPV HR Result` == "NA"~ "NO"))
table(HPV_cohortdata_28yo_20240126_newcytnewDNA.1$DNA_results)
sum(is.na(HPV_cohortdata_28yo_20240126_newcytnewDNA.1$DNA_results))

table(HPV_cohortdata_28yo_20240126_newcytnewDNA.1$Cytology1)
sum(is.na(HPV_cohortdata_28yo_20240126_newcytnewDNA.1$Cytology1))

#general information 
N_arm28 <- ddply(HPV_cohortdata_28yo_20240126_newcytnewDNA.1, .((ARM_HPV004)), summarize,
                 n = length(FID),
                nCIT = sum(Cytology1 != "NA", na.rm = TRUE),##CYT is not exactly correct. It doesn't count possible "missing" and "inadequate" (see #478) 
                 nDNA = sum(DNA_results =="YES", na.rm = TRUE))
N_arm28 <-N_arm28[-4,]
N_arm28A1A2 <- N_arm28[-3,] 
N_arm28%>%dplyr::mutate(misDNA=n-nDNA)

#HPV DNA findings 28 (number and prevalence) 
types.list28 <- c("16", "18", "31", "45", "51", "52", "P1|33|58", "P2|56|59|66","P3|35|39|68")
typessumarylist28 <- list()
for (i in types.list28){
  HPVtype.28 <- HPV_cohortdata_28yo_20240126_newcytnewDNA.1 %>% filter(grepl(i, `COR Positive Genotypes`)==TRUE)
  typessumarylist28[[i]] <- ddply(HPVtype.28, .(ARM_HPV004), summarize,
                                  nHPVtype = length(FID))
}
typessumarylist_prct28 <- list()
for (i in names(typessumarylist28)){
  typessumarylist_prct28[[i]] <- cbind(typessumarylist28[[i]],N_arm28A1A2["nDNA"]) %>%  cbind(N_arm28A1A2["nCIT"])
  typessumarylist_prct28[[i]] <- typessumarylist_prct28[[i]] %>% dplyr::mutate(percentDNA = nHPVtype/nDNA*100)
}
#total hrHPV
pos28 <- HPV_cohortdata_28yo_20240126_newcytnewDNA.1 %>% dplyr::filter(!is.na(`COR Positive Genotypes`)) %>%dplyr::filter(!`COR Positive Genotypes` =="NA") 
totalHPV_28yo <- as.data.frame(table(pos28$ARM_HPV004)) %>% cbind(N_arm28A1A2["nDNA"]) %>% dplyr::mutate(precent_totalHPV = Freq/nDNA*100)

#Bar plot with table of prevalence of HPV findings 28yo
HPVDNAprev_forbarplot
labelstypesX <- c("HPV16","HPV18", "HPV31","HPV45", "HPV51", "HPV52" , "HPV33/58", "HPV35/39/68", "HPV56/59/66", "Total hrHPV")
HPVDNAprev_forbarplot %>% ggplot(aes(x = HPVtype, y = prevalence28, fill = Arm)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  labs(title = "hrHPV prevalence at 3rd screening visit", y = "Prevalence (%)", x = "HPV Type") +
  scale_fill_manual(name = NULL , labels = c("Frequent Screening arm", "Infrequent screening arm", "Safety arm"),
                    values = c( '#800080', '#a0e8eb', '#f0b7ef')) +
  scale_x_discrete(limits = labelstypesX, labels = labelstypesX) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45), text = element_text(family = "serif", face=NULL, size=15))


# 3 MERGING ALL 3 DATASETS 22-25-28 ----

HPV_cohortdata_allages_2024.1 <- merge(HPV_cohortdata_22yo_newcyt,HPV_cohortdata_25yo_20240126_newcytnewDNA, by="FID", all.x = TRUE )
HPV_cohortdata_allages_20240201 <- merge(HPV_cohortdata_allages_2024.1,HPV_cohortdata_28yo_20240126_newcytnewDNA, by="FID", all.x = TRUE )
dim(HPV_cohortdata_allages_20240201)

#Histology: diagnosis variables, CIN1 and CIN2.3+
HPV_cohortdata_allages_20240201.x <- HPV_cohortdata_allages_20240201 %>% 
  mutate(Diagnosis = case_when(!is.na(Histology1)==TRUE ~ Histology1,
                               !is.na(Histology1.y)==TRUE ~ Histology1.y,
                               !is.na(Histology1.x) == TRUE ~ Histology1.x,
                               TRUE ~ "NA"))  %>% 
 dplyr::mutate(CIN2.3 = case_when(grepl("CIN2", Diagnosis)==TRUE ~"YES",
                                  grepl("CIN3", Diagnosis)==TRUE ~"YES",
                                  grepl("CIN 2", Diagnosis)==TRUE ~"YES",
                                  grepl("HSIL", Diagnosis)==TRUE ~"YES",
                                  grepl("CIN 1-2", Diagnosis)==TRUE ~"YES",
                                  grepl("CIN II", Diagnosis)==TRUE ~"YES",
                                  grepl("CIN 3", Diagnosis)==TRUE ~"YES",
                                  TRUE ~ "NO")) %>%
  dplyr::mutate(CIN1 = case_when(grepl("CIN1", Diagnosis)==TRUE ~"YES",
                                 grepl("CIN 1", Diagnosis)==TRUE ~"YES",
                                 grepl("LSIL", Diagnosis)==TRUE ~"YES",
                                 TRUE ~ "NO")) %>%
  dplyr::select(!ARM_HPV004, ARM_HPV004.y) #OBS!: ARM_HPV004.x is the real arm information!

#Accounting for the CYT that are missing or inadequate or fake NAs (var: CYTvisit)
HPV_cohortdata_allages_20240201.yA1A2A3 <- HPV_cohortdata_allages_20240201.x %>% 
  dplyr::mutate(CYTvisit22 = case_when(!is.na(Cytology1.x)==TRUE ~ "YES",
                                       TRUE ~ "NO")) %>%
  dplyr::mutate(CYTvisit25 = case_when(Cytology1.y =="Inadequate" ~ "NO",
                                       Cytology1.y !="NA" ~ "YES",
                                       TRUE ~ "NO")) %>%
  dplyr::mutate(CYTvisit28 = case_when(Cytology1 =="missing" ~ "NO",
                                       Cytology1 =="Inadequate" ~ "NO",
                                       Cytology1 !="NA" ~ "YES",
                                       TRUE ~ "NO")) %>%
  dplyr::mutate(Histologyvisit = case_when(!is.na(`Date _of_histology.x`) ~ "Hist_v22" ,   
                                           !is.na(`Date _of_histology.y`) ~ "Hist_v25",  
                                           !is.na(`Date _of_histology`) ~ "Hist_v28",
                                           TRUE ~ NA)) 

ddply(HPV_cohortdata_allages_20240201.yA1A2A3, .(ARM_HPV004.x),summarise,
      CYT22 = sum(CYTvisit22 == "YES"),
      CYT25 = sum(CYTvisit25 == "YES"),
      CYT28 = sum(CYTvisit28 == "YES"))


## Combined age: prevalence and prevalence ratio CIN2+ (binomial regression model) ----
HPV_cohortdata_allages_20240201.xA1A2 <- HPV_cohortdata_allages_20240201.yA1A2A3 %>% dplyr::filter(!ARM_HPV004.x == "A3") %>% 
  dplyr::mutate(binomial_ARM = if_else(ARM_HPV004.x=="A1", 0, 1)) %>% #A1 = 0 / A2 = 1
  dplyr::mutate(binomial_CIN2.3 = if_else(CIN2.3 == "YES", 1, 0))      #CIN2+ = 1 / NO = 0

table(HPV_cohortdata_allages_20240201.yA1A2A3$Cytology1)

ddply(HPV_cohortdata_allages_20240201.xA1A2, .(binomial_ARM, Histologyvisit), summarise,
      n.cases= sum(binomial_CIN2.3 ==1),
      n = length(FID))

model_combinedages <- glm(HPV_cohortdata_allages_20240201.xA1A2$binomial_CIN2.3 ~ HPV_cohortdata_allages_20240201.xA1A2$binomial_ARM, family=binomial(link="log")) #log gives prevalence ratio if i used logit it would give me odds ratio
summary(model_combinedages)
ci.exp(model_combinedages, alpha = 0.05)

## Per visit: Prevalence and prevalence ratio CIN2+ (binomial regression model) -----
HPV_cohortdata_allages_20240201.yA1A2.2 <- HPV_cohortdata_allages_20240201.xA1A2
HPV_cohortdata_allages_20240201.yA1A2.2[which(HPV_cohortdata_allages_20240201.yA1A2.2$FID == "A21159"),"Histologyvisit"] <- "Hist_v22" #missing date
HPV_cohortdata_allages_20240201.yA1A2.2[which(HPV_cohortdata_allages_20240201.yA1A2.2$FID == "A41531"),"Histologyvisit"] <- "Hist_v28" #missing date
HPV_cohortdata_allages_20240201.yA1A2.2[which(HPV_cohortdata_allages_20240201.yA1A2.2$FID == "B20837"),"Histologyvisit"] <- "Hist_v25" #missing date
HPV_cohortdata_allages_20240201.yA1A2.2[which(HPV_cohortdata_allages_20240201.yA1A2.2$FID == "B31244"),"Histologyvisit"] <- "Hist_v28" #missing date
HPV_cohortdata_allages_20240201.yA1A2x <- HPV_cohortdata_allages_20240201.yA1A2.2 %>% dplyr::mutate(CIN2.3_22 = if_else(CIN2.3=="YES"& Histologyvisit=="Hist_v22", 1, 0)) %>%
  dplyr::mutate(CIN2.3_25 = if_else(CIN2.3=="YES" & Histologyvisit=="Hist_v25", 1, 0)) %>%
  dplyr::mutate(CIN2.3_28 = if_else(CIN2.3=="YES" & Histologyvisit=="Hist_v28", 1, 0))
HPV_cohortdata_allages_20240201.yA1A2x %>% dplyr::filter(CIN2.3=="YES" & is.na(Histologyvisit)) %>% select(FID, `Date _of_histology.x`, `Date _of_histology.y`, `Date _of_histology`, CYTvisit22, CYTvisit25, CYTvisit28, Histologyvisit, CIN2.3)

HPV_cohortdata_allages_20240201.yA1A2_22yo <- HPV_cohortdata_allages_20240201.yA1A2x %>% dplyr::filter(CYTvisit22 == "YES")
HPV_cohortdata_allages_20240201.yA1A2_25yo <- HPV_cohortdata_allages_20240201.yA1A2x %>% dplyr::filter(CYTvisit25 == "YES")
HPV_cohortdata_allages_20240201.yA1A2_28yo <- HPV_cohortdata_allages_20240201.yA1A2x %>% dplyr::filter(CYTvisit28 == "YES") #Selecting by these this will not count the one CIN2+s outside of the study 


model_22yo <- glm(HPV_cohortdata_allages_20240201.yA1A2_22yo$CIN2.3_22 ~ HPV_cohortdata_allages_20240201.yA1A2_22yo$binomial_ARM, family=binomial(link="log")) #log gives prevalence ratio if i used logit it would give me odds ratio
summary(model_22yo)
ci.exp(model_22yo, alpha = 0.05)

model_25yo <- glm(HPV_cohortdata_allages_20240201.yA1A2_25yo$CIN2.3_25 ~ HPV_cohortdata_allages_20240201.yA1A2_25yo$binomial_ARM, family=binomial(link="log")) #log gives prevalence ratio if i used logit it would give me odds ratio
summary(model_25yo)
ci.exp(model_25yo, alpha = 0.05)

model_28yo <- glm(HPV_cohortdata_allages_20240201.yA1A2_28yo$CIN2.3_28 ~ HPV_cohortdata_allages_20240201.yA1A2_28yo$binomial_ARM, family=binomial(link="log")) #log gives prevalence ratio if i used logit it would give me odds ratio
summary(model_28yo)
ci.exp(model_28yo, alpha = 0.05)

#Calculate the prevalence ratio and CI for CIN2+ by arm and by age of diagnosis
sum22 <- ddply(HPV_cohortdata_allages_20240201.yA1A2_22yo, .(ARM_HPV004.x), summarise,
               CIN2= sum(na.omit(CIN2.3_22=="1")),
               n=length(FID)) %>% dplyr::mutate(visit="22yo")
sum25 <- ddply(HPV_cohortdata_allages_20240201.yA1A2_25yo, .(ARM_HPV004.x), summarise,
               CIN2= sum(na.omit(CIN2.3_25=="1")),
               n=length(FID)) %>% dplyr::mutate(visit="25yo")
sum28 <- ddply(HPV_cohortdata_allages_20240201.yA1A2_28yo, .(ARM_HPV004.x), summarise,
               CIN2= sum(na.omit(CIN2.3_28=="1")),
               n=length(FID)) %>% dplyr::mutate(visit="28yo")
sumary22.25.28 <- rbind(sum22,sum25, sum28)
sumary22.25.28 <- sumary22.25.28%>% dplyr::mutate(prev.CIN2.3 = CIN2/n) %>% mutate(CI = CIdel(0.95, n, prev.CIN2.3))
cat("PR 22yo:", sumary22.25.28$prev.CIN2.3[2]/sumary22.25.28$prev.CIN2.3[1],"/",
    "PR 25yo:",sumary22.25.28$prev.CIN2.3[4]/sumary22.25.28$prev.CIN2.3[3], "/",
    "PR 28yo:", sumary22.25.28$prev.CIN2.3[6]/sumary22.25.28$prev.CIN2.3[5])

#Calculate the prevalence ratio and CI by arm (combined age of diagnosis)
HPV_cohortdata_allages_20240201.yA1A2x %>% dplyr::filter(CIN2.3=="YES"&Histologyvisit=="Hist_v28") %>% 
  select(FID, `Date _of_histology.x`, `Date _of_histology.y`, `Date _of_histology`, CYTvisit22, CYTvisit25, CYTvisit28, 
         Histologyvisit, CIN2.3) ##One individual that did not attend visit 28 but outside of the study (7). Thats why donesnt mach with all combined arms. 

HPV_forPRyA1A2x <- HPV_cohortdata_allages_20240201.yA1A2x %>% dplyr::mutate(CIN2allvisits = CIN2.3_22+CIN2.3_25+CIN2.3_28) 
CIN_perARMallages <- ddply(HPV_forPRyA1A2x, .(ARM_HPV004.x), summarise,
                           CIN2= sum(na.omit(CIN2allvisits==1)),
                           n=length(FID)) 
CIN_perARMallages <- CIN_perARMallages %>% mutate(prev.CIN2.3 = CIN2/n) %>% mutate(CI = CIdel(0.95, n, prev.CIN2.3))
cat("PR combined:", CIN_perARMallages$prev.CIN2.3[2]/CIN_perARMallages$prev.CIN2.3[1])


## Table of all histology cases (with HPV DNA and Cyt findings) ----
Allhistologtcases_A1A2A3 <- HPV_cohortdata_allages_20240201.yA1A2A3 %>% mutate(All.CYT = case_when(!is.na(Cytology1.x) ~ Cytology1.x ,   
                                                                                                           !is.na(Cytology1.y) ~ Cytology1.y,  
                                                                                                           !is.na(Cytology1) ~ Cytology1,
                                                                                                           TRUE ~ "Not Available")) %>% dplyr::mutate(Diagnosis.date = case_when(!is.na(`Date _of_histology.x`) ~ `Date _of_histology.x`,   
                                                                                                                                                                                 !is.na(`Date _of_histology.y`) ~ `Date _of_histology.y`,  
                                                                                                                                                                                 !is.na(`Date _of_histology`) ~ `Date _of_histology`)) %>%
  dplyr::select(FID, ARM_HPV004.x, Birth_year.x, dose_1, `COR Positive Genotypes.x`, `COR Positive Genotypes.y`, `COR Positive Genotypes`,Cytology1.x,Cytology1.y ,Cytology1, All.CYT, Diagnosis, Diagnosis.date, CIN2.3, CIN1) %>% 
  dplyr::filter(!Diagnosis == "NA")

colnames(Allhistologtcases_A1A2A3) <- c("FID", "Arm", "Birth_year", "Date.of.VAccination", "Genotyping22", "Genotyping25", "Genotyping28",
                                                "Cytology22", "Cytology25", "Cytology28", "AllCytology", "Diagnosis", "Diagnosis.date", "CIN2+", "CIN1")

onlycin2_HISA3 <- Allhistologtcases_A1A2A3 %>% dplyr::filter(`CIN2+`=="YES")
write.csv(onlycin2_HISA3, file = "HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/CIN2plusHSIL_findings_A1A2A3_20250505.csv")


## Combined age: Prevalence and prevalence ratios of Cytology findings (binomial regression model) -----

VisitN <- ddply(HPV_cohortdata_allages_20240201.xA1A2, .(ARM_HPV004.x), summarize,
                Visit22 = sum(CYTvisit22 =="YES"),
                Visit25 = sum(CYTvisit25 =="YES"),
                Visit28 = sum(CYTvisit28 == "YES")) 
CYTperarm_prev22 <- ddply(HPV_cohortdata_allages_20240201.xA1A2, .(ARM_HPV004.x), summarize,
                          ASCUS_22 = sum(grepl("ASCUS", Cytology1.x)),
                          LSIL_22 = sum(grepl("LSIL", Cytology1.x)),
                          HSIL_22 = sum(grepl("HSIL", Cytology1.x))) %>% dplyr::mutate(N=c(2713,2865))%>%
  dplyr::mutate(prevASCUS=ASCUS_22/N)%>%dplyr::mutate(prevLSIL=LSIL_22/N)%>%dplyr::mutate(prevHSIL=HSIL_22/N) %>%
  dplyr::mutate(CI_ASCUS=(CIdel(0.95, N, prevASCUS)))%>%dplyr::mutate(CI_LSIL=CIdel(0.95, N, prevLSIL))%>%dplyr::mutate(CI_HSIL=CIdel(0.95, N, prevHSIL))

CYTperarm_prev25 <- ddply(HPV_cohortdata_allages_20240201.xA1A2, .(ARM_HPV004.x), summarize,
                          ASCUS_25 = sum(grepl("ASCUS", Cytology1.y)),
                          LSIL_25 = sum(grepl("LSIL", Cytology1.y)),
                          HSIL_25 = sum(grepl("HSIL", Cytology1.y))) %>% dplyr::mutate(N=c(2423,2601))%>%
  dplyr::mutate(prevASCUS=ASCUS_25/N)%>%dplyr::mutate(prevLSIL=LSIL_25/N)%>%dplyr::mutate(prevHSIL=HSIL_25/N) %>%
  dplyr::mutate(CI_ASCUS=(CIdel(0.95, N, prevASCUS)))%>%dplyr::mutate(CI_LSIL=CIdel(0.95, N, prevLSIL))%>%dplyr::mutate(CI_HSIL=CIdel(0.95, N, prevHSIL))

CYTperarm_prev28 <- ddply(HPV_cohortdata_allages_20240201.xA1A2, .(ARM_HPV004.x), summarize,
                          ASCUS_28 = sum(grepl("ASCUS", Cytology1)),
                          LSIL_28 = sum(grepl("LSIL", Cytology1)),
                          HSIL_28 = sum(grepl("HSIL", Cytology1))) %>% dplyr::mutate(N=c(2245,2428))%>%
  dplyr::mutate(prevASCUS=ASCUS_28/N)%>%dplyr::mutate(prevLSIL=LSIL_28/N)%>%dplyr::mutate(prevHSIL=HSIL_28/N) %>%
  dplyr::mutate(CI_ASCUS=(CIdel(0.95, N, prevASCUS)))%>%dplyr::mutate(CI_LSIL=CIdel(0.95, N, prevLSIL))%>%dplyr::mutate(CI_HSIL=CIdel(0.95, N, prevHSIL))

#Combined ages: prevalence 
HPVcohort_allcyt.forcombinedprev <- HPV_cohortdata_allages_20240201.xA1A2 %>% dplyr::mutate(All.CYT.mesh = paste(Cytology1.x,Cytology1.y, Cytology1))
CYTperarm_allages <- ddply(HPVcohort_allcyt.forcombinedprev, .(ARM_HPV004.x), summarise,
                           ASCUS = sum(grepl("ASCUS", All.CYT.mesh)),
                           LSIL = sum(grepl("LSIL", All.CYT.mesh)),
                           HSIL = sum(grepl("HSIL", All.CYT.mesh))) %>% dplyr::mutate(N=c(2713,2865))%>%
  dplyr::mutate(prevASCUS=ASCUS/N)%>%dplyr::mutate(prevLSIL=LSIL/N)%>%dplyr::mutate(prevHSIL=HSIL/N) %>%
  dplyr::mutate(CI_ASCUS=(CIdel(0.95, N, prevASCUS)))%>%dplyr::mutate(CI_LSIL=CIdel(0.95, N, prevLSIL))%>%dplyr::mutate(CI_HSIL=CIdel(0.95, N, prevHSIL))

#Combined ages prevalence ratio
HPVcohort_allcyt.forcombinedprev <- HPVcohort_allcyt.forcombinedprev %>% dplyr::mutate(ASCUSyesno = case_when(grepl("ASCUS", All.CYT.mesh)==TRUE~"YES",
                                                                                                              TRUE~"NO")) %>%
  dplyr::mutate(LSILyesno = case_when(grepl("LSIL", All.CYT.mesh)==TRUE ~"YES",
                                      TRUE~"NO")) %>%
  dplyr::mutate(HSILyesno = case_when(grepl("HSIL", All.CYT.mesh)==TRUE~"YES",
                                      TRUE~"NO"))
HPVcohort_allcyt.forcombinedprev.bin <- HPVcohort_allcyt.forcombinedprev %>% dplyr::mutate(bin_ASCUS = if_else(ASCUSyesno =="YES", 1, 0)) %>%    #ASCUS LSIL and HSIL yes/no to 1/0 (binomial for model)
  dplyr::mutate(bin_HSIL = if_else(HSILyesno =="YES", 1, 0)) %>%
  dplyr::mutate(bin_LSIL = if_else(LSILyesno =="YES", 1, 0))
ASCUSmodel <- glm(HPVcohort_allcyt.forcombinedprev.bin$bin_ASCUS ~ HPVcohort_allcyt.forcombinedprev.bin$binomial_ARM, family=binomial(link="log")) #log gives prevalence ratio if i used logit it would give me odds ratio
summary(ASCUSmodel)
ci.exp(ASCUSmodel, alpha = 0.05)
LSILmodel <- glm(HPVcohort_allcyt.forcombinedprev.bin$bin_LSIL ~ HPVcohort_allcyt.forcombinedprev.bin$binomial_ARM, family=binomial(link="log")) #log gives prevalence ratio if i used logit it would give me odds ratio
summary(LSILmodel)
ci.exp(LSILmodel, alpha = 0.05)
HSILmodel <- glm(HPVcohort_allcyt.forcombinedprev.bin$bin_HSIL ~ HPVcohort_allcyt.forcombinedprev.bin$binomial_ARM, family=binomial(link="log")) #log gives prevalence ratio if i used logit it would give me odds ratio
summary(HSILmodel)
ci.exp(HSILmodel, alpha = 0.05)

## Per visit: Prevalence and prevalence ratios of Cytology findings (binomial regression model) -----

HPV_cohortdata_cytfindings_bin22 <- HPV_cohortdata_allages_20240201.xA1A2 %>% dplyr::filter(CYTvisit22 == "YES")
HPV_cohortdata_cytfindings_bin25 <- HPV_cohortdata_allages_20240201.xA1A2 %>% dplyr::filter(CYTvisit25 == "YES")
HPV_cohortdata_cytfindings_bin28 <- HPV_cohortdata_allages_20240201.xA1A2 %>% dplyr::filter(CYTvisit28 == "YES")
ddply(HPV_cohortdata_cytfindings_bin28, .(binomial_ARM), summarise, #Check
      finding = sum(grepl("LSIL", Cytology1)==TRUE),
      n=length(FID))

#22 year old
HPV_cohortdata_cytfindings_bin22 <- HPV_cohortdata_cytfindings_bin22 %>% dplyr::mutate(ASCUS22 = if_else(grepl("ASCUS", Cytology1.x)==TRUE, 1, 0)) %>%
  dplyr::mutate(LSIL22 = if_else(grepl("LSIL", Cytology1.x)==TRUE, 1, 0)) %>%
  dplyr::mutate(HSIL22 = if_else(grepl("HSIL", Cytology1.x)==TRUE, 1, 0))
#ASCUS 22
ASCUSmodel_22 <- glm(HPV_cohortdata_cytfindings_bin22$ASCUS22 ~ HPV_cohortdata_cytfindings_bin22$binomial_ARM, family=binomial(link="log")) 
summary(ASCUSmodel_22)
ci.exp(ASCUSmodel_22, alpha = 0.05)
#LSIL 22
LSILmodel_22 <- glm(HPV_cohortdata_cytfindings_bin22$LSIL22 ~ HPV_cohortdata_cytfindings_bin22$binomial_ARM, family=binomial(link="log")) 
summary(LSILmodel_22)
ci.exp(LSILmodel_22, alpha = 0.05)
#HSIL 22
HSILmodel_22 <- glm(HPV_cohortdata_cytfindings_bin22$HSIL22 ~ HPV_cohortdata_cytfindings_bin22$binomial_ARM, family=binomial(link="log")) 
summary(HSILmodel_22)
ci.exp(HSILmodel_22, alpha = 0.05)

#25 year old
HPV_cohortdata_cytfindings_bin25 <- HPV_cohortdata_cytfindings_bin25 %>% dplyr::mutate(ASCUS25 = if_else(grepl("ASCUS", Cytology1.y)==TRUE, 1, 0)) %>%
  dplyr::mutate(LSIL25 = if_else(grepl("LSIL", Cytology1.y)==TRUE, 1, 0)) %>%
  dplyr::mutate(HSIL25 = if_else(grepl("HSIL", Cytology1.y)==TRUE, 1, 0))
#ASCUS 25
ASCUSmodel_25 <- glm(HPV_cohortdata_cytfindings_bin25$ASCUS25 ~ HPV_cohortdata_cytfindings_bin25$binomial_ARM, family=binomial(link="log")) 
summary(ASCUSmodel_25)
ci.exp(ASCUSmodel_25, alpha = 0.05)
#LSIL 25
LSILmodel_25 <- glm(HPV_cohortdata_cytfindings_bin25$LSIL25 ~ HPV_cohortdata_cytfindings_bin25$binomial_ARM, family=binomial(link="log")) 
summary(LSILmodel_25)
ci.exp(LSILmodel_25, alpha = 0.05)
#HSIL 25
HSILmodel_25 <- glm(HPV_cohortdata_cytfindings_bin25$HSIL25 ~ HPV_cohortdata_cytfindings_bin25$binomial_ARM, family=binomial(link="log")) 
summary(HSILmodel_25)
ci.exp(HSILmodel_25, alpha = 0.05)

#28 year old
HPV_cohortdata_cytfindings_bin28 <- HPV_cohortdata_cytfindings_bin28 %>% dplyr::mutate(ASCUS28 = if_else(grepl("ASCUS", Cytology1)==TRUE, 1, 0)) %>%
  dplyr::mutate(LSIL28 = if_else(grepl("LSIL", Cytology1)==TRUE, 1, 0)) %>%
  dplyr::mutate(HSIL28 = if_else(grepl("HSIL", Cytology1)==TRUE, 1, 0))
#ASCUS 28
ASCUSmodel_28 <- glm(HPV_cohortdata_cytfindings_bin28$ASCUS28 ~ HPV_cohortdata_cytfindings_bin28$binomial_ARM, family=binomial(link="log")) 
summary(ASCUSmodel_28)
ci.exp(ASCUSmodel_28, alpha = 0.05)
#LSIL 28
LSILmodel_28 <- glm(HPV_cohortdata_cytfindings_bin28$LSIL28 ~ HPV_cohortdata_cytfindings_bin28$binomial_ARM, family=binomial(link="log")) 
summary(LSILmodel_28)
ci.exp(LSILmodel_28, alpha = 0.05)
#HSIL 28
HSILmodel_28 <- glm(HPV_cohortdata_cytfindings_bin28$HSIL28 ~ HPV_cohortdata_cytfindings_bin28$binomial_ARM, family=binomial(link="log")) 
summary(HSILmodel_28)
ci.exp(HSILmodel_28, alpha = 0.05)


# 4 HAZARD RATIOS AND NON INFERIORITY ----

HPV_cohortdata_allsamplingsA1A2A3_forHazards <- HPV_cohortdata_allages_20240201.yA1A2A3 %>%
  dplyr::mutate(Histologydate = case_when(!is.na(`Date _of_histology.x`) ~ `Date _of_histology.x` ,   
                                          !is.na(`Date _of_histology.y`) ~ `Date _of_histology.y`,  
                                          !is.na(`Date _of_histology`) ~ `Date _of_histology`,
                                          TRUE ~ NA)) %>%
  dplyr::mutate(Histologyvisit = case_when(!is.na(`Date _of_histology.x`) ~ "Hist_v22" ,   
                                           !is.na(`Date _of_histology.y`) ~ "Hist_v25",  
                                           !is.na(`Date _of_histology`) ~ "Hist_v28",
                                           TRUE ~ NA)) %>%
  dplyr::select(FID, ARM_HPV004.x, CIN2.3, CYTvisit22, CYTvisit25, CYTvisit28, `Dates_of_sampling_(cervix)...2.x`, 
                `Dates_of_sampling_(cervix)_25yo`, `Dates_of_sampling_(cervix)_28yo`, Histologydate, Histologyvisit, `COR Positive Genotypes.x`)

ddply(HPV_cohortdata_allsamplingsA1A2A3_forHazards, .(ARM_HPV004.x), summarize, 
      CYT22 = sum(CYTvisit22 == "YES"),
      CYT25 = sum(CYTvisit25 == "YES"),
      CYT28 = sum(CYTvisit28 == "YES")) #CYT visit: the inadequate/missing cyt results are accounted for and the numbers are correct

#Find min and max date and calculate person-time
excluded_col_cin23min <- HPV_cohortdata_allsamplingsA1A2A3_forHazards[, -which(names(HPV_cohortdata_allsamplingsA1A2A3_forHazards) %in% c("FID", "ARM_HPV004.x", "CIN2.3", "CYTvisit22", "CYTvisit25", "CYTvisit28", "Histologydate", "Histologyvisit", "COR Positive Genotypes.x"))]
excluded_col_cin23max <- HPV_cohortdata_allsamplingsA1A2A3_forHazards[, -which(names(HPV_cohortdata_allsamplingsA1A2A3_forHazards) %in% c("FID", "ARM_HPV004.x", "CIN2.3", "CYTvisit22", "CYTvisit25", "CYTvisit28", "Histologyvisit", "COR Positive Genotypes.x"))]
HPV_A1A2A3_forHazards <- HPV_cohortdata_allsamplingsA1A2A3_forHazards %>% 
  dplyr::mutate(Min_date = apply(excluded_col_cin23min, 1, function(x) safe_min(x))) %>%
  dplyr::mutate(Max_date = apply(excluded_col_cin23max, 1, function(x) safe_max(x))) %>%
  dplyr::mutate(Persontime = time_length(difftime(Max_date, Min_date), unit = "year"))

range(HPV_A1A2A3_forHazards$Persontime, na.rm = TRUE)

ddply(HPV_A1A2A3_forHazards, .(ARM_HPV004.x), summarize, 
      Lost_to_follow_up = length(which(Persontime == 0)),
      CIN2plus_infirstvisit = length(which(Persontime == 0 & CIN2.3 == "YES")),
      
      PTime_NA = length(which(is.na(Persontime))),
      
      visit22 = sum(!is.na(`Dates_of_sampling_(cervix)...2.x`)),
      visit25 = sum(!is.na(`Dates_of_sampling_(cervix)_25yo`)),
      visit28 = sum(!is.na(`Dates_of_sampling_(cervix)_28yo`)))

ddply(HPV_A1A2A3_forHazards, .(ARM_HPV004.x), summarize, 
      lost22 = length(which(CYTvisit22=="YES"&CYTvisit25=="NO"&CYTvisit28=="NO")),
      lost25 = length(which(CYTvisit22=="NO"&CYTvisit25=="YES"&CYTvisit28=="NO")),
      lost28 = length(which(CYTvisit22=="NO"&CYTvisit25=="NO"&CYTvisit28=="YES")))

## Cumulative hazards fit (plot and risk table) ----
fit<-survfit(Surv(Persontime,CIN2.3=="YES")~ARM_HPV004.x,data=HPV_A1A2A3_forHazards)
plot(fit,fun="cumhaz")
survplot <- ggsurvplot(fit,data=HPV_A1A2A3_forHazards,fun="cumhaz", risk.table=TRUE, fontsize = 4, conf.int = TRUE) 
survplot_obj <- survplot$plot
survplot_color <- survplot_obj + labs(title="Cumulative Hazard of CIN2+", x="Person-time (years)", y="Cumulative Hazard") +
  scale_color_manual(name = "Arm" , labels = c("A1", "A2", "A3"),
                     values = c( '#800080', '#f0b7ef' ,'#a0e8eb')) +
  scale_fill_manual(name=NULL, labels = NULL, values = c( '#ea9999', '#C3D2BD' ,'#F9CB9C'))+
  guides(fill = FALSE)    +
  theme_minimal() + 
  theme(text = element_text(family = "serif", face=NULL, size=15))
survtable_obj <- survplot$table
survtable_obj <- survtable_obj + scale_y_discrete(name="Arm", labels = c("A3", "A2", "A1")) + scale_x_continuous(name = "Time (years)") +
  labs(title = "Number of women at risk") + theme_classic() + theme(text = element_text(family = "serif", face=NULL, size=15)) 
gridExtra::grid.arrange(survplot_color, survtable_obj)


## Cumulative Hazards with negative for HPV16/18/other at study start ----

#Merge with 18 yo data
data18YO_20250522 <- readxl::read_excel("HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/data_HPV_Vaccination_cohort_18yrsold_Mónica_ver_17.4.2024_kopia.xlsx")
data18YO <- data18YO_20250522 %>% mutate(HPV16.18pos= case_when(HPV_16=="1"~"pos16/18",  HPV_18=="1"~"pos16/18",
                                                       HPV_6=="1"~"otherHPV", HPV_11=="1"~"otherHPV",
                                                       HPV_31=="1"~"otherHPV", HPV_33=="1"~"otherHPV",
                                                       HPV_35=="1"~"otherHPV",HPV_39=="1"~"otherHPV",
                                                       HPV_45=="1"~"otherHPV", HPV_51=="1"~"otherHPV",
                                                       HPV_52=="1"~"otherHPV", HPV_56=="1"~"otherHPV",
                                                       HPV_58=="1"~"otherHPV",HPV_59=="1"~"otherHPV",
                                                       HPV_66=="1"~"otherHPV",
                                                       TRUE~"neg")) %>% dplyr::select(FID,HPV16.18pos)
HPV_A1A2A3_forHazards_plus18yo <- merge(HPV_A1A2A3_forHazards, data18YO, by="FID", all.x = TRUE)

HPV_A1A2A3_forHazards_plus18yo$arm<-Relevel.factor(HPV_A1A2A3_forHazards$ARM_HPV004.x,ref="A1")
HPV_A1A2A3_forHazards_plus18yo$outcome<-if_else(HPV_A1A2A3_forHazards$CIN2.3=="YES",1,
                                       if_else(HPV_A1A2A3_forHazards$CIN2.3=="NO",0,NA))
#A3 split hazards according to HPV negative / HPV16/18 / other HPV
A3_hazards <- HPV_A1A2A3_forHazards_plus18yo %>% dplyr::filter(arm=="A3")
fitA3<-survfit(Surv(Persontime,CIN2.3=="YES")~HPV16.18pos,data=A3_hazards)
survplotA3 <- ggsurvplot(fitA3,data=A3_hazards,fun="cumhaz", con.fint=TRUE, risk.table=TRUE, conf.int = TRUE) 

#(2)Cumulative Hazards with only negative for HPV16/18 at the start of the trial at 18yo
Negatstart_HPV_A1A2A3_forHazards_plus18yo <- HPV_A1A2A3_forHazards_plus18yo %>% dplyr::filter(HPV16.18pos != "pos16/18")
ddply(Negatstart_HPV_A1A2A3_forHazards_plus18yo, .(ARM_HPV004.x), summarize, 
      visit22 = sum(!is.na(`Dates_of_sampling_(cervix)...2.x`)),
      visit25 = sum(!is.na(`Dates_of_sampling_(cervix)_25yo`)),
      visit28 = sum(!is.na(`Dates_of_sampling_(cervix)_28yo`)))
fit2<-survfit(Surv(Persontime,CIN2.3=="YES")~arm,data=Negatstart_HPV_A1A2A3_forHazards_plus18yo)
survplot2 <- ggsurvplot(fit2,data=Negatstart_HPV_A1A2A3_forHazards_plus18yo,fun="cumhaz", risk.table=TRUE, fontsize = 4, conf.int=TRUE)
survplot_obj2 <- survplot2$plot
survplot_color2 <- survplot_obj2 + labs(title="Cumulative hazards of CIN2+", x="Person-time (years)", y="Cumulative Hazard") +
  scale_color_manual(name = "Arm" , labels = c("A1", "A2", "A3"),
                     values = c( '#800080', '#f0b7ef' ,'#a0e8eb')) +
  scale_fill_manual(name=NULL, labels = NULL, values = c( '#ea9999', '#C3D2BD' ,'#F9CB9C'))+
  guides(fill = FALSE)    +
  theme_minimal() + 
  theme(text = element_text(family = "serif", face=NULL, size=15))
survtable_obj2 <- survplot2$table
survtable_obj2 <- survtable_obj2 + scale_y_discrete(name="Arm", labels = c("A3", "A2", "A1")) + scale_x_continuous(name = "Time (years)") +
  labs(title = "Number of women at risk") + theme_classic() + theme(text = element_text(family = "serif", face=NULL, size=15)) 
gridExtra::grid.arrange(survplot_color2, survtable_obj2)

#(3)Cumulative Hazards with negative for HPV16/18/other at the start of the trial at 18yo
Negatstart_HPV_A1A2A3_forHazards_plus18yo <- HPV_A1A2A3_forHazards_plus18yo%>% dplyr::filter(HPV16.18pos != "pos16/18") %>%dplyr::filter(HPV16.18pos != "otherHPV")

ddply(Negatstart_HPV_A1A2A3_forHazards_plus18yo, .(ARM_HPV004.x), summarize, 
      visit22 = sum(!is.na(`Dates_of_sampling_(cervix)...2.x`)),
      visit25 = sum(!is.na(`Dates_of_sampling_(cervix)_25yo`)),
      visit28 = sum(!is.na(`Dates_of_sampling_(cervix)_28yo`)))

fit3<-survfit(Surv(Persontime,CIN2.3=="YES")~arm,data=Negatstart_HPV_A1A2A3_forHazards_plus18yo)
survplot3 <- ggsurvplot(fit3,data=Negatstart_HPV_A1A2A3_forHazards_plus18yo,fun="cumhaz", risk.table=TRUE, fontsize = 4, conf.int=TRUE)
survplot_obj3 <- survplot3$plot
survplot_color3 <- survplot_obj3 + labs(title="Cumulative hazards of CIN2+", x="Person-time (years)", y="Cumulative Hazard") +
  scale_color_manual(name = "Arm" , labels = c("A1", "A2", "A3"),
                     values = c( '#800080', '#f0b7ef' ,'#a0e8eb')) +
  scale_fill_manual(name=NULL, labels = NULL, values = c( '#ea9999', '#C3D2BD' ,'#F9CB9C'))+
  guides(fill = FALSE)    +
  theme_minimal() + 
  theme(text = element_text(family = "serif", face=NULL, size=15))
survtable_obj3 <- survplot3$table
survtable_obj3 <- survtable_obj3 + scale_y_discrete(name="Arm", labels = c("A3", "A2", "A1")) + scale_x_continuous(name = "Time (years)") +
  labs(title = "Number of women at risk") + theme_classic() + theme(text = element_text(family = "serif", face=NULL, size=15)) 
gridExtra::grid.arrange(survplot_color3, survtable_obj3)

##(4)Cumulative Hazards with the negative at both 18-22 for HPV16/18
table(NEG1822_HPV_A1A2A3_forHazards_plus18yo$ARM_HPV004.x)
NEG1822_HPV_A1A2A3_forHazards_plus18yo <- HPV_A1A2A3_forHazards_plus18yo %>%
  dplyr::filter(HPV16.18pos != "pos16/18") %>% 
  dplyr::filter(!grepl("16", `COR Positive Genotypes.x`)) %>% dplyr::filter(!grepl("18", `COR Positive Genotypes.x`))
ddply(NEG1822_HPV_A1A2A3_forHazards_plus18yo, .(ARM_HPV004.x), summarize, 
      visit22 = sum(!is.na(`Dates_of_sampling_(cervix)...2.x`)),
      visit25 = sum(!is.na(`Dates_of_sampling_(cervix)_25yo`)),
      visit28 = sum(!is.na(`Dates_of_sampling_(cervix)_28yo`)))
fitNEG1618<-survfit(Surv(Persontime,CIN2.3=="YES")~arm,data=NEG1822_HPV_A1A2A3_forHazards_plus18yo)
survplotNEG1618 <- ggsurvplot(fitNEG1618,data=NEG1822_HPV_A1A2A3_forHazards_plus18yo,fun="cumhaz", risk.table=TRUE,, fontsize = 4, conf.int=TRUE)
survplot_obj4 <- survplotNEG1618$plot
survplot_color4 <- survplot_obj4 + labs(title="Cumulative hazards of CIN2+", x="Person-time (years)", y="Cumulative Hazard") +
  scale_color_manual(name = "Arm" , labels = c("A1", "A2", "A3"),
                     values = c( '#800080', '#f0b7ef' ,'#a0e8eb')) +
  scale_fill_manual(name=NULL, labels = NULL, values = c( '#ea9999', '#C3D2BD' ,'#F9CB9C'))+
  guides(fill = FALSE)    +
  theme_minimal() + 
  theme(text = element_text(family = "serif", face=NULL, size=15))
survtable_obj4 <- survplotNEG1618$table
survtable_obj4 <- survtable_obj4 + scale_y_discrete(name="Arm", labels = c("A3", "A2", "A1")) + scale_x_continuous(name = "Time (years)") +
  labs(title = "Number of women at risk") +  theme_classic() + theme(text = element_text(family = "serif", face=NULL, size=15)) 
gridExtra::grid.arrange(survplot_color4, survtable_obj4)

##(5)Cumulative Hazards with the negative at 18-22 for HPV16/18/other
NEG_HPV_A1A2A3_forHazards_plus18yo <- HPV_A1A2A3_forHazards_plus18yo %>%
  dplyr::filter(HPV16.18pos != "pos16/18") %>% dplyr::filter(HPV16.18pos != "otherHPV") %>%
  dplyr::filter(is.na(`COR Positive Genotypes.x`)) 
ddply(NEG_HPV_A1A2A3_forHazards_plus18yo, .(ARM_HPV004.x), summarize, 
      visit22 = sum(!is.na(`Dates_of_sampling_(cervix)...2.x`)),
      visit25 = sum(!is.na(`Dates_of_sampling_(cervix)_25yo`)),
      visit28 = sum(!is.na(`Dates_of_sampling_(cervix)_28yo`)))
fitNEG<-survfit(Surv(Persontime,CIN2.3=="YES")~arm,data=NEG_HPV_A1A2A3_forHazards_plus18yo)
survplotNEG <- ggsurvplot(fitNEG,data=NEG_HPV_A1A2A3_forHazards_plus18yo,fun="cumhaz",  risk.table=TRUE, fontsize = 4, conf.int=TRUE)
survplot_obj5 <- survplotNEG$plot
survplot_color5 <- survplot_obj5 + labs(title="Cumulative hazards of CIN2+", x="Person-time (years)", y="Cumulative Hazard") +
  scale_color_manual(name = "Arm" , labels = c("A1", "A2", "A3"),
                     values = c( '#800080', '#f0b7ef' ,'#a0e8eb')) +
  scale_fill_manual(name=NULL, labels = NULL, values = c( '#ea9999', '#C3D2BD' ,'#F9CB9C'))+
  guides(fill = FALSE)    +
  theme_minimal() + 
  theme(text = element_text(family = "serif", face=NULL, size=15))
survtable_obj5 <- survplotNEG$table
survtable_obj5 <- survtable_obj5 + scale_y_discrete(name="Arm", labels = c("A3", "A2", "A1")) + scale_x_continuous(name = "Time (years)") +
  labs(title = "Number of women at risk") + theme_classic() + theme(text = element_text(family = "serif", face=NULL, size=15)) 
gridExtra::grid.arrange(survplot_color5, survtable_obj5)


##Hazard ratios: Complementary log log binomial regression link function (all at start) ----

##change the Hist_time for the 1 out of study 
HPV_A1A2A3_forHazards.1 <-HPV_A1A2A3_forHazards_plus18yo
HPV_A1A2A3_forHazards.1[which(HPV_A1A2A3_forHazards.1$FID == "B31753"),"Histologyvisit"] <- "Hist_v22"

##v25 (attended both visits + 22 and had outcome at 22)
##v28 (attended three visits + 22 and had outcome at 22 + 22and25 and had outcome at 25)
HPV_A1A2A3_forHazards_perprotocol <- HPV_A1A2A3_forHazards.1 %>% 
  dplyr::mutate(perprotocol25=case_when(CYTvisit22=="YES" & Histologyvisit =="Hist_v22" & CIN2.3=="YES" ~"YES",
                                        CYTvisit22=="YES" & CYTvisit25=="YES"~"YES",
                                        TRUE~"NO")) %>%
  dplyr::mutate(perprotocol28=case_when(CYTvisit22=="YES" & Histologyvisit =="Hist_v22"& CIN2.3=="YES"~"YES",
                                        CYTvisit22=="YES" & CYTvisit25=="YES" & Histologyvisit =="Hist_v25" & CIN2.3=="YES" ~"YES",
                                        CYTvisit22=="YES" & CYTvisit25=="YES" & CYTvisit28=="YES"~"YES",
                                        TRUE~"NO")) %>%
  dplyr::mutate(Histvisit = case_when(Histologyvisit == "Hist_v22" ~"H22",
                                      Histologyvisit == "Hist_v25" ~"H25",
                                      Histologyvisit == "Hist_v28" ~"H28",
                                      TRUE~"NO"))

ddply(HPV_A1A2A3_forHazards_perprotocol, .(arm), summarize, 
      v25 = sum(perprotocol25 == "YES"),
      v28 = sum(perprotocol28 == "YES"))

#comparing A1 and A2
HPV_A1A2_forHazards_perprotocol <- HPV_A1A2A3_forHazards_perprotocol %>% dplyr::filter(arm != "A3") %>%
  dplyr::mutate(binomialarm=if_else(arm=="A2", 1, 0)) %>% dplyr::filter(Histvisit!="H22")
protocol25 <- HPV_A1A2_forHazards_perprotocol %>% dplyr::filter(perprotocol25 =="YES") 
modelA1A2_25<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol25)
modelA1A2_25
ci.exp(modelA1A2_25) 
protocol28 <- HPV_A1A2_forHazards_perprotocol %>% dplyr::filter(perprotocol28 =="YES")
modelA1A2_28<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol28)
modelA1A2_28
ci.exp(modelA1A2_28) 

#comparing A1 and A3
HPV_A1A3_forHazards_perprotocol <- HPV_A1A2A3_forHazards_perprotocol %>% dplyr::filter(arm != "A2") %>%
  dplyr::filter(Histvisit!="H22") %>%  dplyr::mutate(binomialarm=if_else(arm=="A3", 1, 0))
protocol25A3 <- HPV_A1A3_forHazards_perprotocol %>% dplyr::filter(perprotocol25 =="YES") 
modelA1A3_25<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol25A3)
modelA1A3_25
ci.exp(modelA1A3_25) 

##Hazard ratios: Complementary log log binomial regression link function (Negative at start) ----

##(1) Neg for HPV 16/18 at 18
HPV_A1A2_forHazards_perprotocol <- HPV_A1A2A3_forHazards_perprotocol %>% dplyr::filter(arm != "A3") %>% dplyr::filter(HPV16.18pos!="pos16/18") %>%
  dplyr::filter(Histvisit!="H22") %>% dplyr::mutate(binomialarm=if_else(arm=="A2", 1, 0))
protocol25 <- HPV_A1A2_forHazards_perprotocol %>% dplyr::filter(perprotocol25 =="YES") 
modelA1A2_25<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol25)
modelA1A2_25
ci.exp(modelA1A2_25) 

protocol28 <- HPV_A1A2_forHazards_perprotocol %>% dplyr::filter(perprotocol28 =="YES")
modelA1A2_28<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol28)
modelA1A2_28
ci.exp(modelA1A2_28) 

HPV_A1A3_forHazards_perprotocol <- HPV_A1A2A3_forHazards_perprotocol %>% 
  dplyr::filter(arm != "A2") %>% 
  dplyr::filter(HPV16.18pos!="pos16/18") %>%
  dplyr::filter(Histvisit!="H22")%>%  
  dplyr::mutate(binomialarm=if_else(arm=="A3", 1, 0))

protocol25A3 <- HPV_A1A3_forHazards_perprotocol %>% dplyr::filter(perprotocol25 =="YES") 
modelA1A3_25<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol25A3)
modelA1A3_25
ci.exp(modelA1A3_25)

##(2) Neg for HPV 16/18/other at 18
HPV_A1A2_forHazards_perprotocol <- HPV_A1A2A3_forHazards_perprotocol %>% dplyr::filter(arm != "A3") %>% dplyr::filter(HPV16.18pos=="neg") %>%
  dplyr::filter(Histvisit!="H22") %>% dplyr::mutate(binomialarm=if_else(arm=="A2", 1, 0))
protocol25 <- HPV_A1A2_forHazards_perprotocol %>% dplyr::filter(perprotocol25 =="YES") 
modelA1A2_25<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol25)
modelA1A2_25
ci.exp(modelA1A2_25) 

protocol28 <- HPV_A1A2_forHazards_perprotocol %>% dplyr::filter(perprotocol28 =="YES")
modelA1A2_28<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol28)
modelA1A2_28
ci.exp(modelA1A2_28) 

HPV_A1A3_forHazards_perprotocol <- HPV_A1A2A3_forHazards_perprotocol %>% dplyr::filter(arm != "A2") %>% dplyr::filter(HPV16.18pos!="neg") %>%
  dplyr::filter(Histvisit!="H22")%>%  dplyr::mutate(binomialarm=if_else(arm=="A3", 1, 0))
protocol25A3 <- HPV_A1A3_forHazards_perprotocol %>% dplyr::filter(perprotocol25 =="YES") 
modelA1A3_25<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol25A3)
modelA1A3_25
ci.exp(modelA1A3_25)

##(3) Neg for HPV 16/18 at 18+22
HPV_A1A2_forHazards_perprotocol <- HPV_A1A2A3_forHazards_perprotocol %>% dplyr::filter(arm != "A3") %>% dplyr::filter(HPV16.18pos != "pos16/18") %>% 
  dplyr::filter(Histvisit!="H22") %>% dplyr::mutate(binomialarm=if_else(arm=="A2", 1, 0)) %>% dplyr::filter(!grepl("16", `COR Positive Genotypes.x`)) %>% dplyr::filter(!grepl("18", `COR Positive Genotypes.x`))
protocol25 <- HPV_A1A2_forHazards_perprotocol %>% dplyr::filter(perprotocol25 =="YES") 
modelA1A2_25<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol25)
modelA1A2_25
ci.exp(modelA1A2_25) 

protocol28 <- HPV_A1A2_forHazards_perprotocol %>% dplyr::filter(perprotocol28 =="YES")
modelA1A2_28<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol28)
modelA1A2_28
ci.exp(modelA1A2_28) 

HPV_A1A3_forHazards_perprotocol <- HPV_A1A2A3_forHazards_perprotocol %>% dplyr::filter(arm != "A2") %>% dplyr::filter(HPV16.18pos != "pos16/18") %>% 
  dplyr::filter(Histvisit!="H22") %>% dplyr::mutate(binomialarm=if_else(arm=="A3", 1, 0)) %>% dplyr::filter(!grepl("16", `COR Positive Genotypes.x`)) %>% dplyr::filter(!grepl("18", `COR Positive Genotypes.x`))
protocol25A3 <- HPV_A1A3_forHazards_perprotocol %>% dplyr::filter(perprotocol25 =="YES") 
modelA1A3_25<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol25A3)
modelA1A3_25
ci.exp(modelA1A3_25)

##(4) Neg for HPV 16/18/other at 18+22
HPV_A1A2_forHazards_perprotocol <- HPV_A1A2A3_forHazards_perprotocol %>% dplyr::filter(arm != "A3") %>% dplyr::filter(HPV16.18pos != "pos16/18") %>% dplyr::filter(HPV16.18pos != "otherHPV") %>%
  dplyr::filter(Histvisit!="H22") %>% dplyr::mutate(binomialarm=if_else(arm=="A2", 1, 0)) %>% dplyr::filter(is.na(`COR Positive Genotypes.x`)) 
protocol25 <- HPV_A1A2_forHazards_perprotocol %>% dplyr::filter(perprotocol25 =="YES") 
modelA1A2_25<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol25)
modelA1A2_25
ci.exp(modelA1A2_25) 

protocol28 <- HPV_A1A2_forHazards_perprotocol %>% dplyr::filter(perprotocol28 =="YES")
modelA1A2_28<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol28)
modelA1A2_28
ci.exp(modelA1A2_28) 

HPV_A1A3_forHazards_perprotocol <- HPV_A1A2A3_forHazards_perprotocol %>% dplyr::filter(arm != "A2") %>% dplyr::filter(HPV16.18pos != "pos16/18") %>% dplyr::filter(HPV16.18pos != "otherHPV") %>%
  dplyr::filter(Histvisit!="H22") %>% dplyr::mutate(binomialarm=if_else(arm=="A3", 1, 0)) %>% dplyr::filter(is.na(`COR Positive Genotypes.x`)) 
protocol25A3 <- HPV_A1A3_forHazards_perprotocol %>% dplyr::filter(perprotocol25 =="YES") 
modelA1A3_25<-glm(outcome~binomialarm, family=binomial(link="cloglog"), data=protocol25A3)
modelA1A3_25
ci.exp(modelA1A3_25)


# 5 QUESTIONAIRE DATA: Baseline characteristics ----
#Load questionnaire data (age of 22)
Questionnairedata22_16022024 <- readxl::read_excel("HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/BQ to Monica and Penny 16.2.2024 - kopia.xlsx")
Questionnairedata22_16022024 <- Questionnairedata22_16022024[-1,] #eliminate first row with names

#Age of sexual debut                                ##OBS: There's many non-numerical answers, make a new variable for adjusting, same with smoking and alcohol intoxications
Questionnairedata22_16022024.x <- Questionnairedata22_16022024 %>% mutate(sexual.debut = case_when(`kys23#` == "13-14" ~ "13.5", `kys23#` == "13714"~ "13.5",
                                                                                                   `kys23#` == "14?"~ "14",
                                                                                                   `kys23#` == "14-15"~ "14.5",
                                                                                                   `kys23#` == "15-16"~ "15.5",
                                                                                                   `kys23#` == "16-17"~ "16.5", `kys23#` == "16 tai 17"~ "16.5", `kys23#` == "16/17"~ "16.5",
                                                                                                   `kys23#` == "16?"~ "16",
                                                                                                   `kys23#` == "18-19"~ "18.5", `kys23#` == "18Â½"~ "18.5",
                                                                                                   `kys23#` == "19-20"~ "19.5",
                                                                                                   TRUE ~ `kys23#`)) %>%
  mutate(smoke.habit = case_when(`kys13#` == "0" ~ "Missing?",
                                 `kys13#` == "1" ~ "Never",
                                 `kys13#` == "2" ~ "Quiter",
                                 `kys13#` == "3" ~ "Smoker",
                                 `kys13#` == "4" ~ "Other than cigarettes")) %>%
  mutate(alcohol.intoxication = case_when(`kys17#` == "0" ~ "Missing?",
                                          `kys17#` == "1" ~ "Never",
                                          `kys17#` == "2" ~ "<Once a month",
                                          `kys17#` == "3" ~ "1/2 a month",
                                          `kys17#` == "4" ~ ">Once a week"))

Questionnairedata22_16022024.x$sexual.debut <- as.numeric(Questionnairedata22_16022024.x$sexual.debut)
sum(is.na(Questionnairedata22_16022024.x$`kys23#`))
sum(is.na(Questionnairedata22_16022024.x$sexual.debut))

#no. of lifetime sexual partners
sumdata <- as.data.frame(table(Questionnairedata22_16022024$kys26)) #no need to change anything 
sum(sumdata$Freq)
as.data.frame(table(Questionnairedata22_16022024$`kys22#`))

#smoking
as.data.frame(table(Questionnairedata22_16022024$`kys13#`)) 
#alcohol use
as.data.frame(table(Questionnairedata22_16022024$`kys17#`))
#condom use
as.data.frame(table(Questionnairedata22_16022024$kys59_1)) #When do you use condom: always
as.data.frame(table(Questionnairedata22_16022024$kys59_2)) #beginning of a sexual relation
as.data.frame(table(Questionnairedata22_16022024$kys59_3)) #random sexual relationship
as.data.frame(table(Questionnairedata22_16022024$kys59_4)) #after forgetting the pill

#summaries
ddply(Questionnairedata22_16022024.x, .(ARM_HPV004), summarise, #for age of sexual debut, mean and sd
      n = length(ARM_HPV004),
      Age.sexual.debut.mean = mean(sexual.debut, na.rm=TRUE),
      Age.sexual.debut.sd = sd(sexual.debut, na.rm=TRUE))

ddply(Questionnairedata22_16022024.x, .(ARM_HPV004, `kys22#`), summarise, #for no. of lifetime partners (OBS:SEX no = 0 lifetime sexual partners // if sex yes = then how many partners)
      n = length(ARM_HPV004),
      sex.NO = sum(`kys22#` == "1"),
      sex.YES = sum(`kys22#` == "2"),
      partners.missing = sum(kys26 == "0"),
      partners.0 = sum(kys26 == "1"),
      partners.1 = sum(kys26 == "2"),
      partners.2 = sum(kys26 == "3"),
      partners.3 = sum(kys26 == "4"),
      partners.4 = sum(kys26 == "5"),
      partners.5plus = sum(kys26 == "6")) %>% dplyr::filter(`kys22#` == "2") %>% mutate(n=c(2255, 2372,1047,1)) %>%
  mutate(part1.pr = partners.1/n) %>% mutate(part2.pr = partners.2/n) %>% mutate(part3.pr = partners.3/n) %>%
  mutate(part4.pr = partners.4/n) %>% mutate(part5pl.pr = partners.5plus/n) 

misA1 <- 2255-(93+369+256+272+256+950) #the actual missing answers for each arm 
misA2 <- 2372-(90+371+278+284+295+1001)
misA3 <- 1047-(38+200+132+121+125+397)
misA1/2255
misA2/2372
misA3/1047

ddply(Questionnairedata22_16022024.x, .(ARM_HPV004), summarise, #smoking habits 
      n = length(ARM_HPV004),
      Smoking.never = sum(smoke.habit == "Never"),
      Smoking.Quiter = sum(smoke.habit == "Quiter"),
      Smoking.yes = sum(smoke.habit == "Smoker"),
      Other.than.cig = sum(smoke.habit == "Other than cigarettes" ),
      Smoking.missing = sum(smoke.habit == "Missing?")) %>% mutate(never.pr = Smoking.never/n) %>% mutate(quiter.pr = Smoking.Quiter/n) %>%
  mutate(yes.pr = Smoking.yes/n) %>% mutate(other.pr = Other.than.cig/n) %>% mutate(missing.pr = Smoking.missing/n)

ddply(Questionnairedata22_16022024.x, .(ARM_HPV004), summarise, #For alcohol intoxication
      n = length(ARM_HPV004),
      Intox.never = sum(alcohol.intoxication == "Never"),
      Intox.1month = sum(alcohol.intoxication == "<Once a month"),
      Intox.12month = sum(alcohol.intoxication == "1/2 a month"),
      Intox.week = sum(alcohol.intoxication == ">Once a week"),
      Intox.missing = sum(alcohol.intoxication == "Missing?")) %>% mutate(never.pr = Intox.never/n) %>% mutate(month1.p = Intox.1month/n) %>%
  mutate(month12.pr = Intox.12month/n) %>% mutate(week.pr = Intox.week/n) %>% mutate(miss.pr = Intox.missing/n)

ddply(Questionnairedata22_16022024.x, .(ARM_HPV004), summarise, 
      n = length(ARM_HPV004),
      C.always = sum(kys59_1 == "1"),
      C.beginningofrel = sum(kys59_2 == "1"),
      C.randomsex = sum(kys59_3 == "1"),
      C.nopill = sum(kys59_4 == "1"),
      C.never = sum(kys60 == "1"),
      C.missing = sum(kys59_1 == "0"& kys59_2 == "0" & kys59_3 == "0" & kys59_4 == "0"& kys59_5 == "0" & kys60 == "0")) %>% 
  mutate(Cnever.pr = C.never/n) %>% mutate(C.always.pr = C.always/n) %>%
  mutate(C.beginningofrel.pr = C.beginningofrel/n) %>% mutate(C.randomsex.pr = C.randomsex/n) %>% mutate(C.nopill.pr = C.nopill/n)

## Date of vaccination ----
Vaccinationdatesandages22 <- HPV_cohortdata_22yo_newcyt %>% 
  dplyr::mutate(A1A2vaccyear= format(dose_1, "%Y")) %>% 
  dplyr::mutate(A3vaccyear= format(I_dose, "%Y")) %>% 
  dplyr::select(FID, ARM_HPV004, A1A2vaccyear, A3vaccyear, Birth_year) 
Vaccinationdatesandages22$A1A2vaccyear <- as.numeric(Vaccinationdatesandages22$A1A2vaccyear)
Vaccinationdatesandages22$A3vaccyear <- as.numeric(Vaccinationdatesandages22$A3vaccyear)

Vaccinationdatesandages22 <- Vaccinationdatesandages22 %>% dplyr::mutate(age_at_vaccination = case_when(ARM_HPV004 =="A1" ~ A1A2vaccyear-Birth_year , 
                                               ARM_HPV004 =="A2" ~ A1A2vaccyear-Birth_year , 
                                               ARM_HPV004 =="A3" ~ A3vaccyear-Birth_year ))

ddply(Vaccinationdatesandages22, .(ARM_HPV004), summarise, 
             n = length(ARM_HPV004),
             Age.vac.mean = mean(age_at_vaccination, na.rm=TRUE),
             Age.vac.sd = sd(age_at_vaccination, na.rm=TRUE))


# 6 EXCLUSION/LOST TO FOLLOWUP FOR LEEP CIRURGY ----

#Load data on LEEP (individuals that got leep with dates)
LEEPdata_1582024 <- readxl::read_excel("HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/data/To Monica leep_15.8.2024_kopia.xlsx")
arminfo <- HPV_cohortdata_22yo_newcyt %>% dplyr::select(FID, ARM_HPV004)

LEEPdata_witharm <- merge(LEEPdata_1582024, arminfo, all.x = TRUE, by= "FID")
ddply(LEEPdata_witharm, .(ARM_HPV004, WHEN_LEEP), summarize,
      n = length(FID))

#lost to followup - attendance
CYT22yes <- HPV_cohortdata_allages_20240201.yA1A2A3 %>% dplyr::filter(grepl("Female|female", Gender.x)) %>% dplyr::mutate(list22="YES")%>% select(FID, ARM_HPV004.x, list22) 
CYT25yes <- HPV_cohortdata_allages_20240201.yA1A2A3 %>% dplyr::filter(grepl("Female|female", Gender.y)) %>% dplyr::mutate(list25="YES") %>% select(FID, list25)
CYT28yes <- HPV_cohortdata_allages_20240201.yA1A2A3 %>% dplyr::filter(grepl("Female|female", Gender)) %>% dplyr::mutate(list28="YES") %>% select(FID,list28)
ddply(CYT22yes, .(ARM_HPV004.x), summarize,
      n = length(FID))

CYT22to25 <- merge(CYT22yes, CYT25yes, by="FID", all.x = TRUE )
CYT22to25to28 <- merge(CYT22to25, CYT28yes, by="FID", all.x = TRUE )
CYT22to25to28<- CYT22to25to28 %>% dplyr::mutate(Losttofollowup_CYT = case_when(list22 =="YES" & is.na(list25)==TRUE & is.na(list28)==TRUE ~ "No25No28",
                                                                               list22 =="YES" & is.na(list25)==TRUE & list28 =="YES" ~ "Miss25",
                                                                               list22 =="YES" & list25 =="YES" & is.na(list28)==TRUE ~ "Miss28",
                                                                               list22 =="YES" & list25 =="YES"& list28 =="YES" ~ "All_visits"))
ddply(CYT22to25to28, .(ARM_HPV004.x), summarize,
      All_visits = sum(Losttofollowup_CYT=="All_visits"),     
      Miss25 =  sum(Losttofollowup_CYT=="Miss25"),      
      Miss28 =  sum(Losttofollowup_CYT=="Miss28"),   
      No25No28=sum(Losttofollowup_CYT=="No25No28"))

                                 
#DEPENDENCIES, PACKAGES AND VERSIONS-----
packINFO <- as.data.frame(installed.packages())[,c("Package", "Version")]
rownames(packINFO) <- NULL
dependencies.pkg <- dependencies("C:/Users/monort/OneDrive - Karolinska Institutet/Dokument/HPV004-Effectiveness of Cervical Screening in HPV vaccinated women/RCodes_004/004_CLEANCODE_FREQUENTVSINFREQUENT_20250220.R")
dependencies.pkg <- as.data.frame(dependencies.pkg$Package)
colnames(dependencies.pkg) <- c("Package")
depencencies.versions <- merge(dependencies.pkg,packINFO, by="Package", all.x = TRUE )

#----
