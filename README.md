# Finnish screening trial

This repository compiles the R code used to analyze the data from the Finnish screening trial NCT02149030 aiming to assess the effectiveness of cervical screening frequency among women HPV vaccinated as early adolescents.

### The code is divided in different sections:
First data loading (participant information, HPV DNA results and Cytology results) and merging; HPV DNA and cytlological findings per screening visit; prevalence of CIN2+ calculations and histology cases; Hazard calculations and finally Analysis of questionnaire-based data (baseline characteristics and risk behaviour in the different study arms). 


For this project RStudio 2024.04.1+748 version was used.

## Package dependencies and verison:

caret	6.0-94
data.table	1.15.4
dplyr	1.1.4
Epi	2.51
ggplot2	3.5.1
ggpubr	0.6.0
grid	4.4.0
gridExtra	2.3
gtable	0.3.5
lubridate	1.9.3
plyr	1.8.9
purrr	1.0.2
readr	2.1.5
readxl	1.4.3
rlist	0.4.6.2
stringr	1.5.1
survival	3.8-3
survminer	0.5.0
systemfonts	1.1.0
tibble	3.2.1
tidyr	1.3.1
tidyverse	2.0.0

