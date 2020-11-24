## Table of contents
* [General Info](#general-info)
* [Functional Data Analysis](#introduction)
* [Data](#data)

## General Info

### Contents:
* **main.Rmd** is the main RMarkdown file where all the analysis is conducted.
* The Figures folder contains all the figures that are produced in main.Rmd.
* The Data folder contains the Population and Mortality data from the Central Agency for Statistics in the Netherlands, as well as metadata.
* **Report.pdf** is a report in the form of a scientific article, which includes a detailed explanation of FDA and all methods that were applied in this project. Additionally, it contains an exploratory analysis, the modeling proccess and references. 

### Libraries:
- gridExtra 2.3
- ftsa 5.8
- tidyr 1.1.0
- magrittr 1.5
- rainbow 3.6
- fda 5.1.4
- tidyverse 1.3.0
- dplyr 1.0.0
- fda.usc 2.0.2


## Functional Data Analysis
Functional Data Analysis (FDA) has received increased attention within the demographic forecasting framework since it was first introduced by Hyndman and Ullah (2007) as a robust approach on forecasting mortality and fertility rates. In this project we apply this methodology on the Dutch mortality rates for the years 1950-2018. We make a model comparison by fitting Functional Principal Components Analysis (FPCA), weighted FPCA
and Functional Partial Least Squares models on the age-specific mortality rates from 1950 to 1999, for the male and female populations separately. We obtain the results of the evaluation by calculating the Root Mean Squared Error (RMSE) between the forecasted and true mortality rates from 2000 to 2018. They show that for the male population, the weighted FPCA outperformed the other models, while for the female population the FPCA and weighted FPCA models had a very similar performance.


## Data
In this project, we use population and mortality data by age and sex for the Netherlands. The data were collected from the online database of the Central Agency for Statistics in the Netherlands
(Centraal Bureau voor de Statistiek, or CBS). More specifically, two datasets were combined in order to derive the mortality rates by age and sex. These include the population data by age and sex, and the total deaths by age and sex. According to CBS, the population data become available at the 1st of January, whereas the total number of deaths, are being published on December of the same year. The two databases consist of population and mortality data for males and females
of all ages across a 68 year period, starting from 1950 up to 2018. Its also important to note that data are recorded for 21, 5-year age groups (0, 1-5, 5-10, 10-15, ..., 90-95, 95+) with the only exceptions being the newborns and the age group of 95+ years old due to their distinctiveness on displaying higher mortality rates. Therefore, the mortality rates are calculated by dividing the number of deaths for a particular age group within a calendar year, by the total number of individuals in this age group for the same year.

