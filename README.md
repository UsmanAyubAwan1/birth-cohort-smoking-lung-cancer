# Smoking Prevalence by Birth Cohort, Cigarettes per Smoker, and Lung Cancer Mortality, 2000-2022

This repository contains the R code for the statistical analysis reported in the
manuscript of the same title. The study is a longitudinal ecological analysis of
lung cancer mortality in 64 countries from 2000 to 2022, with lung cancer incidence
in 43 countries as a secondary outcome.

## Requirements

- R 4.4 or later
- R packages `readxl` and `MASS`

```r
install.packages("readxl")   # MASS is included with R
```

## Data

The script reads `JAMA_Replication_Dataset.xlsx` from the working directory and uses
three of its sheets: `Mortality`, `Control_Outcomes`, and `CI5_Incidence`.

The dataset was assembled from publicly available sources:

- WHO Mortality Database (deaths by cause, sex, and age group)
- United Nations World Population Prospects 2024 and United Nations Statistics
  Division (population estimates and death registration coverage)
- Institute for Health Metrics and Evaluation, Global Burden of Disease 2023
  (smoking prevalence, cigarettes per smoker per day, and Sociodemographic Index)
- World Bank (country classification and fine particulate matter)
- International Agency for Research on Cancer, Cancer Incidence in Five
  Continents (CI5plus)

Details of each source are given in the Methods section of the manuscript. The
assembled analysis dataset is available from the corresponding authors on request.

## Running the analysis

1. Put `Final Code.R` and `JAMA_Replication_Dataset.xlsx` in the same folder.
2. Set that folder as the working directory, for example with `setwd()` or the
   RStudio menu Session > Set Working Directory.
3. Run the script:

```r
source("Final Code.R")
```

Each set of estimates is printed to the console as it is computed, and all
estimates are saved to `analysis_results.csv`. A full run takes about 1 hour on a
desktop computer.

## Contents of the script

| Section | Analysis | Manuscript items |
|---|---|---|
| 1 | Data and model set-up | |
| 2 | Mortality trends within countries | Figure 1, eTable 7 |
| 3 | Primary associations, absolute rate differences, likelihood ratio test | Table 2, eTables 8 and 12 |
| 4 | Effect modification with false discovery rate control | Table 2, eTable 9 |
| 5 | Birth-cohort-aligned prevalence | Table 2, Figure 2, eTables 14 and 15 |
| 6 | Lung cancer incidence in CI5plus registries | eTable 17 |
| 7 | Control causes of death | Figure 3A, eTable 18 |
| 8 | Histologic types | Figure 3B, eTables 19 and 20 |
| 9 | Sensitivity analyses | eTable 22 |

All models are negative binomial (NB2) regressions with a log link and the
logarithm of population (or person-years) as an offset. Confidence intervals use
country-clustered sandwich variances with the CR1 correction and a t distribution
with (number of countries - 1) degrees of freedom.

## License

The code is released under the MIT License. See `LICENSE`.
