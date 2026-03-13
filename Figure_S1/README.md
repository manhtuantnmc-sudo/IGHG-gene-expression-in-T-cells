
# README – RNA-seq Analysis Code (FigureS1)

## Overview
This repository contains additional R scripts used to generate figures from the same RNA-seq dataset obtained from the ImmPort database.

The data used in this study are available from ImmPort under the study accession: **SDY998**.

Because the dataset is publicly available, the raw data are **not included in this repository**.  
Please download the data directly from ImmPort and run the analysis scripts provided here.

---

## Data Source

**Database:** ImmPort (Immunology Database and Analysis Portal)  
**Study accession:** SDY998  

Study link:  
https://www.immport.org/shared/study/SDY998

The following files from the study were used in this analysis:

- Aggregated count matrix  
- Sample metadata  

---

## Analysis Order

These scripts use the **same dataset** as the Figure 1 analysis but generate **different figures**.

To avoid unnecessary duplication of code, the preprocessing and initial analysis steps are **not repeated here**.

Therefore:

1. First run the analysis scripts used for **Figure 1**.
2. After completing the Figure 1 analysis, run the scripts provided in this directory.

These scripts assume that the objects or processed data generated in the Figure 1 analysis are already available.

---

## Reproducibility

To reproduce the analysis:

1. Download the dataset from ImmPort (SDY998)  
2. Run the analysis scripts used for **Figure 1**  
3. Then run the scripts provided in this directory  

The scripts will generate the additional figures reported in the manuscript without repeating the earlier preprocessing steps.
