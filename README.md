# MLfactor

This repository contains the code and resources for a simulation study designed to demonstrate under which circumstances violations of measurement invariance can produce item-level predictive gains. The code is structured to enable reproducibility and transparency of the results.

---

## Table of Contents
1. [Overview](#overview)
2. [Folder Structure](#folder-structure)
3. [Getting Started](#getting-started)
4. [File Descriptions](#file-descriptions)
5. [Usage](#usage)
6. [Contact Information](#contact-information)

---

## Overview

In personality psychology, prediction models that use individual items as predictors (item-level models) often outperform models that use aggregate scale scores as predictors (e.g., sum scores or factor scores representing facets or trait domains). This pattern is commonly attributed to stable, outcome-relevant item-specific variance (“nuances”). Using a Monte Carlo simulation, we demonstrate that violations of measurement invariance can also result in item-level predictive gains even when no substantive item-specific effects are present in the data-generating process. We generated item responses under reflective factor models and manipulated key aspects of the measurement model (metric vs. scalar non-invariance, proportion and source of item-specific variance, number of affected items), the explained variance in the outcome, and sample size. Predictive performance was evaluated out of sample using regularized regression (elastic net regression; with ridge and OLS regression for comparison) and quantified as difference in $R^2$ between item-level and score-level prediction models ($\Delta R^2$). Under specific conditions, MI violations produced item-level predictive gains of similar magnitude as induced by direct item-specific effects. In particular, scalar non-invariance yielded gains when the grouping variable that generated non-invariance was correlated with the outcome, whereas metric non-invariance yielded gains mainly when the latent variable was strongly related to the outcome. These results suggest that, under certain conditions, apparent item-level predictive gains can also arise from group-dependent measurement and need not reflect nuances. We discuss implications for interpreting item-level predictive gains and for evaluating predictive models across levels of the personality hierarchy.

---

## Folder Structure

project/\
  ├── 01_simulateData.R # Script to simulate raw data\
  ├── 02_fitData.R # Script to fit models and save results\
  ├── 03_joinData.R # Script to merge data across conditions\
  ├── 04a_analysePerformance.R # runs ANOVA\
  ├── 04b_plotPerformance_ViolinPlots.R # plot performance plots\
  ├── MLfactor.Rproj # R project\
  ├── README.md # Documentation (this file)\
  ├── utils/ # Utility functions for simulation and analysis\
  │ └── [utility scripts, e.g., fitOLS.R, fitPenalizedReg.R, ..., setParam.R]\
  ├── onlineMaterial/ # code and plots for supplementary analyses, etc.\
  ├── results/ # Folder for fitted models and dependent measure files\
  │ ├── dependentMeasures/ # Contains pre-processed result files\
  ├── plots/ # Folder for result plots\
  │ ├── ANOVA/ # Contains plots for the ANOVA results\
  │ ├── checkGraphics/ # Contains most of the Online Material plots\
  │ ├── hyperPlots/ # Contains histograms for hyperparameter choices\
  │ ├── performance/ # Contains Delta R² plots\
  │ ├── ... paper png-files # Contains plots that made it into the paper\
  ├── data/ # Folder containing simulated data\
  │   └── [simulated data files; this data is not provided due to memory limitations]\
  └── log/ # log files from data simulation and fitting\

---

## Getting Started

1. **Dependencies**:
    - R (>= 4.5)
    - Required R packages: `mvtnorm`, `parallel`, `lavaan`, `semTools`, etc.

2. **Installation**:
    Clone this repository:
    ```bash
    git clone https://github.com/kimSpeck/MLfactor.git
    ```

3. **Usage**:
    Run the scripts in the specified order:
    1. Simulate data: 01_simulateData.R
    2. Fit models to the simulated data: 02_fitData.R
    3. Merge data: 03_joinData.R
    4. Analyse data: {04...} files

---
## File Descriptions

### Main Scripts

1. **01_simulateData.R**:
   - **Purpose**: Simulates raw data using a multigroup SEM, sampling from a multivariate normal distribution.
   - **Key Features**:
     - uses utility functions from the `utils/` folder (e.g., `utils/setParam.R`, ...)
     - Outputs raw data stored in the `data/` folder.

2. **02_fitData.R**:
   - **Purpose**: Fits models to the simulated data and stores results for each simulated condition.
   - **Output**: 
     - uses utility functions from the `utils/` folder (e.g., `fitOLS.R`, `fitPenalizedReg.R`, ...) 
     - Results saved as `.rds` files in the `results/` folder.

3. **03_joinData.R**:
   - **Purpose**: Merges data from all simulated conditions for analysis.
   - **Features**:
     - Creates subfiles for dependent measures to avoid RAM overflow when analysing results.
     - Outputs merged data files stored in the `results/dependentMeasures` folder.

4. **04a_analysePerformance.R**:
   - **Purpose**: Analyzes \(\Delta R^2\) results to determine the effect of experimental manipulations.
   - **Process**:
     - Performs ANOVA to assess model performance and plots results of the generalized $\eta^2$. 
     - Prepares data for visualizations and generates ANOVA results plots. 
     
5. **04b_plotPerformance_ViolinPlots.R**:
    - **Purpose**: Plot \(\Delta R^2\) result graphics
    - **Process**:
      - Plots results graphics (paper-style).
      - Plots results graphics for the full simulation design (provided in online material).

---

### Utility Files

- Located in the `utils/` folder.
- Files and functions to facilitate...
  - ... overviewing simulation conditions and parameter setup (`setParam.R`).
  - ... simulation of the data and analysis of the results (`simTools.R`, `utilityFUNS.R`)
  - ... model fitting (`fitOLS.R`, `fitPenalizedReg.R`, `fitTools.R`)

---

### Data Folder

- The `data/` folder contains raw simulated data produced by `01_simulateData.R`.

---

### Results Folder

- The `results/` folder stores:
  - Fitted model outputs.
  - Subfiles with dependent measures.

---

## Usage

Follow the script order for reproducibility:
1. Simulate data using `01_simulateData.R`.
2. Fit models with `02_fitData.R`.
3. Merge and analyze data:
   - Run `03_joinData.R`.
   - Use the {04...} scripts to replicate the results.

---

## Contact Information

If you have any questions or encounter issues, please contact:
- **Name**: Kim-Laura Speck
- **Email**: kim.speck@uni-kassel.de
- **Affiliation**: Universität Kassel

---