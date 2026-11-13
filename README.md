# Wild boar habitat suitability drove African swine fever transmission dynamics in South Korea

This repository contains the official R code and datasets for the study:\
**"Wild boar habitat suitability drove African swine fever transmission dynamics in South Korea"**

This project is archived on Zenodo to ensure long-term accessibility and reproducibility of the research findings.

## Project Structure

Based on the repository organization, the folders are categorized as follows:

- **`model_code/`**: R scripts for the candidate transmission models, including model fitting, simulation, and post-hoc analysis.
- **`model_code/legacy/`**: Archived model scripts retained for reference.
- **`simulation_function/`**: Shared R functions for simulation, model initialization, summary statistics, and supporting calculations.
- **`data/`**: Datasets including habitat suitability indices and ASF transmission records. (Access Restricted)
- **`results/`**: Model-specific outputs. Create this directory before running the scripts. This directory is excluded from version control.

## Getting Started

### Prerequisites

- **R (version 4.0 or higher)** is required.
- This project primarily utilizes the **Approximate Bayesian Computation (ABC)** framework.

### Running a Model

Before running a model, replace `/path/to/github_code` in the R scripts under `model_code/` and `simulation_function/` with the absolute path to your copy of this repository. Create a `results/` directory in the repository root before running.

Then, from the repository root, run:

```sh
Rscript model_code/Model_A.R
```

Alternatively, in R:

```r
source("/path/to/github_code/model_code/Model_A.R")
```

The scripts use explicit file paths to access `simulation_function/` and `data/`. Outputs are saved under `results/` in a subdirectory named after the model.

## License

This project applies different licenses for the source code and the research data.

### Source Code

The source code in this repository is licensed under the **MIT License**.

- You are free to use, modify, and distribute the code as long as the original license and copyright notice are included.

### Research Data

**Copyright (c) 2025 Epidesa. All Rights Reserved.**

- The data provided in this repository (typically located in the `data/` directory) is for **archival and peer-review/validation purposes only**.

### Caution

**No unauthorized use, reproduction, redistribution, or modification** of the data is permitted without explicit written consent from the authors and the Korean government.
