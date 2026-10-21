# Cancer Outpatient 'Did Not Attend' (DNA) Risk Predictor Pipeline

This repository hosts the scripts necessary to execute the model training and prediction pipeline for the outpatient DNA (Did Not Attend) risk tool.

## Installation

Clone the repository and navigate into the project directory:

```bash
git clone git@github.com:niklz/dna_predictor.git
```

*(Note: There is no project environment for this project at this stage. A list of required packages can be found in the `00_libraries_and_utils.R` script, and this code will run using R version 4.5.1).*

## Training Data

Data is ingested in the `01_data_prep.R` script. Due to information governance, patient data is not included in the repository. 

To run this pipeline locally, you must provide your own dataset and reference the correct file path on **line 11** of `01_data_prep.R`.

## Usage

The analytical process is broken into enumerated modules (`00_..R`, `01_..R`, etc.) and executed via two primary orchestration scripts. 

You can run these from the command line:

**1. Train the model** (Trains and caches the model)
```bash
Rscript master_train.R
```

**2. Generate predictions** (Generates predictions on new data)
```bash
Rscript master_predict.R
```

---

## Contact
**Author:** Nick Howlett, nick.howlett5@nhs.net