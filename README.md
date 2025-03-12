[![DOI](https://zenodo.org/badge/756931849.svg)](https://zenodo.org/doi/10.5281/zenodo.11067889)
[![DOI](https://zenodo.org/badge/756931849.svg)](https://zenodo.org/doi/10.5281/zenodo.15011335)

# Epidemicity indices and reproduction numbers from infectious disease data in connected human populations

This repository contains code and data described in

C. Trevisin, L. Mari, M. Gatto, A. Rinaldo.
Epidemicity indices and reproduction numbers from infectious disease data in connected human populations
*Infectious Disease Modelling* (2024), **9**(3), pp. 875-891

## Repository structure

- The top-level directory contains MATLAB® scripts for running all the analyses described in the above paper and generating the relevant figures:

	- `experiment_synthetic1.m` (generates Figures 1 and 2 of the above paper)
	- `experiment_synthetic2.m` (generates Figure 3 of the above paper)
	- `experiment_italy.m` (generates Figures 4-6 of the above paper)

  Please refer to the comments inside the scripts for more information.

- The `private` directory contains the MATLAB® implementation of the models described in the paper along with auxiliary functions.

- The `data` directory contains input data for the models and routines for data ingestion from the primary sources. Inside this repository, the `fetch_[...].py` files fetch the Google Community Reports, the baseline mobility in Italy as well as the epidemiological data from the relevant sources. Running the `extract_[...].m` prepares the MATLAB® `.mat` files that may be used for the simulations referring to the COVID-19 pandemic in Italy.

# Epidemiological indices with multiple circulating pathogen strains

This repository contains code and data described in

C. Trevisin, L. Mari, M. Gatto, V. Colizza, A. Rinaldo.
Epidemicity indices and reproduction numbers from infectious disease data in connected human populations

## Repository structure

- The top-level directory contains MATLAB® scripts for running all the analyses described in the above paper and generating the relevant figures:

	- `experiment_variants.m` (generates all figures of the above paper)

- The `data` directory, in addition to the files indicated above, also contains input data regarding the variants prevalence in Italy. Running the `prepare_variants.m` prepares the MATLAB® `.mat` files that may be used for the simulations referring to the COVID-19 pandemic in Italy.
