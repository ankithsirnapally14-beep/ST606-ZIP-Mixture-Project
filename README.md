
# ST606-ZIP-Mixture-Project

## Zero-Inflated Poisson Mixture Models for Clustering Count Data

**Author:** Ankith Sirnapally
**Student ID:** 25253247
**Supervisor:** Dr. Keefe Murphy
**Module:** ST606 MSc Thesis Project
**Institution:** Maynooth University

## Project Overview

This project investigates clustering count data using Poisson Mixture Models and Zero-Inflated Poisson (ZIP) Mixture Models.

The aim is to determine whether explicitly modelling excess zeros through a ZIP mixture model improves clustering performance compared to the traditional Anscombe Transformation followed by Gaussian mixture clustering using mclust.

## Methods

* Poisson Mixture EM Algorithm
* Zero-Inflated Poisson (ZIP) Mixture EM Algorithm
* Anscombe Transformation + mclust
* Simulation Study
* Adjusted Rand Index (ARI) Evaluation

## Key Contributions

* Implemented Poisson Mixture EM algorithm.
* Extended the model to a ZIP Mixture EM algorithm.
* Added estimation of the latent variable U-hat.
* Updated lambda and omega estimation using U-hat.
* Conducted a simulation study across multiple settings.

## Main Results

| Method            | ARI   |
| ----------------- | ----- |
| Anscombe + mclust | 0.069 |
| Poisson EM        | 0.858 |
| ZIP EM            | 0.896 |

The ZIP EM model achieved the highest clustering accuracy and consistently outperformed Anscombe + mclust across all simulation settings.

## Repository Contents

* `final_code.R` – Main R implementation
* `ZIP_Mixture_Thesis_Presentation.pptx` – Thesis presentation slides

## Future Work

* Larger simulation studies
* Real-world count datasets
* Poisson Log-Normal mixture models
* Negative Binomial mixture models

## Acknowledgements

This work was completed as part of the ST606 MSc Thesis Project at Maynooth University.
