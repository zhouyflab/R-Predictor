# R-Predictor
## Table of Contents
- [Introduction](#Introduction)
- [Installation](#Installation)
  - [Quick installation using Docker](#QuickinstallationusingDocker)
  - [Manual installation](#Manualintallation)
- [Inputs](#Inputs)
- [Outputs](#Outputs)
- [R-Predictor usage](#R-Predictorusage)
- [Citations](#Citations)
- [Acknowledgements](#Acknowledgements)
## Introduction
This pipeline is designed to automatically annotate 15 distinct domain topologies of disease resistance genes across the entire plant genome.

R-Predictor, designed for the de novo annotation of various R genes integrate four modules for data pre-processing and the identification of different types of proteins. Each module incorporates customized filtering scripts and the best-performing methods identified through benchmarking.

![示例图片](images/pipeline.png)
## Installation
There are several ways to install R-Predictor. You just need to find the best that works for your system.

### Quick installation using Docker
The Docker version of the R-Predictor will be released shortly!!!

### Manual installation
Manually installing R-Predictor can be cumbersome, but fortunately, the tools it depends on are easy to work with.
- [PfamScan]
'conda create -n pfam_scan'
