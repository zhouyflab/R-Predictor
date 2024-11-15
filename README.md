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
Manually installing R-Predictor can be cumbersome, but fortunately, these tools it depends on are easy to work with.
- [Pfam36.0](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam36.0/)
~~~
conda create -n pfam_scan
source activate pfam_scan
conda install -c bioconda pfam_scan hmmer hmmer2 -y
~~~
- [SignalP 6.0](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md)
~~~
conda create -n signalp
source activate signalp
conda install -c predector signalp6
signalp6-register ./signalp-6.0h.fast.tar.gz
~~~
- [TMHMM-2.0](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=tmhmm&version=2.0c&packageversion=2.0c&platform=Linux)
~~~
The installation details are as follows:
1.Click the hyperlink to download the installation package.
2.Modify the first line of the `tmhmm` in the bin folder to use your own perl path, which can be found by running `which perl`. On line 33, change `$opt_basedir="./tmhmm-2.0c"`.
3.Modify the first line of the `tmhmmformat.pl` to use your own perl path.
~~~

