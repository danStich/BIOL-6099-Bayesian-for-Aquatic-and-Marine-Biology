# BIOL-6099-Bayesian-for-Aquatic-and-Marine-Biology
## Overview
This repository holds data, code, and models used for BIOL-6099 (Bayesian Hierarchical Modeling for Aquatic and Marine Biology). Data may be from published or ongoing studies, so please contact us directly if you are interested in using them. This site is currently under development and may change without notice.

### Directory structure
The main directory contains R scripts from class examples using JAGS or Stan. These are named using the module number, a data set name, and the sampler used for analysis in each script, which are separated for modularization. The main directory contains several sub-directories used to organize data and models in addition to some organizational files such as a `.gitignore`, an R project file to make life easy, and this `README.md` file.

### Subdirectories
`data` Data files used for R scripts in the main directory and in `toy_examples`

`models` JAGS and Stan model files used in the example scripts

`toy_examples` Scripts for toy examples to demonstrate variety of applications and analysis tools ([brms](https://paulbuerkner.com/brms/), [rstanarm](https://mc-stan.org/rstanarm/articles/rstanarm.html), etc.) that are not covered in detail in this class (but they are awesome tools and you should use them when you can). Most examples in this directory are derived from silly disagreements or questions that required "rigorous" sampling and analysis with Bayesian statistics of course.

## Modules
R scripts in the main directory for class activities.

`01_cray_*.R` Bayesian hierarchical linear model. Uses log-10 transformed length and mass of rusty crayfish *Phoxonius rusticus* collected from 7 sites in the upper Susquehanna River to fit length-weight regressions while accounting for variation within and among sites. Used to demonstrate the 'bayesics' of writing and wrangling Bayesian hierarchical models "by hand" in JAGS and Stan.

`02_otsego_sav.R` Bayesian logistic regression for simulated aquatic plant survey at Otsego Lake, NY using [rstanarm](https://mc-stan.org/rstanarm/articles/rstanarm.html) and [rstan](https://mc-stan.org/rstan/) for comparing syntax and functionality.
