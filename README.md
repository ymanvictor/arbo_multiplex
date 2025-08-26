# A repository to analyse multiplex arbovirus serological data

## Data

Individual-level data are subject to participant confidentiality requirements under GDPR and national regulations in France, Peru, and Senegal, as well as the terms of local ethics approvals. As a result, raw individual-level datasets cannot be made publicly available.
De-identified summary data sufficient to reproduce the main analyses is included within this repository.

## Code

This code was written in R (version 4.4.3) running on MacOS 15.0.1. Bayesian finite mixture models were executed from R on the Stan platform (Stan Development Team. 2024, Stan Modeling Language version 2.35, https://mc-stan.org) using the cmdstanr package (cmdstanr: R Interface to 'CmdStan', version 0.7.1) as the R interface to Stan.

Before running this repository, ensure Stan’s command-line toolchain (CmdStan) is installed and available through the R interface CmdStanR, and that you have a working C/C++ compiler. Use CmdStanR to verify your toolchain and perform the CmdStan installation—the first build compiles locally and can take several minutes. Confirm that CmdStanR detects a valid CmdStan version and path before executing any analysis scripts.

* File naming implies order
* Code to run "cross-reactivity models" is available at https://github.com/nathoze/Mayaro
