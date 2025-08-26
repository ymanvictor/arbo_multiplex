# A repository to analyse multiplex arbovirus serological data

## Data

Individual-level data are protected by participant confidentiality requirements under the GDPR and national regulations in France, Peru, and Senegal, as well as local ethics approvals. Consequently, raw individual-level datasets cannot be shared publicly. De-identified summary data sufficient to reproduce the main analyses are included in this repository.

## Code

All analyses were performed in R (v4.4.3) on macOS (v15.0.1). Bayesian finite mixture models were executed from R on the Stan platform (Stan Development Team. 2024, Stan Modeling Language version 2.35, https://mc-stan.org) using the cmdstanr package (cmdstanr: R Interface to 'CmdStan', version 0.7.1) as the R interface to Stan.

### Setup requirements

Before running any analysis scripts, ensure that Stan’s command-line toolchain (CmdStan) is installed and available through CmdStanR, and that a working C/C++ toolchain is present. Use CmdStanR to verify the toolchain and complete the CmdStan installation; the initial build compiles locally and may take several minutes. Confirm that CmdStanR detects a valid CmdStan version and path prior to execution.

### Execution order

Filenames are prefixed to indicate the intended order of execution

### Stan code

Stan code to run Bayesian finite mixture models is included within the /src directory

### Related code

Code to run the "cross-reactivity models" is available at https://github.com/nathoze/Mayaro
