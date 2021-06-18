# ubcc1pi
A package for the selection and analysis of charged current single pion muon neutrino interactions in MicroBooNE

## Overview
This repository contains code which can be run in LArSoft or as a standalone package. The LArSoft component is used to exctract relevant
information from art ROOT files. This information can either be persisted to a lightweight analysis ROOT file for standalone analysis, or be
analysed within the LArSoft framework. The actual physics logic of the analysis is factored out and doesn't depend on LArSoft at all.

## Installation as a LArSoft package
Setup your uBooNE suite LArSoft development area using `mrb` as usual (see
[here](https://cdcvs.fnal.gov/redmine/projects/uboonecode/wiki/Uboone_guide) for more details). Then follow the steps below:

```
# Clone this package in your sources directory
# Note that we here use the recurse-submodules option because we also want to pick up the ubsmear submodule
cd $MRB_SOURCE
git clone --recurse-submodules https://github.com/a-d-smith/ubcc1pi.git

# Update your CMakeLists.txt so that mrb knows to build this package
mrb updateDepsCM

# Set up your development environment and check for dependency issues
mrbsetenv

# Install the package
mrb install -j4

# Setup your local products
mrbslp
```

## Installation for standalone analysis

```
# Clone this repository in your desired area
git clone --recurse-submodules https://github.com/a-d-smith/ubcc1pi.git

# See instructions in ubcc1pi_standalone for compilation
```

## Documentation

Please see <https://a-d-smith.github.io/ubcc1pi/> for documentation. Or look in [/docs/](docs/) to build your own.
