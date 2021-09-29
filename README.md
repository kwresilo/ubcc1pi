# ubcc1pi
A package for the selection and analysis of charged current single pion muon neutrino interactions in MicroBooNE

## Overview
This repository contains code which can be run in LArSoft or as a standalone package. The LArSoft component is used to exctract relevant
information from art ROOT files. This information can either be persisted to a lightweight analysis ROOT file for standalone analysis, or be
analysed within the LArSoft framework. The actual physics logic of the analysis is factored out and doesn't depend on LArSoft at all.

## Installation as a LArSoft package
Setup your uBooNE suite LArSoft development area using `mrb` as usual (see
[here](https://cdcvs.fnal.gov/redmine/projects/uboonecode/wiki/Uboone_guide) for more details). Then follow the steps below:

```shell
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

### Additional steps
A small number of overlay files produce an error. To avoid this, it is necessary to also add `ubana` to your local products
and modify a few lines to catch an exception.

```shell
# Get a copy of ubana
cd $MRB_SOURCE
mrb g ubana
cd ubana
git checkout $UBANA_VERSION
```

In the file `/ubana/Calibrations/CalibrationdEdX_module.cc` add a try-catch block to the method `void ub::CalibrationdEdX::produce(art::Event & evt)`: 
```cpp
for (size_t j = 0; j<vdQdx.size(); ++j){
  float yzcorrection{1.0};
  float xcorrection{1.0};
  try{yzcorrection = energyCalibProvider.YZdqdxCorrection(planeID.Plane, vXYZ[j].Y(), vXYZ[j].Z());}
  catch(const lariov::IOVDataError &){std::cout << "Error - Unable to obtain yzcorrection in ubana/Calibrations/CalibrationdEdX_module.cc\n";}
  try{ xcorrection  = energyCalibProvider.XdqdxCorrection(planeID.Plane, vXYZ[j].X());}
  catch(const lariov::IOVDataError &){std::cout << "Error - Unable to obtain xcorrection in ubana/Calibrations/CalibrationdEdX_module.cc\n";}
  float elifetime  = elifetimeCalibProvider.Lifetime(); // [ms]
  float driftvelocity = detprop->DriftVelocity(); // [cm/us]
```

Repeat the initial installation steps starting with `mrbsetenv`. It might be necessary to do `unsetup cetbuildtools` and `unsetup cmake` first.
If there are incompatibilities with cetbuildtools between `ubcc1pi` and `ubana`, change `ubana/ups/product_deps` from this:
```
cetbuildtools           v7_xx_xx 	-       only_for_build
```
to this:
```
cetbuildtools           v7_04_00 	-       only_for_build
```

## Installation for standalone analysis

```shell
# Clone this repository in your desired area
git clone --recurse-submodules https://github.com/a-d-smith/ubcc1pi.git

# See instructions in ubcc1pi_standalone for compilation
```

## Documentation

Please see <https://a-d-smith.github.io/ubcc1pi/> for documentation. Or look in [/docs/](docs/) to build your own.
