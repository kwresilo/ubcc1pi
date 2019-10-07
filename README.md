# ubcc1pi
A LArSoft package for the selection and analysis of charged current single pion muon neutrino interactions in MicroBooNE

## Installation
Setup your uBooNE suite LArSoft development area using `mrb` as usual (see [here](https://cdcvs.fnal.gov/redmine/projects/uboonecode/wiki/Uboone_guide) for more details). Then follow the steps below:

```
# Clone this package in your sources directory
cd $MRB_SOURCE
git clone https://github.com/a-d-smith/ubcc1pi.git

# Update your CMakeLists.txt so that mrb knows to build this package
mrb updateDepsCM

# Set up your development environment and check for dependency issues
mrbsetenv

# Install the package
mrb install -j4

# Setup your local products
mrbslp
```

## Documentation

Please see <https://a-d-smith.github.io/ubcc1pi/> for documentation (may be out of date). Or look in [/docs/](docs/) to build your own.
