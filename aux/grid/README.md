# Running on the grid

This directory contains the XML settings file for running this analysis code on the FERMI grid system using `project.py`.

## Getting set up
```
# Make sure you have your local products set up
mrbslp

# Update the tarball of the local products
./make_tar

# Update GridSettings.xml to run on the SAM input definition of your choosing - be sure to update the number of jobs too!
vim GridSettings.xml

# Choose the study you wish to run (or choose "all" to run all studies together). You can see the available stages by running:
project.py --xml GridSettings.xml --status
```

## Submitting and checking on jobs
```
# Run your stage (replacing <my stage> as required)
project.py --xml GridSettings.xml --stage <my stage> --submit

# Check on the status of your jobs
project.py --xml GridSettings.xml --stage <my stage> --status

# Once your jobs start to finish you can find the output directory by running
project.py --xml GridSettings.xml --stage <my stage> --outdir
```
