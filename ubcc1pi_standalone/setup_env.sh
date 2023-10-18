#!/bin/bash
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
# See https://cdcvs.fnal.gov/redmine/projects/uboonecode/wiki/Mrb_v5_breaking_change
setup uboonecode v08_00_00_57 -q e17:prof
unsetup larbatch
setup larbatch v01_51_15
source /uboone/app/users/jdetje/cc1pi_handover_57_2/localProducts_larsoft_v08_05_00_17_e17_prof/setup
mrbslp
ource /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
# See https://cdcvs.fnal.gov/redmine/projects/uboonecode/wiki/Mrb_v5_breaking_change
setup uboonecode v08_00_00_57 -q e17:prof
unsetup larbatch
setup larbatch v01_51_15
source /uboone/app/users/jdetje/cc1pi_handover_57_2/localProducts_larsoft_v08_05_00_17_e17_prof/setup
mrbslp
echo '-------------------Setup Done---------------------'cho '-------------------Setup Done---------------------'
