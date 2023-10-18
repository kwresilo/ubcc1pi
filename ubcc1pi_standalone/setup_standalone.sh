#!/bin/bash

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
# See https://cdcvs.fnal.gov/redmine/projects/uboonecode/wiki/Mrb_v5_breaking_change
setup uboonecode v08_00_00_57 -q e17:prof
#unsetup larbatch
#setup larbatch v01_51_15
#source /uboone/app/users/jdetje/cc1pi_handover_57_2/localProducts_larsoft_v08_05_00_17_e17_prof/setup
#mrbslp
echo '-------------------Setup Done---------------------'

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if [ ! -d "${SCRIPT_DIR}" ]; then
    echo "ERROR - Can't locate the ubcc1pi_standalone dir."
    return 1
fi

export UBCC1PI_STANDALONE_DIR=$SCRIPT_DIR
echo "Set UBCC1PI_STANDALONE_DIR=${UBCC1PI_STANDALONE_DIR}"

alias ubcc1pi="root -l -b $UBCC1PI_STANDALONE_DIR/compile.cxx"
echo "To run the standalone analysis code, start root with:"
echo " $ ubcc1pi"

return 0
