#!/bin/bash

THIS_DIR=`basename $(pwd)`

if [ $THIS_DIR != "ubcc1pi_standalone" ]; then
    echo "ERROR - you must source this setup script while in the directory ubcc1pi_standalone"
    return 1
fi

export UBCC1PI_STANDALONE_DIR=`pwd`
alias ubcc1pi="root -l $UBCC1PI_STANDALONE_DIR/compile.C"

echo "To run the standalone analysis code, start root with:"
echo " $ ubcc1pi"

return 0
