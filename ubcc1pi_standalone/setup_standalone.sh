#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if [ ! -d "${SCRIPT_DIR}" ]; then
    echo "ERROR - Can't locate the ubcc1pi_standalone dir."
    return 1
fi

export UBCC1PI_STANDALONE_DIR=$SCRIPT_DIR
echo "Set UBCC1PI_STANDALONE_DIR=${UBCC1PI_STANDALONE_DIR}"

alias ubcc1pi="root -l -b $UBCC1PI_STANDALONE_DIR/compile.C"
echo "To run the standalone analysis code, start root with:"
echo " $ ubcc1pi"

return 0
