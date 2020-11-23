#!/bin/bash

# Recursively find all files with the supplied extensions
for EXT in {h,cxx,cc,fcl,xml}; do
    for FILE in $(find . -type f -name "*.${EXT}"); do

        # Remove any trailing whitespace
        sed --in-place 's/[[:space:]]\+$//' $FILE
    done
done
