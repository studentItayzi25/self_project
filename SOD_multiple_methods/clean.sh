#!/usr/bin/env bash

# Files to keep
KEEP_FILES=("reimann.f90" "shock_tube_1D_multiple_methods.f90" "variables.f90" "calc_sub.f90" "plotpy.py" "clean.sh" "run.sh" "build.sh")

# Remove all files except the ones in KEEP_FILES
echo "Cleaning up files..."
for file in *; do
    if [[ ! " ${KEEP_FILES[@]} " =~ " ${file} " ]]; then
        rm -f "$file"
    fi
done

echo "Cleanup complete. Kept files: ${KEEP_FILES[@]}"