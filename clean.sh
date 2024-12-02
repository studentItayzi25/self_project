#!/usr/bin/env bash

# Files to keep
KEEP_FILES=("variables.f90" "grid.f90" "boundary_conditions.f90" "interpolations.f90" "prim_con.f90" "fluxes.f90" "residuals.f90" "time_int.f90" "shock_tube_1D.f90" "plotpy.py" "clean.sh" "run.sh" "build.sh" "Makefile")

# Remove all files except the ones in KEEP_FILES
echo "Cleaning up files..."
for file in *; do
    if [[ ! " ${KEEP_FILES[@]} " =~ " ${file} " ]]; then
        rm -f "$file"
    fi
done

echo "Cleanup complete. Kept files: ${KEEP_FILES[@]}"