#!/usr/bin/env bash

# Define the Fortran compiler
FC=gfortran

# Define the compiler flags
FFLAGS="-O2 -Wall"

# Define the source files
SRC="variables.f90 grid.f90 boundary_conditions.f90 interpolations.f90 prim_con.f90 fluxes.f90 residuals.f90 time_int.f90 shock_tube_1D.f90"

# Define the executable name
EXE="shock_tube"

# Compile the source files into object files
echo "Compiling source files..."
for file in $SRC; do
    $FC $FFLAGS -c $file
done

# Link the object files to create the executable
echo "Linking object files..."
$FC $FFLAGS -o $EXE *.o

# Clean up object files
echo "Cleaning up..."
rm -f *.o

echo "Build complete. Executable created: $EXE"