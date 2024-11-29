#!/usr/bin/env bash

# Define the executable name
EXE="shock_tube"

# Check if the executable exists
if [[ ! -f $EXE ]]; then
    echo "Executable $EXE not found. Please run build.sh first."
    exit 1
fi

# Run the executable
echo "Running the executable..."
./$EXE

# Run the Python plotting script
echo "Running the Python plotting script..."
python3 plotpy.py

echo "Execution complete."