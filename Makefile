# Compiler
FC = gfortran

# Compiler flags
FFLAGS = -O2 -Wall

# Source files
SRC = $(wildcard *.f90)

# Object files
OBJ = $(SRC:.f90=.o)

# Executable
EXEC = shock_tube

# Default target
all: $(EXEC)

# Link object files to create executable
$(EXEC): $(OBJ)
    $(FC) $(FFLAGS) -o $@ $^

# Compile source files to object files
%.o: %.f90
    $(FC) $(FFLAGS) -c $< -o $@

# Run the executable
run: $(EXEC)
    ./$(EXEC)

# Clean up object files and executable
clean:
    rm -f $(OBJ) $(EXEC)

# Clean up all files except .f90 files
clean_all:
    find . -type f ! -name '*.f90' -delete

.PHONY: all run clean clean_all