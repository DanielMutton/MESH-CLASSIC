.PHONY: all clean
# Serial compiler
COMPILER = gfortran
# Serial Include Flags
IFLAGS = -O3 -g -I/usr/include
# Serial Library Flags
LFLAGS = -lnetcdff
# Additional flags
LFLAGS += -ldl -lz -lm -g -fdefault-real-8 -ffree-line-length-none -fbacktrace -ffpe-trap=invalid,zero,overflow  -fbounds-check

#-Wall -Wextra


# RECIPES
# Compile object files from .F90 sources
%.o: src/%.F90
	$(COMPILER) $(IFLAGS) -c $< -o $@ $(LFLAGS)

# Compile object files from .f90 sources
%.o: src/%.f90
	$(COMPILER) $(IFLAGS) -c $< -o $@ $(LFLAGS)

# Compile object files from .f (Fortran 77) sources
%.o: src/%.f
	$(COMPILER) $(IFLAGS) -c $< -o $@ $(LFLAGS)

all: fileIOModule.o initFileConverter.o
	$(COMPILER) $(IFLAGS) -o bin/initFileConverter fileIOModule.o initFileConverter.o $(LFLAGS)

# "make clean" removes all object files and executables
clean:
	rm -f *.o *.mod *~ core bin/initFileConverter