# Automatically detect environment and available libraries
OBJECTS=evalsurf.o
LIBOBJS=hddata.o combinatorial.o progdata.o libutil.o libsym.o libinternal.o \
	io.o potlib.o getver.o

# Allow environment to override compiler
FC ?= gfortran
FCOPT ?= -g -fdefault-integer-8 -fopenmp

# Detect CI environment
CI ?= $(shell printenv CI || echo "")

# Auto-detect library paths if not specified
LIBRARY_DIRS ?= $(shell if [ -d "/usr/lib/x86_64-linux-gnu" ]; then echo "/usr/lib/x86_64-linux-gnu"; fi)
INCLUDE_DIRS ?= $(shell if [ -d "/usr/include" ]; then echo "/usr/include"; fi)

# Additional paths to try for BLAS/LAPACK
POSSIBLE_LAPACK_LIBS := -llapack -lopenblas -lflexiblas
POSSIBLE_BLAS_LIBS := -lblas -lopenblas -lflexiblas

# Choose the first available LAPACK/BLAS library
LAPACK_LIB ?= $(shell for lib in $(POSSIBLE_LAPACK_LIBS); do \
                 if ldconfig -p | grep -q $$lib; then \
                   echo $$lib; \
                   break; \
                 fi; \
               done)

BLAS_LIB ?= $(shell for lib in $(POSSIBLE_BLAS_LIBS); do \
              if ldconfig -p | grep -q $$lib; then \
                echo $$lib; \
                break; \
              fi; \
            done)

# Fallback to standard names if detection fails
LAPACK_LIB ?= -llapack
BLAS_LIB ?= -lblas

# Add library directories if specified
LDFLAGS = $(if $(LIBRARY_DIRS),-L$(LIBRARY_DIRS),)

all:	lib $(OBJECTS)
	@echo "Building evalsurf.x with FC=$(FC), LAPACK_LIB=$(LAPACK_LIB), BLAS_LIB=$(BLAS_LIB)"
	$(FC) $(FCOPT) -c evalsurf.f90
	$(FC) $(FCOPT) -o evalsurf.x evalsurf.o libsurfgen.a $(LDFLAGS) $(LAPACK_LIB) $(BLAS_LIB)

lib:	$(LIBOBJS)
	$(AR) -r -v libsurfgen.a $(LIBOBJS)
	@echo " Cleaning up..."
	rm -f $(LIBOBJS)
	rm -f ./*.mod

./%.o : ./lib/%.f90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

./%.o : ./lib/%.F90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

./%.o : ./%.f90
	@echo ' Building file: $<'
	$(FC) $(FCOPT) -c -o $@ $<
	@echo ' Finished building $<'
	@echo ''

clean :
	rm -f ./*.mod ./*.o