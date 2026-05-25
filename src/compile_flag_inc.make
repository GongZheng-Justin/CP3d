ifeq ($(CMP),intel_MPI)
  FortC = mpiifort
  CFLAG = -free -fpp -i4 -r8 -traceback -Ofast -mcmodel=large -g -heap-arrays 20 -fp-model fast#-CB -CU
else ifeq ($(CMP),intel_serial)
  FortC = ifort
  CFLAG = -free -fpp -i4 -r8 -traceback -Ofast -mcmodel=large -g -heap-arrays 20 -fp-model fast#-CB -CU
else ifeq ($(CMP),gcc_MPI)
  FortC = mpif90
  CFLAG = -ffree-line-length-none -cpp -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fbacktrace -Wall \
          -lm -g -O3 -mcmodel=large -funroll-loops -floop-optimize #-march=native #-fcheck=all
else ifeq ($(CMP),gcc_serial) 
  FortC = gfortran
  CFLAG = -ffree-line-length-none -cpp -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fbacktrace -Wall \
          -lm -g -O3 -mcmodel=large -funroll-loops -floop-optimize #-march=native #-fcheck=all
endif
