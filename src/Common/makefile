#=======================================================================
# Makefile example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================
CMP     =  gcc_MPI#intel_MPI
exeName =  MPI_Decomp2d
CFD_DEFS_Add =#-DFFTW_1D
IsSolveScalar= 0

# ============================ CFD part ============================
# Choose pre-processing options
#   -DOverWriteFFT
#   -DSaveNode
#   -DFFTW_1D

CFD_dir    =  ./
Commdir    =  ./
CommDEFS   = -DSAVE_SINGLE

CFD_inc =
CFD_lib =

# CFD source files
SrcT   := m_Decomp2d.f90
Commsrc:= $(addprefix $(Commdir), ${SrcT})

SrcT   := testDecomp2d.f90
CFD_src:= $(addprefix $(CFD_dir), ${SrcT})

#-----------------------------------------------------------------------
# Normally no need to change anything below
#-----------------------------------------------------------------------
ifeq ($(CMP),intel_MPI)
  FortC = mpiifort
  CFLAG = -fpp -O3 -mcmodel=large -fp-model fast=2#-ipo -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2
else ifeq ($(CMP),intel_serial)
  FortC = ifort
  CFLAG = -fpp -O3 -mcmodel=large -fp-model fast=2#-ipo -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2
else ifeq ($(CMP),gcc_MPI)
  FortC = mpif90
  CFLAG = -cpp -O3 -Wall -fbacktrace -g -ffree-line-length-none -mcmodel=large -funroll-loops -floop-optimize
else ifeq ($(CMP),gcc_serial) 
  FortC = gfortran
  CFLAG = -cpp -O3 -Wall -fbacktrace -g -ffree-line-length-none -mcmodel=large -funroll-loops -floop-optimize
endif

all: $(exeName)
Commobj  = $(Commsrc:%.f90=%.o)
CFD_obj  = $(CFD_src:%.f90=%.o)
$(exeName):$(Commobj) $(CFD_obj)
	         $(FortC) $(CFLAG) -o $@ $(Commobj) $(CFD_obj) $(CFD_lib)
$(Commobj):$(Commdir)%.o :$(Commdir)%.f90
	         $(FortC) $(CFLAG) $(CommDEFS) -c $<
$(CFD_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(CFD_inc) $(CFD_DEFS) -c $<

.PHONY: clean
clean:
	rm -fr *.o *.mod $(exeName) $(CFD_dir)*.o $(Commdir)*.o
