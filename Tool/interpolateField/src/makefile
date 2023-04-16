#=======================================================================
# Makefile example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================
CMP          =  gcc_MPI#intel_MPI
Post_exe     =  interpolateField
IsSolveScalar=  1
POST_DEFS_Add=

# ============================ Post part ============================
# Post options

Post_lib  = -lm
post_Def :=

# Scalar or not
ifeq ($(IsSolveScalar),1)
  post_Def += -DScalarFlow
endif
post_Def += $(POST_DEFS_Add)

# Post source files
SrcF  = m_Decomp2d.f90 m_TypeDef.f90 m_LogInfo.f90 m_Parameters.f90 m_MeshAndMetries.f90 m_BC_and_Halo.f90 m_Interp.f90 main.f90

#-----------------------------------------------------------------------
# Normally no need to change anything below
#-----------------------------------------------------------------------
ifeq ($(CMP),intel_MPI)
  FortC = mpiifort
  CFLAG = -fpp -O3 -mcmodel=large -fp-model fast=2 -g# -heap-arrays 20 #-ipo -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2
else ifeq ($(CMP),intel_serial)
  FortC = ifort
  CFLAG = -fpp -O3 -mcmodel=large -fp-model fast=2 -g# -heap-arrays 20 #-ipo -xSSE4.2 -axAVX,CORE-AVX-I,CORE-AVX2
else ifeq ($(CMP),gcc_MPI)
  FortC = mpif90
  CFLAG = -cpp -O3 -Wall -fbacktrace -g -ffree-line-length-none -mcmodel=large -funroll-loops -floop-optimize
else ifeq ($(CMP),gcc_serial) 
  FortC = gfortran
  CFLAG = -cpp -O3 -Wall -fbacktrace -g -ffree-line-length-none -mcmodel=large -funroll-loops -floop-optimize
endif

all: $(Post_exe)
Post_obj   = $(SrcF:%.f90=%.o)
$(Post_exe): $(Post_objC) $(Post_obj)
	        $(FortC) $(CFLAG) -o $@  $(Post_obj) $(Post_lib)
$(Post_obj):%.o : %.f90
	        $(FortC) $(CFLAG) $(post_Def) -c $<

.PHONY: clean
clean:
	rm -fr  *.o *.mod $(Post_exe)  *o
