#=======================================================================
# Makefile example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================
CMP      = gcc_MPI#intel_MPI #intel_serial #gcc_serial
exeName  = cfd_LPT
CFD_DEFS_Add=
LPT_DEFS_Add=
CFDOrder = 2
Coupling = 2

# ============================  CFDLPT  ============================
# Choose pre-processing options
#   -DDOUBLE_PREC - use double-precision
#   -DSAVE_SINGLE - Save 3D data in single-precision
# CFDOrder: 1=SecondOrder, 2=FourthOrder
# Coupling: 1=One-way coupling; 2:Two-way coupling

CFDLPT_dir =  ./CFD_LPT/
CommDEFS   = -DCompiled_With_MPI -DVectorOperator -DSAVE_SINGLE

CFD_DEFS   = -DSaveNode
LPT_DEFS   =
CFD_DEFS  += $(CFD_DEFS_Add)
LPT_DEFS  += $(LPT_DEFS_Add)

CFD_inc    = -I./ThirdParty/ -I./Common/
CFD_lib    = -L./ThirdParty/fftw/ -lfftw3 -lm

LPT_inc    = -I./Common/

# CFDLPT source files
ifeq ($(Coupling),2) 
  CFD_DEFS  += -DCFDLPT_TwoWay
  LPT_DEFS  += -DCFDLPT_TwoWay
endif

ifeq ($(CFDOrder),1) 
  LPT_DEFS += -DCFDSecondOrder
  
  CFD_dir:=  ./CFD_2nd/
  Commdir:=  ./Common/
  SrcT   := f2_Parameters.f90 f2_Variables.f90 f2_MeshAndMetries.f90 f2_BC_and_Halo.f90 f2_Tools.f90  f2_TScheme.f90  \
            f2_FlowType_Channel.f90 f2_FlowType_Duct.f90 f2_FlowType_TGVortex.f90 f2_FlowType_HIT.f90 f2_FlowType_AddedNew.f90 \
            f2_FlowCase.f90 f2_Poisson.f90 f2_IOAndVisu.f90 f2_DumpPlane.f90 f2_Stat_User.f90 f2_CFDSystem.f90
  CFD_src:= $(addprefix $(CFD_dir),${SrcT})
else ifeq ($(CFDOrder),2)
  LPT_DEFS += -DCFDFourthOrder
  
  CFD_dir:=  ./CFD_4th/
  Commdir:=  ./Common/
  SrcT   := f4_Parameters.f90 f4_Variables.f90 f4_WritePlane.f90 f4_MeshAndMetries.f90 f4_BC_and_Halo.f90 f4_Tools.f90 f4_TScheme.f90 \
            f4_FlowCase.f90 f4_Poisson.f90 f4_IOAndVisu.f90 f4_CFDSystem.f90
  CFD_src:= $(addprefix $(CFD_dir),${SrcT})
endif
SrcT   := mc_Decomp2d.f90 mc_TypeDef.f90 mc_Timer.f90 mc_LogInfo.f90 mc_FileOperator.f90
Commsrc:= $(addprefix $(Commdir),${SrcT})
  
SrcT   := lp_Parameters.f90 lp_Decomp_2d.f90 lp_Property.f90 lp_Geometry.f90 lp_Variables.f90 \
          lp_ContactSearchPW.f90 lp_Integration.f90 lp_Comm.f90 lp_Statistics.f90 lp_IOAndVisu.f90 \
          lp_System.f90 lp_Fpforce.f90 main_CFDLPT.f90
LPT_src:= $(addprefix $(CFDLPT_dir),${SrcT})

#-----------------------------------------------------------------------
# Normally no need to change anything below
#-----------------------------------------------------------------------
include compile_flag_inc.make

all: $(exeName)
Commobj   = $(Commsrc:%.f90=%.o)
CFD_obj   = $(CFD_src:%.f90=%.o)
CFDLPT_obj= $(LPT_src:%.f90=%.o)
$(exeName):$(Commobj) $(CFD_obj) $(CFDLPT_obj)
	         $(FortC) $(CFLAG) -o $@ $(Commobj) $(CFD_obj) $(CFDLPT_obj) $(CFD_lib)
$(Commobj):$(Commdir)%.o :$(Commdir)%.f90
	         $(FortC) $(CFLAG) $(CommDEFS) -c $<
	         @ mv $(@F) ${Commdir}
$(CFD_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(CFD_inc) $(CFD_DEFS) -c $<
	         @ mv $(@F) ${CFD_dir}
$(CFDLPT_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(LPT_inc) $(LPT_DEFS) -c $<
	         @ mv $(@F) ${CFDLPT_dir}

.PHONY: clean
clean:
	rm -fr *.o *.mod $(exeName) $(CFD_dir)*.o $(Commdir)*.o $(CFDLPT_dir)*.o
