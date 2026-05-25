#=======================================================================
# Makefile example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================
CMP            =  gcc_MPI#intel_MPI #intel_serial #gcc_serial
exeName        =  cfd_Porous_IBM
CFD_DEFS_Add   =
IBM_DEFS_Add=#-DIBMDistributeLinear -DSeveralSphereInfo -DChanBraunJFM2011 -DTest_IBM_MPI

# ============================  CFDIBM  ============================
# Choose pre-processing options
#   -DDOUBLE_PREC  - use double-precision
#   -DSAVE_SINGLE  - Save 3D data in single-precision
#   -DChanBraunJFM2011 - Test the work of Chan-Braun et al., J.Fluid Mech. (2011), vol.684,pp.441-474
# CFDIBM options

CFD_dir     =  ./CFD_2nd/
Commdir     =  ./Common/
CFDIBM_dir  =  ./CFD_Porous_IBM/
CommDEFS    = -DCompiled_With_MPI -DVectorOperator -DSAVE_SINGLE

CFD_DEFS    = -DCFDACM#-DSaveNode
CFDIBM_DEFS = #-DSeveralSphereInfo#-DIBMDistributeLinear
CFD_DEFS   += $(CFD_DEFS_Add)
CFDIBM_DEFS+= $(IBM_DEFS_Add)

CFD_inc    = -I./ThirdParty/ -I./Common/
CFD_lib    = -L./ThirdParty/fftw/ -lfftw3 -lm

CFDIBM_inc = -I./Common/

# CFDIBM source files
SrcT   := mc_Decomp2d.f90 mc_TypeDef.f90 mc_Timer.f90 mc_LogInfo.f90 mc_FileOperator.f90 mc_EqualSphere.f90
Commsrc:= $(addprefix $(Commdir), ${SrcT})

SrcT   := f2_Parameters.f90 f2_Variables.f90 f2_MeshAndMetries.f90 f2_BC_and_Halo.f90 f2_Tools.f90  f2_TScheme.f90  \
          f2_FlowType_Channel.f90 f2_FlowType_Duct.f90 f2_FlowType_TGVortex.f90 f2_FlowType_HIT.f90 f2_FlowType_AddedNew.f90 \
          f2_FlowCase.f90 f2_Poisson.f90 f2_IOAndVisu.f90 f2_DumpPlane.f90 f2_Stat_User.f90 f2_CFDSystem.f90         
CFD_src:= $(addprefix $(CFD_dir), ${SrcT})

SrcT   := pIBM_BC_and_Halo.f90 pIBM_IBM.f90 pIBM_System.f90 main_Porous_IBM.f90
CA_src := $(addprefix $(CFDIBM_dir), ${SrcT})

#-----------------------------------------------------------------------
# Normally no need to change anything below
#-----------------------------------------------------------------------
include compile_flag_inc.make

all: $(exeName)
Commobj   = $(Commsrc:%.f90=%.o)
CFD_obj   = $(CFD_src:%.f90=%.o)
CFDIBM_obj= $(CA_src:%.f90=%.o)
$(exeName):$(Commobj) $(CFD_obj) $(CFDIBM_obj)
	         $(FortC) $(CFLAG) -o $@ $(Commobj) $(CFD_obj) $(CFDIBM_obj) $(CFD_lib)
$(Commobj):$(Commdir)%.o :$(Commdir)%.f90
	         $(FortC) $(CFLAG) $(CommDEFS) -c $<
	         @ mv $(@F) ${Commdir}
$(CFD_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(CFD_inc) $(CFD_DEFS) -c $<
	         @ mv $(@F) ${CFD_dir}
$(CFDIBM_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(CFDIBM_inc) $(CFDIBM_DEFS) -c $<
	         @ mv $(@F) ${CFDIBM_dir}

.PHONY: clean
clean:
	rm -fr *.o *.mod $(exeName) $(CFD_dir)*.o $(Commdir)*.o $(CFDIBM_dir)*.o
