#=======================================================================
# Makefile example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================
CMP            =  gcc_MPI#intel_MPI #intel_serial #gcc_serial
exeName        =  cfd_ACM
CFD_DEFS_Add   =
ACM_DEFS_Add   =#-DObliqueWallTest -DRotateOnly -DChanBraunJFM2011
CFDACM_DEFS_Add=#-DIBMDistributeLinear -DSeveralSphereInfo -DChanBraunJFM2011 -DTest_IBM_MPI

# ============================  CFDACM  ============================
# Choose pre-processing options
#   -DDOUBLE_PREC  - use double-precision
#   -DSAVE_SINGLE  - Save 3D data in single-precision
#   -DChanBraunJFM2011 - Test the work of Chan-Braun et al., J.Fluid Mech. (2011), vol.684,pp.441-474
# CFDACM options

CFD_dir     =  ./CFD_2nd/
Commdir     =  ./Common/
ACM_dir     =  ./DEM/
CFDACM_dir  =  ./CFD_ACM/
CommDEFS    = -DCompiled_With_MPI -DVectorOperator -DSAVE_SINGLE

CFD_DEFS    = -DCFDACM -DSaveNode 
ACM_DEFS    = -DCFDACM #-DOnlyDumpFpForce
CFDACM_DEFS = #-DSeveralSphereInfo#-DIBMDistributeLinear
CFD_DEFS   += $(CFD_DEFS_Add)
ACM_DEFS   += $(ACM_DEFS_Add)
CFDACM_DEFS+= $(CFDACM_DEFS_Add)

CFD_inc    = -I./ThirdParty/ -I./Common/
CFD_lib    = -L./ThirdParty/fftw/ -lfftw3 -lm

ACM_inc    = -I./Common/
CFDACM_inc = -I./Common/

# CFDACM source files
SrcT   := mc_Decomp2d.f90 mc_TypeDef.f90 mc_Timer.f90 mc_LogInfo.f90 mc_FileOperator.f90 mc_EqualSphere.f90
Commsrc:= $(addprefix $(Commdir), ${SrcT})

SrcT   := f2_Parameters.f90 f2_Variables.f90 f2_MeshAndMetries.f90 f2_BC_and_Halo.f90 f2_Tools.f90  f2_TScheme.f90  \
          f2_FlowType_Channel.f90 f2_FlowType_Duct.f90 f2_FlowType_TGVortex.f90 f2_FlowType_HIT.f90 f2_FlowType_AddedNew.f90 \
          f2_FlowCase.f90 f2_Poisson.f90 f2_IOAndVisu.f90 f2_DumpPlane.f90 f2_Stat_User.f90 f2_CFDSystem.f90         
CFD_src:= $(addprefix $(CFD_dir), ${SrcT})

SrcT   := sp_Parameters.f90 sp_Decomp_2d.f90 sp_Property.f90 sp_Geometry.f90 sp_Variables.f90     \
          sp_CL_and_CF.f90 sp_ContactSearchPW.f90 ac_Integration.f90 sp_Comm.f90 sp_IOAndVisu.f90 \
          sp_NBS_Munjiza.f90 sp_Hrchl_Munjiza.f90 sp_ContactSearch.f90 sp_DumpPrtcl.f90 ac_System.f90
ACM_src:= $(addprefix $(ACM_dir), ${SrcT})

SrcT   := ca_BC_and_Halo.f90 ca_IBM.f90 ca_IBM_implicit.f90 ca_Statistics.f90 ca_System.f90 main_CFDACM.f90
CA_src := $(addprefix $(CFDACM_dir), ${SrcT})

#-----------------------------------------------------------------------
# Normally no need to change anything below
#-----------------------------------------------------------------------
include compile_flag_inc.make

all: $(exeName)
Commobj   = $(Commsrc:%.f90=%.o)
CFD_obj   = $(CFD_src:%.f90=%.o)
ACM_obj   = $(ACM_src:%.f90=%.o)
CFDACM_obj= $(CA_src:%.f90=%.o)
$(exeName):$(Commobj) $(CFD_obj) $(ACM_obj) $(CFDACM_obj)
	         $(FortC) $(CFLAG) -o $@ $(Commobj) $(CFD_obj) $(ACM_obj) $(CFDACM_obj) $(CFD_lib)
$(Commobj):$(Commdir)%.o :$(Commdir)%.f90
	         $(FortC) $(CFLAG) $(CommDEFS) -c $<
	         @ mv $(@F) ${Commdir}
$(CFD_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(CFD_inc) $(CFD_DEFS) -c $<
	         @ mv $(@F) ${CFD_dir}
$(ACM_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(ACM_inc) $(ACM_DEFS) -c $<
	         @ mv $(@F) ${ACM_dir}
$(CFDACM_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(CFDACM_inc) $(CFDACM_DEFS) -c $<
	         @ mv $(@F) ${CFDACM_dir}

.PHONY: clean
clean:
	rm -fr *.o *.mod $(exeName) $(CFD_dir)*.o $(Commdir)*.o $(ACM_dir)*.o $(CFDACM_dir)*.o
