#=======================================================================
# Makefile example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================
CMP       =  gcc_MPI#intel_MPI #intel_serial #gcc_serial
exeName   =  cfd_DEM
CFD_DEFS_Add   =
DEM_DEFS_Add   =
CFDDEM_DEFS_Add=#-DSeveralSphereInfo

# ============================  CFDDEM  ============================
# Choose pre-processing options
#   -DDOUBLE_PREC - use double-precision
#   -DSAVE_SINGLE - Save 3D data in single-precision
# CFDDEM options

CFD_dir    =  ./CFD_2nd/
Commdir    =  ./Common/
DEM_dir    =  ./DEM/
CFDDEM_dir =  ./CFD_DEM/
CommDEFS   = -DCompiled_With_MPI -DVectorOperator -DSAVE_SINGLE

CFD_DEFS   = -DSaveNode -DCFDDEM
DEM_DEFS   = -DCFDDEM
CFDDEM_DEFS=
CFD_DEFS   += $(CFD_DEFS_Add)
DEM_DEFS   += $(DEM_DEFS_Add)
CFDDEM_DEFS+= $(CFDDEM_DEFS_Add)

CFD_inc    = -I./ThirdParty/ -I./Common/
CFD_lib    = -L./ThirdParty/fftw/ -lfftw3 -lm

DEM_inc    = -I./Common/
CFDDEM_inc = -I./Common/

# CFDDEM source files
SrcT   := mc_Decomp2d.f90 mc_TypeDef.f90 mc_Timer.f90 mc_LogInfo.f90 mc_FileOperator.f90
Commsrc:= $(addprefix $(Commdir), ${SrcT})

SrcT   := f2_Parameters.f90 f2_Variables.f90 f2_MeshAndMetries.f90 f2_BC_and_Halo.f90 f2_Tools.f90  f2_TScheme.f90  \
          f2_FlowType_Channel.f90 f2_FlowType_Duct.f90 f2_FlowType_TGVortex.f90 f2_FlowType_HIT.f90 f2_FlowType_AddedNew.f90 \
          f2_FlowCase.f90 f2_Poisson.f90 f2_IOAndVisu.f90 f2_DumpPlane.f90 f2_Stat_User.f90 f2_CFDSystem.f90
CFD_src:= $(addprefix $(CFD_dir), ${SrcT})

SrcT   := sp_Parameters.f90 sp_Decomp_2d.f90 sp_Property.f90 sp_Geometry.f90 sp_Variables.f90 \
          sp_CL_and_CF.f90 sp_ContactSearchPW.f90 sp_Integration.f90 sp_Comm.f90 sp_IOAndVisu.f90 \
          sp_NBS_Munjiza.f90 sp_Hrchl_Munjiza.f90 sp_ContactSearch.f90 sp_DumpPrtcl.f90 sp_System.f90
DEM_src:= $(addprefix $(DEM_dir), ${SrcT})

SrcT   := cd_FpForce.f90 cd_Statistics.f90 cd_System.f90 main_CFDDEM.f90
CD_src := $(addprefix $(CFDDEM_dir), ${SrcT})

#-----------------------------------------------------------------------
# Normally no need to change anything below
#-----------------------------------------------------------------------
include compile_flag_inc.make

all: $(exeName)
Commobj   = $(Commsrc:%.f90=%.o)
CFD_obj   = $(CFD_src:%.f90=%.o)
DEM_obj   = $(DEM_src:%.f90=%.o)
CFDDEM_obj= $(CD_src:%.f90=%.o)
$(exeName):$(Commobj) $(CFD_obj) $(DEM_obj) $(CFDDEM_obj)
	         $(FortC) $(CFLAG) -o $@ $(Commobj) $(CFD_obj) $(DEM_obj) $(CFDDEM_obj) $(CFD_lib)
$(Commobj):$(Commdir)%.o :$(Commdir)%.f90
	         $(FortC) $(CFLAG) $(CommDEFS) -c $<
	         @ mv $(@F) ${Commdir}
$(CFD_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(CFD_inc) $(CFD_DEFS) -c $<
	         @ mv $(@F) ${CFD_dir}
$(DEM_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(DEM_inc) $(DEM_DEFS) -c $<
	         @ mv $(@F) ${DEM_dir}
$(CFDDEM_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(CFDDEM_inc) $(CFDDEM_DEFS) -c $<
	         @ mv $(@F) ${CFDDEM_dir}

.PHONY: clean
clean:
	rm -fr *.o *.mod $(exeName) $(CFD_dir)*.o $(Commdir)*.o $(DEM_dir)*.o $(CFDDEM_dir)*.o
