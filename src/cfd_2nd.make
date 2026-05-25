#=======================================================================
# Makefile example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================
CMP     =  gcc_MPI#intel_MPI #intel_serial #gcc_serial
exeName =  cfd_2nd
CFD_DEFS_Add =

# ============================ CFD part ============================
# Choose pre-processing options
#   -DSaveNode
# CFD options

CFD_dir    =  ./CFD_2nd/
Commdir    =  ./Common/
CommDEFS   = -DCompiled_With_MPI -DSAVE_SINGLE

CFD_DEFS   = -DSaveNode
CFD_DEFS  += $(CFD_DEFS_Add)

CFD_inc    = -I./ThirdParty/ -I./Common/
CFD_lib    = -L./ThirdParty/fftw/ -lfftw3 -lm
    
# CFD source files
SrcT   := mc_Decomp2d.f90 mc_TypeDef.f90 mc_Timer.f90 mc_LogInfo.f90 mc_FileOperator.f90
Commsrc:= $(addprefix $(Commdir), ${SrcT})

SrcT   := f2_Parameters.f90 f2_Variables.f90 f2_MeshAndMetries.f90 f2_BC_and_Halo.f90 f2_Tools.f90  f2_TScheme.f90  \
          f2_FlowType_Channel.f90 f2_FlowType_Duct.f90 f2_FlowType_TGVortex.f90 f2_FlowType_HIT.f90 f2_FlowType_AddedNew.f90 \
          f2_FlowCase.f90 f2_Poisson.f90 f2_IOAndVisu.f90 f2_DumpPlane.f90 f2_Stat_User.f90 f2_CFDSystem.f90 main_CFD2nd.f90
CFD_src:= $(addprefix $(CFD_dir), ${SrcT})

#-----------------------------------------------------------------------
# Normally no need to change anything below
#-----------------------------------------------------------------------
include compile_flag_inc.make

all: $(exeName)
Commobj  = $(Commsrc:%.f90=%.o)
CFD_obj  = $(CFD_src:%.f90=%.o)
$(exeName):$(Commobj) $(CFD_obj)
	         $(FortC) $(CFLAG) -o $@ $(Commobj) $(CFD_obj) $(CFD_lib)
$(Commobj):$(Commdir)%.o :$(Commdir)%.f90
	         $(FortC) $(CFLAG) $(CommDEFS) -c $<
	         @ mv $(@F) ${Commdir}
$(CFD_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(CFD_inc) $(CFD_DEFS) -c $<
	         @ mv $(@F) ${CFD_dir}

.PHONY: clean
clean:
	rm -fr *.o *.mod $(exeName) $(CFD_dir)*.o $(Commdir)*.o
