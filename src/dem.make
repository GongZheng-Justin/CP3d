#=======================================================================
# Makefile example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================
CMP     =  gcc_MPI#intel_MPI #intel_serial #gcc_serial
exeName =  dem
DEM_DEFS_Add =#-DDEMObliqueCollideDry -DTestDEMRestart

# ============================ DEM part ============================
# DEM options
DEM_dir  =  ./DEM/
Commdir  =  ./Common/
CommDEFS = -DCompiled_With_MPI -DVectorOperator

DEM_DEFS =
DEM_DEFS += $(DEM_DEFS_Add)

DEM_inc  = -I./Common/

# DEM source files
SrcT   := mc_TypeDef.f90 mc_Timer.f90 mc_LogInfo.f90 mc_FileOperator.f90
Commsrc:= $(addprefix $(Commdir), ${SrcT})

SrcT   := sp_Parameters.f90 sp_Decomp_2d.f90 sp_Property.f90 sp_Geometry.f90 sp_Variables.f90 \
          sp_CL_and_CF.f90 sp_ContactSearchPW.f90 sp_Integration.f90 sp_Comm.f90 sp_NBS_Munjiza.f90 \
          sp_Hrchl_Munjiza.f90 sp_ContactSearch.f90 sp_IOAndVisu.f90 sp_System.f90 main_DEM.f90
DEM_src:= $(addprefix $(DEM_dir), ${SrcT})

#-----------------------------------------------------------------------
# Normally no need to change anything below
#-----------------------------------------------------------------------
include compile_flag_inc.make

all: $(exeName)
Commobj  = $(Commsrc:%.f90=%.o)
DEM_obj  = $(DEM_src:%.f90=%.o)
$(exeName):$(Commobj) $(DEM_obj)
	         $(FortC) $(CFLAG) -o $@ $(Commobj) $(DEM_obj)
$(Commobj):$(Commdir)%.o :$(Commdir)%.f90
	         $(FortC) $(CFLAG) $(CommDEFS) -c $<
	         @ mv $(@F) ${Commdir}
$(DEM_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(DEM_inc) $(DEM_DEFS) -c $<
	         @ mv $(@F) ${DEM_dir}

.PHONY: clean
clean:
	rm -fr *.o *.mod $(exeName) $(DEM_dir)*.o $(Commdir)*.o
