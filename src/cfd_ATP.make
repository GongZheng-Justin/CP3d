#=======================================================================
# Makefile example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================
CMP      = gcc_MPI#intel_MPI #intel_serial #gcc_serial
exeName  = cfd_ATP
CFD_DEFS_Add=
ATP_DEFS_Add=
CFDOrder = 2


# ============================  CFDATP  ============================
# Choose pre-processing options
#   -DDOUBLE_PREC - use double-precision
#   -DSAVE_SINGLE - Save 3D data in single-precision
# CFDOrder: 1=SecondOrder, 2=FourthOrder
# Coupling: 1=One-way coupling; 2:Two-way coupling

CFDATP_dir =  ./CFD_ATP/
CommDEFS   = -DCompiled_With_MPI -DVectorOperator -DSAVE_SINGLE

CFD_DEFS   = -DSaveNode
ATP_DEFS   =
CFD_DEFS  += $(CFD_DEFS_Add)
ATP_DEFS  += $(ATP_DEFS_Add)

CFD_inc    = -I./ThirdParty/ -I./Common/
CFD_lib    = -L./ThirdParty/fftw/ -lfftw3 -lm

ATP_inc    = -I./Common/

# CFDATP source files
ifeq ($(CFDOrder),1) 
  ATP_DEFS += -DCFDSecondOrder
  
  CFD_dir:=  ./CFD_2nd/
  Commdir:=  ./Common/
  SrcT   := f2_Parameters.f90 f2_Variables.f90 f2_MeshAndMetries.f90 f2_BC_and_Halo.f90 f2_Tools.f90  f2_TScheme.f90  \
            f2_FlowType_Channel.f90 f2_FlowType_Duct.f90 f2_FlowType_TGVortex.f90 f2_FlowType_HIT.f90 f2_FlowType_AddedNew.f90 \
            f2_FlowCase.f90 f2_Poisson.f90 f2_IOAndVisu.f90 f2_DumpPlane.f90 f2_Stat_User.f90 f2_CFDSystem.f90
  CFD_src:= $(addprefix $(CFD_dir),${SrcT})
else ifeq ($(CFDOrder),2)
  ATP_DEFS += -DCFDFourthOrder
  
  CFD_dir:=  ./CFD_4th/
  Commdir:=  ./Common/
  SrcT   := f4_Parameters.f90 f4_Variables.f90 f4_WritePlane.f90 f4_MeshAndMetries.f90 f4_BC_and_Halo.f90 f4_Tools.f90 f4_TScheme.f90 \
            f4_FlowCase.f90 f4_Poisson.f90 f4_IOAndVisu.f90 f4_CFDSystem.f90
  CFD_src:= $(addprefix $(CFD_dir),${SrcT})
endif
  
SrcT   := mc_Decomp2d.f90 mc_TypeDef.f90 mc_Timer.f90 mc_LogInfo.f90 mc_FileOperator.f90
Commsrc:= $(addprefix $(Commdir),${SrcT})
  
SrcT   := ap_Parameters.f90 ap_Decomp_2d.f90 ap_Property.f90 ap_Variables.f90 ap_ContactSearchPW.f90 \
          ap_Integration.f90 ap_Comm.f90 ap_Statistics.f90 ap_IOAndVisu.f90 ap_System.f90 \
          ap_Fpforce.f90 main_CFDATP.f90
ATP_src:= $(addprefix $(CFDATP_dir),${SrcT})

#-----------------------------------------------------------------------
# Normally no need to change anything below
#-----------------------------------------------------------------------
include compile_flag_inc.make

all: $(exeName)
Commobj   = $(Commsrc:%.f90=%.o)
CFD_obj   = $(CFD_src:%.f90=%.o)
CFDATP_obj= $(ATP_src:%.f90=%.o)
$(exeName):$(Commobj) $(CFD_obj) $(CFDATP_obj)
	         $(FortC) $(CFLAG) -o $@ $(Commobj) $(CFD_obj) $(CFDATP_obj) $(CFD_lib)
$(Commobj):$(Commdir)%.o :$(Commdir)%.f90
	         $(FortC) $(CFLAG) $(CommDEFS) -c $<
	         @ mv $(@F) ${Commdir}
$(CFD_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(CFD_inc) $(CFD_DEFS) -c $<
	         @ mv $(@F) ${CFD_dir}
$(CFDATP_obj):%.o :%.f90
	         $(FortC) $(CFLAG) $(ATP_inc) $(ATP_DEFS) -c $<
	         @ mv $(@F) ${CFDATP_dir}

.PHONY: clean
clean:
	rm -fr *.o *.mod $(exeName) $(CFD_dir)*.o $(Commdir)*.o $(CFDATP_dir)*.o
