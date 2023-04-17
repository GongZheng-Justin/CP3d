#ifdef IBMDistribute2
  idxp_interp=ic
  idyp_interp=jc
  idzp_interp=kc
  if(LagrangeP%x>xc(ic)) then
    idxc_interp= ic
  else
    idxc_interp= ic-1
  endif
  if(LagrangeP%y>yc(jc)) then
    idyc_interp= jc
  else
    idyc_interp= jc-1
  endif
  if(LagrangeP%z>zc(kc)) then
    idzc_interp= kc
  else
    idzc_interp= kc-1
  endif
#else
  idxc_interp=ic-1
  idyc_interp=jc-1
  idzc_interp=kc-1
  if(LagrangeP%x>xc(ic)) then
    idxp_interp= ic
  else
    idxp_interp= ic-1
  endif
  if(LagrangeP%y>yc(jc)) then
    idyp_interp= jc
  else
    idyp_interp= jc-1
  endif
  if(LagrangeP%z>zc(kc)) then
    idzp_interp= kc
  else
    idzp_interp= kc-1
  endif
#endif
#ifdef clc_Point_indxyz_IBPFix
  IBPFix_indxyz(1,nIBPFix)=idxc_interp
  IBPFix_indxyz(2,nIBPFix)=idxp_interp
  IBPFix_indxyz(3,nIBPFix)=idyc_interp
  IBPFix_indxyz(4,nIBPFix)=idyp_interp
  IBPFix_indxyz(5,nIBPFix)=idzc_interp
  IBPFix_indxyz(6,nIBPFix)=idzp_interp
#endif
#ifdef clc_Point_indxyz_IBPMove
  IBP_indxyz(1,nIBP)=idxc_interp
  IBP_indxyz(2,nIBP)=idxp_interp
  IBP_indxyz(3,nIBP)=idyc_interp
  IBP_indxyz(4,nIBP)=idyp_interp
  IBP_indxyz(5,nIBP)=idzc_interp
  IBP_indxyz(6,nIBP)=idzp_interp
#endif
