    Vel_w=((ri*Rvei+rj*Rvej).cross.Norm_v)
    Vrij = Veli-Velj + Vel_w ! Relative velocity
    vrn  = Vrij .dot. Norm_v
    if(VelRel_Init(ind)<-1.0_RK)VelRel_Init(ind)=abs(vrn)

    Vij_n= vrn*Norm_v        ! Normal  relative velocity
    Vij_t= Vrij-Vij_n        ! Tangent relative velocity
    
    ! Tangential overlap Vector
    Ovlp_t= TanDelta(ind)
    normTan1=norm(Ovlp_t)
    Ovlp_t= Ovlp_t-(Ovlp_t .dot. Norm_v)*Norm_v
    normTan2=norm(Ovlp_t)
    if(normTan2>1.0E-10_RK) then
      Ovlp_t=(normTan1/normTan2)*Ovlp_t
    else
      Ovlp_t=zero_r3
    endif
    Ovlp_t= (Vij_t*DEM_opt%dt)+Ovlp_t
   
    ! Computing the normal and tangential contact forces
    select case(DEM_opt%CF_Type)
    case(DEM_LSD)
      Vel_in=max(VelRel_Init(ind),1.0E-5_RK); Vel_in=Vel_in**0.2_RK
      k_n= -Prop_ij%StiffnessCoe_n*Vel_in*Vel_in
      d_n= -Prop_ij%DampingCoe_n*Vel_in
      fn =  k_n*ovrlp +d_n* vrn
      k_t= -Prop_ij%StiffnessCoe_t*Vel_in*Vel_in
      d_t= -Prop_ij%DampingCoe_t*Vel_in
      Ftij= k_t*Ovlp_t +d_t*Vij_t
    case(DEM_nLin)
      k_n= -Prop_ij%StiffnessCoe_n
      d_n= -Prop_ij%DampingCoe_n
      fn =  k_n*ovrlp**1.5_RK +d_n*(ovrlp**0.25_RK)*vrn  ! 2.62, P39
      k_t= -Prop_ij%StiffnessCoe_t*sqrt(ovrlp)
      d_t=  0.0_RK  ! No equation is considered for tangential damping yet
      Ftij= k_t*Ovlp_t +d_t*Vij_t                        ! 2.72, P44
    end select
    Fnij= fn*Norm_v

    ! Coulomb's friction law
    ft = norm(Ftij)
    ft_fric = Prop_ij%FrictionCoe_s*abs(fn)
    if(ft>ft_fric) then
      ft_fric = Prop_ij%FrictionCoe_k* abs(fn)
      Ftij =(ft_fric/ft)*Ftij
      Ovlp_t =Ftij/k_t           ! Caution Here !
    endif
    
    ! Computing rolling resistance torque acting on spheres
    !W_hat = Rvei-Rvej
    !w_hat_mag= norm(W_hat) ! 2.123, p56
    !if(w_hat_mag > 1.0E-10_RK) then
    !  W_hat= W_hat/w_hat_mag
    !else
    !  W_hat= zero_r3
    !endif
    !if(DEM_opt%CT_Model == CTM_ConstantTorque) then 
    !  Mrij=(-Prop_ij%FrictionCoe_Roll*abs(fn)*Prop_ij%RadEff)*W_hat  ! 2.122, p56
    !else
    !  Mrij=(-Prop_ij%FrictionCoe_Roll*abs(fn)*Prop_ij%RadEff*norm(Vel_w))*W_hat  ! 2.124, p57
    !endif

    ! Setting the updated contact info pair into the contact list 
    TanDelta(ind) = Ovlp_t

    ! Updating the contact force and torques of particles i and j
    Moment= Norm_v .cross. Ftij
#ifdef ContactForce_PP
    GPrtcl_cntctForce(pid)= GPrtcl_cntctForce(pid) +(Fnij+Ftij)
    GPrtcl_torque(pid)= GPrtcl_torque(pid) +ri*Moment!+Mrij
    GPrtcl_cntctForce(pjd)= GPrtcl_cntctForce(pjd) -(Fnij+Ftij)
    GPrtcl_torque(pjd)= GPrtcl_torque(pjd) +rj*Moment!-Mrij
#endif
#ifdef ContactForce_PPG
    GPrtcl_cntctForce(pid)= GPrtcl_cntctForce(pid) +(Fnij+Ftij)
    GPrtcl_torque(pid)= GPrtcl_torque(pid) +ri*Moment!+Mrij
#endif
#ifdef ContactForce_PGP
    GPrtcl_cntctForce(pid)= GPrtcl_cntctForce(pid) -(Fnij+Ftij)
    GPrtcl_torque(pid)= GPrtcl_torque(pid) +rj*Moment!-Mrij
#endif
#ifdef ContactForce_PPFix_W
    GPrtcl_cntctForce(pid)= GPrtcl_cntctForce(pid) +(Fnij+Ftij)
    GPrtcl_torque(pid)= GPrtcl_torque(pid) +ri*Moment!+Mrij
#endif
