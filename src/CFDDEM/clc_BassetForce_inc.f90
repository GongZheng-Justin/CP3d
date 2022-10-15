!================================== Basset History force ==================================!
! [1] Casas, G. , Ferrer, A. , & OAte, E. . (2017). Journal of Computational Physics
!       Approximating the basset force by optimizing the method of van hinsberg et al.
! [2] Daitche, A. (2013). Journal of Computational Physics, 93-106.
!       Advection of inertial particles in the presence of the history force.

!******************************************************************
! clc_BassetForce
!******************************************************************
subroutine clc_BassetForce(ux,uy,uz,RatioYp_interp,RatioYc_interp)
  implicit none
  real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
  type(real3),dimension(GPrtcl_list%mlocal),intent(in)::RatioYp_interp,RatioYc_interp

  ! locals
  real(RK)::BassetCoe,mu_coe(mWinBasset+1),sqrt_dt,radius
  integer::i,pid,ws0,wsm,wem,tsm,nlocal,HistoryStage,HistStageFix
  type(real3)::IntegralOld,IntegralNew,f_tail_di,vfmvpTail,BassetForce

  nlocal=GPrtcl_list%nlocal 

  sqrt_dt= sqrt(dt) 
  BassetCoe= six*FluidDensity*sqrt(PI*xnu)/dt
  HistoryStage= GPrtcl_BassetSeq%HistoryStage
  HistStageFix= GPrtcl_BassetSeq%HistStageFix

  ws0= GPrtcl_BassetSeq%iWindowStart;  wsm=ws0-1
  wem= GPrtcl_BassetSeq%iWindowEnd-1
  tsm= GPrtcl_BassetSeq%iTailStart-1

  IF(HistoryStage==0) THEN
    DO pid=1,nlocal
      GPrtcl_BassetData(ws0,pid)= GPrtcl_Vfluid(1,pid)- GPrtcl_linVel(1,pid)
      GPrtcl_BassetData(1,pid)= zero_r3
    ENDDO
  ELSEIF(HistoryStage==1) THEN
    mu_coe(1:2)= clc_muCoeBasset(1)
    DO pid=1,nlocal
      GPrtcl_BassetData(ws0+1,pid)= GPrtcl_BassetData(ws0,pid)
      GPrtcl_BassetData(ws0,  pid)= GPrtcl_Vfluid(1,pid)- GPrtcl_linVel(1,pid)
      
      GPrtcl_BassetData(1,pid)= sqrt_dt*(mu_coe(1)*GPrtcl_BassetData(ws0,pid)+mu_coe(2)*GPrtcl_BassetData(ws0+1,pid))
    ENDDO
  ELSEIF(HistoryStage <= mWinBasset) THEN  ! calculate Basset force directly

    mu_coe(1:HistoryStage+1)= clc_muCoeBasset(HistoryStage)
    DO pid=1,nlocal
      radius= GPrtcl_PosR(pid)%w

      ! shift and update vfmvpWin firstly
      do i=wem,ws0,-1
        GPrtcl_BassetData(i+1,pid)= GPrtcl_BassetData(i,pid)
      enddo
      GPrtcl_BassetData(ws0,pid)= GPrtcl_Vfluid(1,pid)- GPrtcl_linVel(1,pid)

      IntegralNew= zero_r3
      do i=1,HistoryStage+1
        IntegralNew= IntegralNew + mu_coe(i)*GPrtcl_BassetData(i+wsm,pid)
      enddo
      IntegralNew= sqrt_dt*IntegralNew
      IntegralOld= GPrtcl_BassetData(1,pid)
      GPrtcl_BassetData(1,pid)= IntegralNew

      BassetForce= BassetCoe*radius*radius*(IntegralNew-IntegralOld)
      GPrtcl_FpForce(pid)= GPrtcl_FpForce(pid)+ BassetForce
    ENDDO
  ELSE  ! calculate Basset force using window-exponential approximation method

    DO pid=1,nlocal
      radius= GPrtcl_PosR(pid)%w

      ! shift and update vfmvpWin firstly
      do i=wem,ws0,-1
        GPrtcl_BassetData(i+1,pid)= GPrtcl_BassetData(i,pid)
      enddo
      GPrtcl_BassetData(ws0,pid)= GPrtcl_Vfluid(1,pid)- GPrtcl_linVel(1,pid)

      IntegralNew= zero_r3
      do i=1,mWinBasset+1  ! windows part
        IntegralNew= IntegralNew + WindowsCoeBasset(i)*GPrtcl_BassetData(i+wsm,pid)
      enddo
      IntegralNew= sqrt_dt*IntegralNew
      do i=1,mTailBasset   ! tail part
        f_tail_di = tailDiCoeBasset(i,1)*GPrtcl_BassetData(wem+1, pid) + &
                    tailDiCoeBasset(i,2)*GPrtcl_BassetData(wem,   pid) + &
                    tailDiCoeBasset(i,3)*GPrtcl_BassetData(wem-1, pid) + &
                    tailDiCoeBasset(i,4)*GPrtcl_BassetData(wem-2, pid)
        vfmvpTail = f_tail_di + exp_tTailBasset(i)  *GPrtcl_BassetData(tsm+i, pid)
        GPrtcl_BassetData(tsm+i,pid)= vfmvpTail
        IntegralNew= IntegralNew+ aTailBasset(i)* vfmvpTail
      enddo
      IntegralOld= GPrtcl_BassetData(1,pid)
      GPrtcl_BassetData(1,pid)= IntegralNew

      BassetForce= BassetCoe*radius*radius*(IntegralNew-IntegralOld)
      GPrtcl_FpForce(pid)= GPrtcl_FpForce(pid)+ BassetForce
    ENDDO
  ENDIF
  GPrtcl_BassetSeq%HistoryStage= HistoryStage+1
  
  if(.not.is_clc_Basset_fixed)return
  IF(HistStageFix==0) THEN
    DO pid=1,mlocalFix
      GPFix_BassetData(ws0,pid)= GPFix_Vfluid(1,pid)
      GPFix_BassetData(1,pid)= zero_r3
    ENDDO
  ELSEIF(HistStageFix==1) THEN
    mu_coe(1:2)= clc_muCoeBasset(1)
    DO pid=1,mlocalFix
      GPFix_BassetData(ws0+1,pid)= GPFix_BassetData(ws0,pid)
      GPFix_BassetData(ws0,  pid)= GPFix_Vfluid(1,pid)
      GPFix_BassetData(1,pid)    = sqrt_dt*(mu_coe(1)*GPFix_BassetData(ws0,pid)+mu_coe(2)*GPFix_BassetData(ws0+1,pid))
    ENDDO
  ELSEIF(HistStageFix <= mWinBasset) THEN  ! calculate Basset force directly

    mu_coe(1:HistStageFix+1)= clc_muCoeBasset(HistStageFix)
    DO pid=1,mlocalFix
      radius= GPFix_PosR(pid)%w

      ! shift and update vfmvpWin firstly
      do i=wem,ws0,-1
        GPFix_BassetData(i+1,pid)= GPFix_BassetData(i,pid)
      enddo
      GPFix_BassetData(ws0,pid)= GPFix_Vfluid(1,pid)

      IntegralNew= zero_r3
      do i=1,HistStageFix+1
        IntegralNew= IntegralNew + mu_coe(i)*GPFix_BassetData(i+wsm,pid)
      enddo
      IntegralOld= GPFix_BassetData(1,pid)
      IntegralNew= sqrt_dt*IntegralNew
      GPFix_BassetData(1,pid)= IntegralNew
 
      BassetForce= BassetCoe*radius*radius*(IntegralNew-IntegralOld)
      GPFix_FpForce(pid)= GPFix_FpForce(pid)+ BassetForce
      ENDDO
  ELSE  ! calculate Basset force using window-exponential approximation method
    DO pid=1,mlocalFix
      radius= GPFix_PosR(pid)%w

      ! shift and update vfmvpWin firstly
      do i=wem,ws0,-1
        GPFix_BassetData(i+1,pid)= GPFix_BassetData(i,pid)
      enddo
      GPFix_BassetData(ws0,pid)= GPFix_Vfluid(1,pid)

      IntegralNew= zero_r3
      do i=1,mWinBasset+1  ! windows part
        IntegralNew= IntegralNew + WindowsCoeBasset(i)*GPFix_BassetData(i+wsm,pid)
      enddo
      IntegralNew= sqrt_dt*IntegralNew
      do i=1,mTailBasset   ! tail part
        f_tail_di = tailDiCoeBasset(i,1)*GPFix_BassetData(wem+1, pid) + &
                    tailDiCoeBasset(i,2)*GPFix_BassetData(wem,   pid) + &
                    tailDiCoeBasset(i,3)*GPFix_BassetData(wem-1, pid) + &
                    tailDiCoeBasset(i,4)*GPFix_BassetData(wem-2, pid)
        vfmvpTail = f_tail_di + exp_tTailBasset(i)  *GPFix_BassetData(tsm+i, pid)
        GPFix_BassetData(tsm+i,pid)= vfmvpTail
        IntegralNew= IntegralNew+ aTailBasset(i)* vfmvpTail
      enddo
      IntegralOld= GPFix_BassetData(1,pid)
      GPFix_BassetData(1,pid)= IntegralNew

      BassetForce= BassetCoe*radius*radius*(IntegralNew-IntegralOld)
      GPFix_FpForce(pid)= GPFix_FpForce(pid)+ BassetForce
    ENDDO
  ENDIF  
  GPrtcl_BassetSeq%HistStageFix= HistStageFix+1

end subroutine clc_BassetForce

!******************************************************************
! clc_TailCoeBasset
!******************************************************************
subroutine clc_TailCoeBasset()
  implicit none

  ! locals
  integer :: j
  real(RK):: tTailBasset_i,alpha_ti,z_tail,tailDiCoeBasset_multi,rmwin,zp1,zp2,zp3,zp4,zp5

  ! ------------------------- aTailBasset, tTailBasset  -------------------------
  IF(BassetTailType==1) THEN  ! G.Casas et al./ Journal of Computational Physics 352(2018), p168
    if(mTailBasset== 1) then
      aTailBasset=(/1.046347992_RK /)
      tTailBasset=(/1.581186674_RK /)
    endif 
    if(mTailBasset== 2) then
      aTailBasset=(/0.566192817_RK,      0.864298391_RK /)
      tTailBasset=(/0.717656182_RK,      8.925153279_RK /)
    endif 
    if(mTailBasset== 3) then 
      aTailBasset=(/0.440072204_RK,      0.538287204_RK,     0.807797346_RK /)
      tTailBasset=(/0.482318894_RK,      3.324763126_RK,    38.928376132_RK /)
    endif 
    if(mTailBasset== 4) then 
      aTailBasset=(/0.374397988_RK,      0.421322343_RK,     0.517872275_RK,      0.761539469_RK /)
      tTailBasset=(/0.365083559_RK,      1.820334739_RK,    11.809488351_RK,    127.109159354_RK /)  
    endif           
    if(mTailBasset== 5) then 
      aTailBasset=(/0.3450551877_RK,    0.3762685526_RK,    0.4383511621_RK,     0.5502868981_RK,     0.7701813938_RK /)
      tTailBasset=(/0.3227320427_RK,    1.4017593843_RK,    7.3543952717_RK,    52.9058339347_RK,   699.4337431732_RK /)
    endif 
    if(mTailBasset== 6) then 
      aTailBasset=(/0.3227460255_RK,    0.3446901326_RK,    0.3924441164_RK,  &
                     0.471576099_RK,    0.5990063177_RK,    0.7695849793_RK  /)
      tTailBasset=(/0.2894856389_RK,    1.1312690586_RK,    5.1207861657_RK,  &
                   29.6345412934_RK,    256.64908268_RK, 4254.1241751139_RK  /)
    endif 
    if(mTailBasset== 7) then 
      aTailBasset=(/ 0.2931405176_RK,   0.3053190176_RK,    0.3394616674_RK,     0.3924532926_RK,  &
                     0.4794140412_RK,   0.5546383969_RK,    0.6207864425_RK                       /)
      tTailBasset=(/ 0.2413624327_RK,   0.8199848671_RK,    3.0838532791_RK,    13.8047974118_RK,  &
                    80.9779742728_RK, 696.8320792921_RK, 6133.2449027098_RK                       /) 
    endif 
    if(mTailBasset== 8) then 
      aTailBasset=(/ 0.2718360249_RK,   0.2685924185_RK,    0.2871214552_RK,     0.3249589764_RK,  &
                     0.3805886345_RK,   0.4469592071_RK,    0.5474439544_RK,     0.7637048975_RK  /)
      tTailBasset=(/ 0.2192620346_RK,    0.662026818_RK,    2.0706383247_RK,     7.2825402363_RK,  &
                    31.0062809826_RK, 169.6857783353_RK,  1226.001409491_RK, 17271.9375778519_RK  /) 
    endif   
    if(mTailBasset== 9) then 
      aTailBasset=(/ 0.2570818336_RK,   0.2610118588_RK,    0.2799238451_RK,     0.3051985477_RK,     0.3418149337_RK, &
                     0.3892337642_RK,   0.4655655296_RK,    0.6107696402_RK,      0.784623916_RK                      /)
      tTailBasset=(/ 0.1878604572_RK,   0.5420260992_RK,    1.6534881587_RK,     5.5204876302_RK,    20.8847203692_RK, &
                    93.9005719593_RK, 532.1532341216_RK, 4683.3937018005_RK, 93277.7129340798_RK                      /)
    endif          
    if(mTailBasset==10) then
      aTailBasset=(/ 0.2520642358_RK,   0.2549130660_RK,    0.2638832071_RK,     0.2666445191_RK,     0.2806268115_RK,   &
                     0.3449146080_RK,   0.4566204962_RK,    0.5663046247_RK,     0.6253574036_RK,     0.6932526975_RK   /)
      tTailBasset=(/ 0.1878604572_RK,   0.5306382498_RK,    1.5524873935_RK,     4.6517443725_RK,    14.2413555446_RK,   &
                    50.7413819742_RK, 263.7561507819_RK,  2146.211201895_RK,  26744.590748687_RK, 348322.670028861_RK   /)           
    endif
  ENDIF

  IF(BassetTailType==2) THEN  ! G.Casas et al./ Journal of Computational Physics 352(2018), p169
    if (mTailBasset==1) then
      aTailBasset= (/0.9384724434_RK /)
      tTailBasset= (/1.4300340551_RK /)
    endif
    if (mTailBasset==2) then
      aTailBasset= (/0.5470597552_RK,   0.8449767491_RK /) 
      tTailBasset= (/0.6666835275_RK,   8.3424872407_RK /)
    endif
    if (mTailBasset==3) then
      aTailBasset=(/  0.430797005_RK,   0.5319402016_RK,    0.8046471493_RK /)
      tTailBasset=(/ 0.4521461414_RK,   3.0597097311_RK,    36.769402335_RK /)
    endif
    if (mTailBasset==4) then
      aTailBasset=(/ 0.3714051613_RK,   0.4221306386_RK,   0.52488276380_RK,     0.7814317902_RK/)
      tTailBasset=(/ 0.3505056162_RK,   1.7525741335_RK,   11.6528756138_RK,   136.8864124598_RK/)
    endif
    if (mTailBasset==5) then
      aTailBasset=(/ 0.3335736291_RK,   0.3629331173_RK,    0.4197252519_RK,      0.520201698_RK,     0.7661038702_RK/)
      tTailBasset=(/ 0.2904610289_RK,   1.2036915740_RK,    5.9370324806_RK,    39.1450598115_RK,   452.8226228869_RK/)
    endif
    if (mTailBasset==6) then
      aTailBasset=(/ 0.3065928563_RK,   0.3243480187_RK,    0.3615932545_RK, &
                      0.418122689_RK,   0.5168085735_RK,    0.7551149413_RK /)
      tTailBasset=(/ 0.2504309713_RK,   0.9103056758_RK,    3.7204976994_RK, &
                    18.2727422761_RK,  119.760302387_RK, 1369.9016377844_RK /)
    endif
    if (mTailBasset==7) then
      aTailBasset=(/ 0.2869667584_RK,   0.2977774360_RK,    0.3249218804_RK,     0.3631687423_RK, &
                      0.420482447_RK,   0.5207711634_RK,    0.7554318595_RK                      /) 
      tTailBasset=(/ 0.2229567355_RK,   0.7379950193_RK,    2.6583099103_RK,    10.9237321521_RK, &
                     54.149026921_RK, 360.6375769122_RK, 4254.1243411105_RK                      /)
    endif
    if (mTailBasset==8) then
      aTailBasset=(/ 0.2695926115_RK,   0.2751628513_RK,    0.2954155007_RK,     0.3228159111_RK, &
                     0.3602702642_RK,   0.4159293673_RK,    0.5121568839_RK,     0.7402280446_RK /) 
      tTailBasset=(/ 0.1998084724_RK,   0.6094217589_RK,    1.9746292865_RK,     7.0527307887_RK, &
                    28.6942745173_RK, 140.1890961279_RK,  911.2555045811_RK, 10263.3419763251_RK /) 
    endif
    if (mTailBasset==9) then
      aTailBasset=(/ 0.2560766303_RK,   0.2580812359_RK,    0.2739535236_RK,     0.2950977014_RK,     0.3224941082_RK, &
                     0.3598005017_RK,   0.4151331109_RK,    0.5104760265_RK,     0.7348997012_RK                      /)
      tTailBasset=(/ 0.1826244466_RK,   0.5231166336_RK,    1.5665809982_RK,     5.0639463936_RK,    18.0664363914_RK, &
                    73.4054449448_RK, 357.9494752882_RK, 2319.7684648904_RK, 25980.6116922192_RK                      /) 
    endif
    if (mTailBasset==10) then
      aTailBasset=(/ 0.2467020831_RK,   0.2464749444_RK,    0.2597178682_RK,     0.2773405882_RK,     0.2995010019_RK, &
                     0.3282822047_RK,   0.3678821811_RK,    0.4276240337_RK,     0.5335800139_RK,     0.7652665389_RK /) 
      tTailBasset=(/ 0.1711374102_RK,   0.4695384557_RK,    1.3336047234_RK,     4.0387298490_RK,    13.2686834339_RK, &
                    48.3505553197_RK, 202.2013044128_RK, 1029.0899279619_RK,  7177.8752909387_RK, 93277.7373733731_RK /) 
    endif
  ENDIF

  ! ------------------------- tailDiCoeBasset -------------------------
  rmwin= real(mWinBasset,RK)
  do j = 1, mTailBasset

    exp_tTailBasset(j)= exp(-0.5_RK/(tTailBasset(j)*real(mWinBasset,RK)))
    tTailBasset_i     = tTailBasset(j)*real(mWinBasset,RK)*dt

    alpha_ti  = sqrt( exp(1.0_RK)/tTailBasset_i )
    z_tail    = 0.5_RK/(tTailBasset(j)*rmWin)
    tailDiCoeBasset_multi  = dt * alpha_ti * exp(-0.5_RK/tTailBasset(j))
    zp1=z_tail; zp2=zp1*zp1; zp3=zp1*zp2; zp4=zp1*zp3; zp5=zp1*zp4

    if(BassetAccuracy==1)then
      if(abs(z_tail) > 0.05_RK) then
        tailDiCoeBasset(j,1) = ( 1.0_RK-(zp1+1.0_RK)*exp(-zp1) )/zp2
        tailDiCoeBasset(j,2) = ( exp(-zp1)+zp1 - 1.0_RK )/zp2
      else
        tailDiCoeBasset(j,1) = 0.5_RK - 1.0_RK/3.0_RK*zp1 + 1.0_RK/8.0_RK*zp2 -  1.0_RK/30.0_RK*zp3 + 1.0_RK/144.0_RK*zp4 -  1.0_RK/840.0_RK*zp5
        tailDiCoeBasset(j,2) = 0.5_RK - 1.0_RK/6.0_RK*zp1 + 1.0_RK/24.0_RK*zp2- 1.0_RK/120.0_RK*zp3 + 1.0_RK/720.0_RK*zp4 - 1.0_RK/5040.0_RK*zp5
      endif
      tailDiCoeBasset(j,3:4) = 0.0_RK
    endif

    if(BassetAccuracy==2)then
      if(abs(z_tail) > 0.05_RK) then
        tailDiCoeBasset(j,1) = (-(2.0_RK*zp2+3.0_RK*zp1+2.0_RK)*exp(-zp1)+zp1+2.0_RK )/zp3*0.5_RK
        tailDiCoeBasset(j,2) = (2.0_RK*(zp1+1.0_RK)*exp(-zp1)+zp2-2.0_RK )/zp3
        tailDiCoeBasset(j,3) = (-(zp1+2.0_RK)*exp(-zp1)-zp1+2.0_RK )/zp3*0.5_RK
      else
        tailDiCoeBasset(j,1) = 5.0_RK/12.0_RK - 7.0_RK/24.0_RK*zp1 + 9.0_RK/80.0_RK*zp2 - 11.0_RK/360.0_RK*zp3 + 13.0_RK/2016.0_RK*zp4 -   1.0_RK/896.0_RK*zp5
        tailDiCoeBasset(j,2) =  2.0_RK/3.0_RK - 1.0_RK/ 4.0_RK*zp1 + 1.0_RK/15.0_RK*zp2 -   1.0_RK/72.0_RK*zp3 +   1.0_RK/420.0_RK*zp4 -  1.0_RK/2880.0_RK*zp5
        tailDiCoeBasset(j,3) =-1.0_RK/12.0_RK + 1.0_RK/24.0_RK*zp1 - 1.0_RK/80.0_RK*zp2 +  1.0_RK/360.0_RK*zp3 -  1.0_RK/2016.0_RK*zp4 + 1.0_RK/13440.0_RK*zp5
      endif
      tailDiCoeBasset(j,4) = 0.0_RK
    endif

    if(BassetAccuracy==3) then
      if(abs(z_tail) > 0.05_RK) then 
        tailDiCoeBasset(j,1)= (-(6.0_RK*zp3 + 11.0_RK*zp2 + 12.0_RK*zp1 + 6.0_RK)*exp(-zp1) + 2.0_RK*zp2 + 6.0_RK*zp1 + 6.0_RK)/(6.0_RK * zp4)
        tailDiCoeBasset(j,2)= ( (6.0_RK*zp2 + 10.0_RK*zp1 + 6.0_RK)*exp(-zp1) + 2.0_RK*zp3 + zp2 -4.0_RK*zp1 - 6.0_RK)/(2.0_RK * zp4)
        tailDiCoeBasset(j,3)= (-(3.0_RK*zp2 + 8.0_RK*zp1  + 6.0_RK)*exp(-zp1) - 2.0_RK*zp2 + 2.0_RK*zp1 + 6.0_RK)/(2.0_RK * zp4)
        tailDiCoeBasset(j,4)= ( (2.0_RK*zp2 + 6.0_RK*zp1  + 6.0_RK)*exp(-zp1) + zp2 - 6.0_RK)/(6.0_RK * zp4)
      else
        tailDiCoeBasset(j,1)= 9.0_RK/24.0_RK -97.0_RK/360.0_RK*zp1 +19.0_RK/180.0_RK*zp2 -73.0_RK/2520.0_RK*zp3 +149.0_RK/24192.0_RK*zp4 - 389.0_RK/362880.0_RK*zp5
        tailDiCoeBasset(j,2)=19.0_RK/24.0_RK -19.0_RK/ 60.0_RK*zp1 + 7.0_RK/ 80.0_RK*zp2 -47.0_RK/2520.0_RK*zp3 +131.0_RK/40320.0_RK*zp4 -   29.0_RK/60480.0_RK*zp5
        tailDiCoeBasset(j,3)=-5.0_RK/24.0_RK +13.0_RK/120.0_RK*zp1 - 1.0_RK/ 30.0_RK*zp2 +19.0_RK/2520.0_RK*zp3 -  11.0_RK/8064.0_RK*zp4 +    5.0_RK/24192.0_RK*zp5
        tailDiCoeBasset(j,4)= 1.0_RK/24.0_RK - 1.0_RK/ 45.0_RK*zp1 + 1.0_RK/144.0_RK*zp2 -  1.0_RK/630.0_RK*zp3 +   1.0_RK/3456.0_RK*zp4 -    1.0_RK/22680.0_RK*zp5
      endif
    endif
    tailDiCoeBasset(j,1:4) = tailDiCoeBasset(j,1:4) * tailDiCoeBasset_multi
  enddo
end subroutine clc_TailCoeBasset

!******************************************************************
! clc_muCoeBasset
!******************************************************************
function clc_muCoeBasset(nc) result(muCoeBasset)
  implicit none
  integer,intent(in) :: nc
  real(RK),dimension(nc+1)::muCoeBasset

  ! locals
  integer :: j
  real(RK):: sqrt2,sqrt3,rn,rn0,rn1,rn2,rn3,rn4,rn5,rj,rjp0,rjp1,rjp2,rjm1,rjm2

  ! Reference : A.Daitche /Journal of Computational Physics 254(2013)93-106, p97-98
  sqrt2= sqrt(2.0_RK); sqrt3= sqrt(3.0_RK)

  rn = real(nc,RK);   rn0 = sqrt(rn)
  if(nc>0) rn1 =sqrt(abs(rn-1.0_RK))
  if(nc>1) rn2 =sqrt(abs(rn-2.0_RK))
  if(nc>2) rn3 =sqrt(abs(rn-3.0_RK)) 
  if(nc>3) rn4 =sqrt(abs(rn-4.0_RK))
  if(nc>4) rn5 =sqrt(abs(rn-5.0_RK))

  IF(BassetAccuracy==1) THEN ! First order
    muCoeBasset(1)=1.0_RK
    do j=1,nc-1
      rj  = real(j,RK)
      rjp0= sqrt(rj); rjp1 = sqrt(rj+1.0_RK); rjm1 = sqrt(abs(rj-1.0_RK))
      if(j<=100) then
        muCoeBasset(j+1) = rjm1**3 +rjp1**3 -2.0_RK*rjp0**3
      else
        muCoeBasset(j+1) = 1.0_RK/rjp0 *(3.0_RK/4.0_RK +3.0_RK/64.0_RK/rj**2 +7.0_RK/512.0_RK/rj**4)
      endif
    enddo
    if(nc<=100) then
      muCoeBasset(nc+1)= rn1**3 -rn0**3 +6.0_RK/4.0_RK*rn0
    else
      muCoeBasset(nc+1)= 1.0_RK/rn0*(3.0_RK/8.0_RK +1.0_RK/16.0_RK/rn +3.0_RK/128.0_RK/rn**2 +3.0_RK/256.0_RK/rn**3 +7.0_RK/1024.0_RK/rn**4 +9.0_RK/2048.0_RK/rn**5)
    endif
    muCoeBasset = muCoeBasset*4.0_RK/3.0_RK
  ENDIF

  IF(BassetAccuracy==2) THEN ! Second order
    if(nc==1) muCoeBasset=(/4.0_RK/3.0_RK, 2.0_RK/3.0_RK /)
    if(nc==2) muCoeBasset=(/12.0_RK, 16.0_RK, 2.0_RK  /) * sqrt2/15.0_RK
    if(nc==3) muCoeBasset=(/4.0_RK*sqrt2,14.0_RK*sqrt3-12.0_RK*sqrt2,-8.0_RK*sqrt3+12.0_RK*sqrt2,4.0_RK*sqrt3-4.0_RK*sqrt2/)/5.0_RK
    if(nc>3) then
      muCoeBasset(1:3)=(/4.0_RK/5.0_RK*sqrt2, 14.0_RK/5.0_RK*sqrt3-12.0_RK/5.0_RK*sqrt2, 176.0_RK/15.0_RK-42.0_RK/5.0_RK*sqrt3+12.0_RK/5.0_RK*sqrt2/)
      do j=3,nc-2
        rj  = real(j,RK);       rjp0 = sqrt(rj); 
        rjp1 = sqrt(rj+1.0_RK); rjp2 = sqrt(rj+2.0_RK)
        rjm1 = sqrt(rj-1.0_RK); rjm2 = sqrt(rj-2.0_RK)
        if(j<=100) then
          muCoeBasset(j+1)= 8.0_RK/15.0_RK *( rjp2**5 -3.0_RK*rjp1**5 +3.0_RK*rjp0**5 -rjm1**5) &
                           + 2.0_RK/3.0_RK *(-rjp2**3 +3.0_RK*rjp1**3 -3.0_RK*rjp0**3 +rjm1**3)
        else
          muCoeBasset(j+1)= 1.0_RK/rjp0 *(1.0_RK +5.0_RK/64.0_RK/rj**3 -7.0_RK/64.0_RK/rj**4 +189.0_RK/1024.0_RK/rj**5)
        endif

      enddo
      if(nc<=100) then
        muCoeBasset(nc)  = 8.0_RK/15.0_RK *(-2.0_RK*rn0**5  +3.0_RK*rn1**5 -rn2**5) &
                         + 2.0_RK/ 3.0_RK *( 4.0_RK*rn0**3  -3.0_RK*rn1**3 +rn2**3)  
        muCoeBasset(nc+1)= 8.0_RK/15.0_RK *(rn0**5 -rn1**5) +2.0_RK/3.0_RK *(-3.0_RK*rn0**3 + rn1**3)+2.0_RK*rn0
      else
        muCoeBasset(nc)  = 1.0_RK/rn0*( 13.0_RK/12.0_RK +23.0_RK/48.0_RK/rn +123.0_RK/320.0_RK/rn**2 +37.0_RK/96.0_RK/rn**3 +677.0_RK/1536.0_RK/rn**4 + 2259.0_RK/4096.0_RK/rn**5)
        muCoeBasset(nc+1)= 1.0_RK/rn0*( 5.0_RK/12.0_RK  + 1.0_RK/16.0_RK/rn +  7.0_RK/320.0_RK/rn**2 + 1.0_RK/96.0_RK/rn**3 +  3.0_RK/ 512.0_RK/rn**4 +       15.0_RK/4096.0_RK/rn**5)    
      endif
    endif
  ENDIF

  IF(BassetAccuracy==3) THEN ! Third order
    if(nc==1) muCoeBasset=(/ 4.0_RK/3.0_RK, 2.0_RK/3.0_RK /)
    if(nc==2) muCoeBasset=(/ 12.0_RK, 16.0_RK, 2.0_RK  /) *  sqrt2/15.0_RK
    if(nc==3) muCoeBasset=(/ 68.0_RK, 90.0_RK, 36.0_RK, 16.0_RK /)   * sqrt3/105.0_RK 
    if(nc==4) then
      muCoeBasset(1) =                      244.0_RK/315.0_RK*sqrt2
      muCoeBasset(2) = 1888.0_RK/315.0_RK - 976.0_RK/315.0_RK*sqrt2
      muCoeBasset(3) = -656.0_RK/105.0_RK + 488.0_RK/105.0_RK*sqrt2
      muCoeBasset(4) =  544.0_RK/105.0_RK - 976.0_RK/315.0_RK*sqrt2
      muCoeBasset(5) = -292.0_RK/315.0_RK + 244.0_RK/315.0_RK*sqrt2
    endif
    if(nc==5) then
      muCoeBasset(1) =  244.0_RK/315.0_RK*sqrt2
      muCoeBasset(2) =  362.0_RK/105.0_RK*sqrt3        -  976.0_RK/315.0_RK*sqrt2
      muCoeBasset(3) =  500.0_RK/63.0_RK *sqrt(5.0_RK) - 1448.0_RK/105.0_RK*sqrt3 + 488.0_RK/105.0_RK*sqrt2
      muCoeBasset(4) = -290.0_RK/21.0_RK *sqrt(5.0_RK) +   724.0_RK/35.0_RK*sqrt3 - 976.0_RK/315.0_RK*sqrt2
      muCoeBasset(5) =  220.0_RK/21.0_RK *sqrt(5.0_RK) - 1448.0_RK/105.0_RK*sqrt3 + 244.0_RK/315.0_RK*sqrt2
      muCoeBasset(6) = -164.0_RK/63.0_RK *sqrt(5.0_RK) +  362.0_RK/105.0_RK*sqrt3
    endif
    if(nc==6) then
      muCoeBasset(1) = 244.0_RK/315.0_RK*sqrt2
      muCoeBasset(2) = 362.0_RK/105.0_RK*sqrt3 - 976.0_RK/315.0_RK*sqrt2
      muCoeBasset(3) =                                    5584.0_RK/315.0_RK - 1448.0_RK/105.0_RK*sqrt3 +  488.0_RK/105.0_RK*sqrt2
      muCoeBasset(4) =   344.0_RK/21.0_RK*sqrt(6.0_RK) - 22336.0_RK/315.0_RK +   724.0_RK/35.0_RK*sqrt3 -  976.0_RK/315.0_RK*sqrt2
      muCoeBasset(5) = -1188.0_RK/35.0_RK*sqrt(6.0_RK) + 11168.0_RK/105.0_RK - 1448.0_RK/105.0_RK*sqrt3 +  244.0_RK/315.0_RK*sqrt2
      muCoeBasset(6) =   936.0_RK/35.0_RK*sqrt(6.0_RK) - 22336.0_RK/315.0_RK +  362.0_RK/105.0_RK*sqrt3
      muCoeBasset(7) = -754.0_RK/105.0_RK*sqrt(6.0_RK) +  5584.0_RK/315.0_RK
    endif
    if(nc>6) then
      muCoeBasset(1) = 244.0_RK/315.0_RK*sqrt2
      muCoeBasset(2) = 362.0_RK/105.0_RK*sqrt3 - 976.0_RK/315.0_RK*sqrt2
      muCoeBasset(3) = 5584.0_RK/315.0_RK - 1448.0_RK/105.0_RK*sqrt3 + 488.0_RK/105.0_RK *sqrt2
      muCoeBasset(4) = 1130.0_RK/63.0_RK*sqrt(5.0_RK) - 22336.0_RK/315.0_RK + 724.0_RK/35.0_RK*sqrt3 - 976.0_RK/315.0_RK*sqrt2
      do j = 4,nc-4
        rj  = real(j,RK)
        rjp0 = sqrt(rj); 
        rjp1 = sqrt(rj+1.0_RK); rjp2 = sqrt(rj+2.0_RK)
        rjm1 = sqrt(rj-1.0_RK); rjm2 = sqrt(rj-2.0_RK)
        if(j<=35) then
          muCoeBasset(j+1)= 16.0_RK/105.0_RK * ( rjp2**7 + rjm2**7 - 4.0_RK*rjp1**7 - 4.0_RK*rjm1**7 + 6.0_RK*rjp0**7) + &
                               2.0_RK/9.0_RK * (-rjp2**3 - rjm2**3 + 4.0_RK*rjp1**3 + 4.0_RK*rjm1**3 - 6.0_RK*rjp0**3)
        else
          muCoeBasset(j+1)= 1.0_RK/rjp0*(1.0_RK- 77.0_RK/768.0_RK/rj**4 -253.0_RK/1024.0_RK/rj**6 -19877.0_RK/32768.0_RK/rj**8)         
        endif
      enddo
      if(nc<=35) then
        muCoeBasset(nc-2)= 16.0_RK/105.0_RK*((rn0**7) -4.0_RK*(rn2**7) + 6.0_RK*(rn3**7) - 4.0_RK*(rn4**7)+ (rn5**7) ) - 8.0_RK/15.0_RK*(rn0**5) + &
                            4.0_RK/9.0_RK*(rn0**3) + 8.0_RK/9.0_RK*(rn2**3) - 4.0_RK/3.0_RK*(rn3**3) + 8.0_RK/9.0_RK*(rn4**3) -2.0_RK/9.0_RK*(rn5**3)
        muCoeBasset(nc-1)= 16.0_RK/105.0_RK*((rn4**7) -4.0_RK*(rn3**7) + 6.0_RK*(rn2**7) - 3.0_RK*(rn0**7)) + 32.0_RK/15.0_RK*(rn0**5) - &
                            2.0_RK*(rn0**3) - 4.0_RK/3.0_RK*(rn2**3) + 8.0_RK/9.0_RK*(rn3**3) -2.0_RK/9.0_RK*(rn4**3)
        muCoeBasset(nc  )= 16.0_RK/105.0_RK * (3.0_RK*(rn0**7) - 4.0_RK*(rn2**7) + (rn3**7)) - &
                            8.0_RK/3.0_RK*(rn0**5) + 4.0_RK*(rn0**3) + 8.0_RK/9.0_RK*(rn2**3) - 2.0_RK/9.0_RK*(rn3**3)
        muCoeBasset(nc+1)= 16.0_RK/105.0_RK * ((rn2**7) - rn0**7)  + 16.0_RK/15.0_RK * (rn0**5) - & 
                           22.0_RK/9.0_RK* (rn0**3) - 2.0_RK/9.0_RK* (rn2**3) + 2.0_RK*rn0
      else
        muCoeBasset(nc-2)= 1.0_RK/rn0*(  25.0_RK/24.0_RK +     1087.0_RK/720.0_RK/rn    +      811.0_RK/240.0_RK/rn**2 +        3781.0_RK/448.0_RK/rn**3 &
                                                         + 406417.0_RK/18432.0_RK/rn**4 + 718211.0_RK/12288.0_RK/rn**5 + 18805457.0_RK/122880.0_RK/rn**6 )
        muCoeBasset(nc-1)= 1.0_RK/rn0*(  20.0_RK/24.0_RK +      173.0_RK/180.0_RK/rn    +        89.0_RK/60.0_RK/rn**2 +         279.0_RK/112.0_RK/rn**3 &
                                                         +   19673.0_RK/4608.0_RK/rn**4 +   21409.0_RK/3072.0_RK/rn**5 +    283393.0_RK/30720.0_RK/rn**6 )
        muCoeBasset(nc  )= 1.0_RK/rn0*(  31.0_RK/24.0_RK +       79.0_RK/144.0_RK/rn    +       97.0_RK/240.0_RK/rn**2 +         145.0_RK/448.0_RK/rn**3 &
                                                         +   3367.0_RK/18432.0_RK/rn**4 -   2465.0_RK/12288.0_RK/rn**5 -   154561.0_RK/122880.0_RK/rn**6 )
        muCoeBasset(nc+1)= 1.0_RK/rn0*(8.0_RK/24.0_RK    +         1.0_RK/45.0_RK/rn    -         1.0_RK/60.0_RK/rn**2 -            1.0_RK/28.0_RK/rn**3 &
                                                         -         1.0_RK/18.0_RK/rn**4 -         1.0_RK/12.0_RK/rn**5 -         121.0_RK/960.0_RK/rn**6 )
      endif
    endif
  ENDIF
end function clc_muCoeBasset
